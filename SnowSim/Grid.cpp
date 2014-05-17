#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* object){
	obj = object;
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	int len = size.product();
	nodes = new GridNode[len];
}
Grid::Grid(const Grid& orig){}
Grid::~Grid(){
	delete[] nodes;
}

//Maps mass and velocity to the grid
void Grid::initializeMass(){
	//Reset the grid
	//If the grid is sparsely filled, it may be better to reset individual nodes
	memset(nodes, 0, sizeof(GridNode)*size.product());
	
	//Map particle data to grid
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//Particle position to grid coordinates
		//This will give errors if the particle is outside the grid bounds
		p.grid_position = (p.position - origin)/cellsize;
		float ox = p.grid_position[0],
			oy = p.grid_position[1];
		
		//Shape function gives a blending radius of two;
		//so we do computations within a 2x2 square for each particle
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			//X-dimension interpolation
			float x_pos = fabs(x-ox),
				wx = Grid::bspline(x_pos),
				dx = Grid::bsplineSlope(x_pos);
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				//Y-dimension interpolation
				float y_pos = fabs(y-oy),
					wy = Grid::bspline(y_pos),
					dy = Grid::bsplineSlope(y_pos);
				
				//Final weight is dyadic product of weights in each dimension
				float weight = wx*wy;
				p.weights[idx] = weight;
				
				//Weight gradient is a vector of partial derivatives
				p.weight_gradient[idx].setPosition(dx*wy, wx*dy);
				
				//Interpolate mass
				nodes[(int) (x*size[0]+y)].mass += weight*p.mass;
			}
		}
	}
}

void Grid::initializeVelocities(){
	//We interpolate velocity after mass, to conserve momentum
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Interpolate velocity
					int n = (int) (x*size[0]+y);
					nodes[n].velocity += p.velocity * w * (p.mass/nodes[n].mass);
				}
			}
		}
	}
	
	//Particles need the velocity gradient to compute stress forces
	//This is only temporary; the actual "resolved" velocity gradient
	//is computed in updateVelocities()
	for (int i=0; i<obj->size; i++){
		//Reset gradient to zero
		Particle& p = obj->particles[i];
		Matrix2f& grad = p.velocity_gradient;
		grad.setData(0.0);
		
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					grad += nodes[(int) (x*size[0]+y)].velocity.trans_product(p.weight_gradient[idx]);
				}
			}
		}
	}
}
//Maps volume from the grid to particles
//This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes() const{
	//Estimate each particles volume (for force calculations)
	//TODO: I multiply by 1.435 because this calculation is underestimating the actual volume; not sure why ????
	float node_volume = (1.435*cellsize).product();
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		//First compute particle density
		p.density = 0;
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Node density is trivial
					float node_density = nodes[(int) (x*size[0]+y)].mass / node_volume;
					//Weighted sum of nodes
					p.density += w * node_density;
				}
			}
		}
		//Volume for each particle can be found from density
		p.volume = p.mass / p.density;
	}
}
//Calculate next timestep velocities for use in explicit/implicit integration
void Grid::calculateVelocities(const Vector2f& gravity){
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//Solve for grid internal forces
		Matrix2f force = p.stressForce();
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Weight the force onto nodes
					int n = (int) (x*size[0]+y);
					nodes[n].force += force*p.weight_gradient[idx];
					nodes[n].has_force = true;
				}
			}
		}
	}
	//Now we have all grid forces, compute velocities (euler integration)
	for (int i=0, l=size.product(); i<l; i++){
		GridNode &node = nodes[i];
		//Check to see if this node needs to be computed
		if (node.has_force){
			node.velocity_new = node.velocity + TIMESTEP*(node.force/node.mass+gravity);
			node.force.setPosition(0);
			node.has_force = false;
		}
	}
}
//Map grid velocities back to particles
void Grid::updateVelocities() const{
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//We calculate PIC and FLIP velocities separately
		Vector2f pic, flip = p.velocity;
		//Also keep track of velocity gradient
		Matrix2f& grad = p.velocity_gradient;
		grad.setData(0.0);
		
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					GridNode &node = nodes[(int) (x*size[0]+y)];
					//Particle in cell
					pic += node.velocity_new*w;
					//Fluid implicit particle
					flip += (node.velocity_new - node.velocity)*w;
					//Velocity gradient
					grad += node.velocity_new.trans_product(p.weight_gradient[idx]);
				}
			}
		}
		//Final velocity is a linear combination of PIC and FLIP components
		p.velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
		std::cout << "Gradient = " << std::endl;
		grad.print();
	}
}
