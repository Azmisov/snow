#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* object){
	obj = object;
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	int len = size.product();
	nodes = new GridNode[len];
	//TODO: this is underestimating? perhaps we can scale by 1.435
	node_volume = cellsize.product();
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
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			//Y-dimension interpolation
			float y_pos = fabs(y-oy),
				wy = Grid::bspline(y_pos),
				dy = Grid::bsplineSlope(y_pos);
			
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				//X-dimension interpolation
				float x_pos = fabs(x-ox),
					wx = Grid::bspline(x_pos),
					dx = Grid::bsplineSlope(x_pos);
				
				//Final weight is dyadic product of weights in each dimension
				float weight = wx*wy;
				p.weights[idx] = weight;
				
				//Weight gradient is a vector of partial derivatives
				p.weight_gradient[idx].setPosition(dx*wy, wx*dy);
				
				//Interpolate mass
				nodes[(int) (y*size[0]+x)].mass += weight*p.mass;
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
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Interpolate velocity
					int n = (int) (y*size[0]+x);
					nodes[n].velocity += p.velocity * w * (p.mass/nodes[n].mass);
				}
			}
		}
	}
	
	/* Implicit update only:
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
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					grad += nodes[(int) (y*size[0]+x)].velocity.trans_product(p.weight_gradient[idx]);
				}
			}
		}
	}
	*/
}
//Maps volume from the grid to particles
//This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes() const{
	//Estimate each particles volume (for force calculations)
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		//First compute particle density
		p.density = 0;
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Node density is trivial
					p.density += w * nodes[(int) (y*size[0]+x)].mass;
				}
			}
		}
		p.density /= node_volume;
		//Volume for each particle can be found from density
		p.volume = p.mass / p.density;
	}
}
//Calculate next timestep velocities for use in implicit integration
void Grid::explicitVelocities(const Vector2f& gravity){
	for (int i=0; i<obj->size; i++){
		Particle& p = obj->particles[i];
		//Solve for grid internal forces
		Matrix2f force = p.cauchyStress();
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					//Weight the force onto nodes
					int n = (int) (y*size[0]+x);
					nodes[n].force += force*p.weight_gradient[idx];
					nodes[n].has_force = true;
				}
			}
		}
	}
	
	//Now we have all grid forces, compute velocities (euler integration)
	Vector2f delta_scale = Vector2f(TIMESTEP);
	delta_scale /= cellsize;
	for (int y=0, idx=0; y<size[1]; y++){
		for (int x=0; x<size[0]; x++, idx++){
			//Get grid node (equivalent to (y*size[0] + x))
			GridNode &node = nodes[idx];
			//Check to see if this node needs to be computed
			if (node.has_force){
				node.velocity_new = node.velocity + TIMESTEP*(node.force/node.mass + gravity);
				node.force.setPosition(0);
				node.has_force = false;

				//Collision response
				//TODO: make this work for arbitrary collision geometry
				const int border_layers = 4;
				Vector2f new_pos = node.velocity_new*delta_scale + Vector2f(x, y);
				//Left border, right border
				if (new_pos[0] < border_layers-1 || new_pos[0] > size[0]-border_layers){
					node.velocity_new[0] = -node.velocity_new[0];
				}
				//Bottom border
				if (new_pos[1] < border_layers-1 || new_pos[1] > size[1]-border_layers)
					node.velocity_new[1] = -node.velocity_new[1];
			}
		}
	}
}
//Solve linear system for implicit velocities
void Grid::implicitVelocities(){
	//IMPLICIT_RATIO
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
		//VISUALIZATION PURPOSES ONLY:
		//Recompute density
		p.density = 0;
		
		int ox = p.grid_position[0],
			oy = p.grid_position[1];
		for (int idx=0, y=oy-1, y_end=oy+2; y<=y_end; y++){
			for (int x=ox-1, x_end=ox+2; x<=x_end; x++, idx++){
				float w = p.weights[idx];
				if (w > BSPLINE_EPSILON){
					GridNode &node = nodes[(int) (y*size[0]+x)];
					//Particle in cell
					pic += node.velocity_new*w;
					//Fluid implicit particle
					flip += (node.velocity_new - node.velocity)*w;
					//Velocity gradient
					grad += node.velocity_new.dyadic_product(p.weight_gradient[idx]);
					//VISUALIZATION ONLY: Update density
					p.density += w * node.mass;
				}
			}
		}
		//Final velocity is a linear combination of PIC and FLIP components
		p.velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
		//VISUALIZATION: Update density
		p.density /= node_volume;
	}
}
