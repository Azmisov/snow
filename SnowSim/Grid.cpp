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
void Grid::initialize(){
	//Reset the grid
	//If the grid is sparsely filled, it may be better to reset individual nodes
	memset(nodes, 0, sizeof(GridNode)*size.product());
	
	//Map particle data to grid
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		//Particle position to grid coordinates
		//This will give errors if the particle is outside the grid bounds
		p[i].grid_position = (p[i].position - origin)/cellsize;
		float ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		
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
				p[i].weights[idx] = weight;
				
				//Weight gradient is a vector of partial derivatives
				p[i].weight_gradient[idx].setPosition(dx*wy, wx*dy);
				
				//Interpolate mass
				nodes[(int) (x*size[0]+y)].mass += weight*p[i].mass;
			}
		}
	}
	
	//We interpolate velocity afterwards, to conserve momentum
	for (int i=0; i<obj->size; i++){
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++){
				int n = (int) (x*size[0]+y);
				if (nodes[n].mass > BSPLINE_EPSILON){
					//Interpolate velocity
					nodes[n].velocity += p[i].velocity * p[i].weights[idx++] * (p[i].mass/nodes[n].mass);
					std::cout << "N" << n << ": " << nodes[n].mass << ", " << nodes[n].velocity[0] << ", " << nodes[n].velocity[1] << std::endl;
				}
			}
		}
	}
}
//Maps volume from the grid to particles
//This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes() const{	
	//Estimate each particles volume (for force calculations)
	float node_volume = cellsize.product();
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		//First compute particle density
		p[i].density = 0;
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++){
				int n = (int) (x*size[0]+y);
				if (nodes[n].mass > BSPLINE_EPSILON){
					//Node density is trivial
					float node_density = nodes[n].mass / node_volume;
					//Weighted sum of nodes
					p[i].density += p[i].weights[idx++] * node_density;
				}
			}
		}
		//Volume for each particle can be found from density
		p[i].volume = p[i].mass / p[i].density;
	}
}
//Calculate next timestep velocities for use in explicit/implicit integration
void Grid::calculateVelocities(const Vector2f& extern_forces){
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		//Solve for grid internal forces
		Matrix2f force = p[i].stressForce();
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				int n = (int) (x*size[0]+y);
				if (nodes[n].mass > BSPLINE_EPSILON){
					//Weight the force onto nodes
					nodes[n].force += force*p[i].weight_gradient[idx];
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
			node.velocity_new = node.velocity + TIMESTEP/node.mass*(node.force+extern_forces);
			node.force.setPosition(0);
			node.has_force = false;
		}
	}
}
//Map grid velocities back to particles
void Grid::updateVelocities() const{
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		//We calculate PIC and FLIP velocities separately
		Vector2f pic, flip = p[i].velocity;
		//Also keep track of velocity gradient
		Matrix2f& grad = p[i].velocity_gradient;
		grad.setData(0.0);
		
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		for (int idx=0, x=ox-1, x_end=ox+2; x<=x_end; x++){
			for (int y=oy-1, y_end=oy+2; y<=y_end; y++, idx++){
				GridNode &node = nodes[(int) (x*size[0]+y)];
				if (node.mass > BSPLINE_EPSILON){
					//Interpolation weight
					float weight = p[i].weights[idx];				
					//Particle in cell
					pic += node.velocity_new*weight;
					//Fluid implicit particle
					flip += (node.velocity_new - node.velocity)*weight;

					//Velocity gradient
					grad += node.velocity_new.trans_product(p[i].weight_gradient[idx]);
				}
			}
		}
		//Final velocity is a linear combination of PIC and FLIP components
		p[i].velocity = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
	}
}
