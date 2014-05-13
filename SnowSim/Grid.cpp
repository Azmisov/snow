#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* object){
	obj = object;
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	nodes = new GridNode[(int) size.product()];
}
Grid::Grid(const Grid& orig){}
Grid::~Grid(){}

//Maps mass and velocity to the grid
//Enable the "volumes" flag for the first iteration only
void Grid::initialize(){
	//Reset the grid
	memset(nodes, 0, sizeof(GridNode)*size.product());
	
	//Map particle data to grid
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		//Particle position to grid coordinates
		//This will give errors if the particle is outside the grid bounds
		p[i].grid_position = (p[i].position - origin)/cellsize;
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		
		//Shape function gives a blending radius of two;
		//so we do computations within a 2x2 square for each particle
		for (int idx=0, x=ox-1, x_end=x+2; x<=x_end; x++){
			float wx = Grid::interpolate(x-ox);
			for (int y=oy-1, y_end=y+2; y<=y_end; y++){
				//Final weight is the dyadic product of weights in each dimension
				float weight = !wx ? 0 : wx * Grid::interpolate(y-oy);
				//Store this weight, because we use it to interpolate forces
				p[i].weights[idx++] = weight;
				
				//Interpolate mass
				nodes[(int) (x*size[0]+y)].mass += weight*p[i].mass;
			}
		}
	}
	
	//We interpolate velocity afterwards, to conserve momentum
	for (int i=0; i<obj->size; i++){
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		for (int idx=0, x=ox-1, x_end=x+2; x<=x_end; x++){
			for (int y=oy-1, y_end=y+2; y<=y_end; y++){
				//Interpolate velocity
				int n = (int) (x*size[0]+y);
				nodes[n].velocity += p[i].weights[idx++] * p[i].velocity * p[i].mass / nodes[n].mass;
			}
		}
	}
}

//Maps volume from the grid to particles
//This should only be called once, at the beginning of the simulation
void Grid::calculateVolumes(){	
	//Estimate each particles volume (for force calculations)
	float node_volume = cellsize.product();
	Particle* p = obj->particles;
	for (int i=0; i<obj->size; i++){
		int ox = p[i].grid_position[0],
			oy = p[i].grid_position[1];
		//First compute particle density
		for (int idx=0, x=ox-1, x_end=x+2; x<=x_end; x++){
			for (int y=oy-1, y_end=y+2; y<=y_end; y++){
				//Node density is trivial
				int n = (int) (x*size[0]+y);
				float node_density = nodes[n].mass / node_volume;
				//Weighted sum of nodes
				p[i].density += p[i].weights[idx++] * node_density;
			}
		}
		//Volume for each particle can be found from density
		p[i].volume = p[i].mass / p[i].density;
	}
}

