#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells){
	origin = pos;
	cellsize = dims/cells;
	size = cells+1;
	nodes = new GridNode[(int) size.product()];
}
Grid::Grid(const Grid& orig){}
Grid::~Grid(){}

void Grid::initialize(PointCloud* obj){
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
				nodes[(int) (x*size[0]+y)] = weight*p[i].mass;
			}
		}
	}
	
	//We interpolate velocity afterwards, to conserve momentum
}

