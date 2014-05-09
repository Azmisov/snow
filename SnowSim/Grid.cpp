#include "Grid.h"

Grid::Grid(Vector2f pos, Vector2f dims, Vector2f cells){
	origin = pos;
	cellsize = dims/cells;
	nodes = new GridNode[(int) (cells[0]*cells[1])];
}
Grid::Grid(const Grid& orig){}
Grid::~Grid(){}

void map(PointCloud* obj);

