#ifndef GRID_H
#define	GRID_H

#include "PointCloud.h"
#include "Vector2f.h"

typedef struct GridNode{
	float mass;
	Vector2f velocity;
} GridNode;

class Grid {
public:
	Vector2f origin, cellsize;
	GridNode* nodes;
	
	Grid(Vector2f pos, Vector2f dims, Vector2f cells);
	Grid(const Grid& orig);
	virtual ~Grid();

	void initialize(PointCloud* obj);
};

#endif

