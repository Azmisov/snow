#ifndef GRID_H
#define	GRID_H

#include "PointCloud.h"
#include "Vector2f.h"
#include <math.h>
#include <cstring>
#include <stdio.h>

typedef struct GridNode{
	float mass;
	Vector2f velocity;
} GridNode;

class Grid {
public:
	Vector2f origin, size, cellsize;
	GridNode* nodes;
	PointCloud* obj;
	
	Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* obj);
	Grid(const Grid& orig);
	virtual ~Grid();

	void initialize();
	void calculateVolumes();
	
	//Cubic B-spline shape/basis/interpolation function
	//A smooth curve from (0,1) to (1,0)
	static float interpolate(float x){
		x = fabs(x);
		if (0 <= x && x < 1)
			return x*x*(x/2 - 1) + 2/3;
		if (1 <= x && x < 2)
			return x*(x*(-x/6 + 1) - 2) + 4/3;
		return 0;
	}
};

#endif

