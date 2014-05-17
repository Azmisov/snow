#ifndef GRID_H
#define	GRID_H

#include "PointCloud.h"
#include "Vector2f.h"
#include "SimConstants.h"
#include <math.h>
#include <cstring>
#include <stdio.h>

const float BSPLINE_EPSILON = 1e-6;

typedef struct GridNode{
	float mass;
	bool has_force;
	Vector2f force,
		velocity,
		velocity_new;
} GridNode;

class Grid {
public:
	Vector2f origin, size, cellsize;
	GridNode* nodes;
	PointCloud* obj;
	
	Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* obj);
	Grid(const Grid& orig);
	virtual ~Grid();

	//Map particles to grid
	void initialize();
	//Map grid volumes back to particles (first timestep only)
	void calculateVolumes() const;
	//Compute grid velocities
	void calculateVelocities(const Vector2f& extern_forces);
	//Map grid velocities back to particles
	void updateVelocities() const;
	
	//Cubic B-spline shape/basis/interpolation function
	//A smooth curve from (0,1) to (1,0)
	static float bspline(float x){
		if (x < 1)
			return x*x*(x/2 - 1) + 2/3.0;
		if (x < 2)
			return x*(x*(-x/6 + 1) - 2) + 4/3.0;
		return 0;
	}
	//Slope of interpolation function
	static float bsplineSlope(float x){
		if (x < 1)
			return x/2*(3*x - 4);
		if (x < 2)
			return x/2*(4 - x) - 2;
		return 0;
	}
};

#endif

