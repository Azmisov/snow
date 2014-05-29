#ifndef GRID_H
#define	GRID_H

#include "PointCloud.h"
#include "Vector2f.h"
#include "SimConstants.h"
#include <math.h>
#include <cstring>
#include <stdio.h>

const float BSPLINE_EPSILON = 1e-4;
const float TWO_THIRDS = 2/3.0;

//Grid node data
typedef struct GridNode{
	float mass;
	bool has_velocity;
	Vector2f force,
		velocity,
		velocity_new;
	
	//All the following variables are used by the implicit linear solver
	bool active;	//are we still solving for vf
	Vector2f err,	//error of estimate
		r,			//residual of estimate
		p,			//residual gradient? squared residual?
		Ep,	Er;		//yeah, I really don't know how this works...
	float rEr;		//r.dot(Er)
} GridNode;

class Grid {
public:
	Vector2f origin, size, cellsize;
	PointCloud* obj;	
	float node_volume;
	//Nodes: use (y*size[0] + x) to index, where zero is the bottom-left corner (e.g. like a cartesian grid)
	int nodes_length;
	GridNode* nodes;
	
	//Grid be at least one cell; there must be one layer of cells surrounding all particles
	Grid(Vector2f pos, Vector2f dims, Vector2f cells, PointCloud* obj);
	Grid(const Grid& orig);
	virtual ~Grid();

	//Map particles to grid
	void initializeMass();
	void initializeVelocities();
	//Map grid volumes back to particles (first timestep only)
	void calculateVolumes() const;
	//Compute grid velocities
	void explicitVelocities(const Vector2f& gravity);
	void implicitVelocities();
	void recomputeImplicitForces();
	//Map grid velocities back to particles
	void updateVelocities() const;
	
	//Collision detection
	void collision();
	
	//Cubic B-spline shape/basis/interpolation function
	//A smooth curve from (0,1) to (1,0)
	static float bspline(float x){
		x = fabs(x);
		float w;
		if (x < 1)
			w = x*x*(x/2 - 1) + TWO_THIRDS;
		else if (x < 2)
			w = x*(x*(-x/6 + 1) - 2) + 2*TWO_THIRDS;
		else return 0;
		//Clamp between 0 and 1... if needed
		if (w < BSPLINE_EPSILON) return 0;
		return w;
	}
	//Slope of interpolation function
	static float bsplineSlope(float x){
		float abs_x = fabs(x), w;
		if (abs_x < 1)
			return 1.5*x*abs_x - 2*x;
		else if (x < 2)
			return -x*abs_x/2 + 2*x - 2*x/abs_x;
		else return 0;
		//Clamp between -2/3 and 2/3... if needed
	}
};

#endif