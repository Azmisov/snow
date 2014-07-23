#ifndef GRID_H
#define	GRID_H

#include <math.h>
#include <cstring>
#include <stdio.h>
#include "PointCloud.h"
#include "Vector2f.h"
#include "SimConstants.h"

const float BSPLINE_EPSILON = 1e-4;
const float TWO_THIRDS = 2/3.0;

//Grid node data
typedef struct GridNode{
	float mass;
	bool active;
	Vector2f force,
		velocity,
		velocity_new;
	
	//All the following variables are used by the implicit linear solver
	bool imp_active;	//are we still solving for vf
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
	float node_area;
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
	void collisionGrid();
	void collisionParticles() const;
	
	//Cubic B-spline shape/basis/interpolation function
	//A smooth curve from (0,1) to (1,0)
	static float bspline(float w){
		/*
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
		 */
		float abs_w=fabs(w);
		
		if(abs_w<=1){
			float abs_w_cubed=abs_w*abs_w*abs_w;
			float abs_w_squared=w*w;
			return (float).5*abs_w_cubed-abs_w_squared+2/3.0;}
		else if(abs_w>1 && abs_w<=2){
			float abs_w_cubed=abs_w*abs_w*abs_w;
			float abs_w_squared=w*w;
			return -1/6.0*abs_w_cubed+abs_w_squared-(float)2*abs_w+4/3.0;}
		return (float)0;
	}
	//Slope of interpolation function
	static float bsplineSlope(float w){
		/*
		float abs_x = fabs(x), w;
		if (abs_x < 1)
			return 1.5*x*abs_x - 2*x;
		else if (x < 2)
			return -x*abs_x/2 + 2*x - 2*x/abs_x;
		else return 0;
		//Clamp between -2/3 and 2/3... if needed
		 */
		float abs_w=fabs(w);
		
		if(abs_w<=1){
			return (float)1.5*abs_w*w-(float)2*w;}
		else if(abs_w>1 && abs_w<=2){
			return -(float).5*w*abs_w+(float)2*w-(float)2*w/abs_w;}
		return (float)0;
	}
};

#endif