#ifndef SIMCONSTANTS_H
#define	SIMCONSTANTS_H

#include "Vector2f.h"

static const float
	PARTICLE_DIAM = .0072,		//Diameter of each particle; smaller = higher resolution
	FRAMERATE = 1/60.0,			//Frames per second
	CFL = .04,					//Adaptive timestep adjustment
	MAX_TIMESTEP = 5e-4,		//Upper timestep limit
	FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC (.95)
	CRIT_COMPRESS = 1-1.9e-2,	//Fracture threshold for compression (1-2.5e-2)
	CRIT_STRETCH = 1+7.5e-3,	//Fracture threshold for stretching (1+7.5e-3)
	HARDENING = 5.0,			//How much plastic deformation strengthens material (10)
	DENSITY = 100,				//Density of snow in kg/m^2 (400 for 3d)
	YOUNGS_MODULUS = 1.5e5,		//Young's modulus (springiness) (1.4e5)
	POISSONS_RATIO = .2,		//Poisson's ratio (transverse/axial strain ratio) (.2)
	IMPLICIT_RATIO = 0,			//Percentage that should be implicit vs explicit for velocity update
	MAX_IMPLICIT_ITERS = 30,	//Maximum iterations for the conjugate residual
	MAX_IMPLICIT_ERR = 1e4,		//Maximum allowed error for conjugate residual
	MIN_IMPLICIT_ERR = 1e-4,	//Minimum allowed error for conjugate residual
	STICKY = .9,				//Collision stickiness (lower = stickier)
	GRAVITY = -9.8;

//Actual timestep is adaptive, based on grid resolution and max velocity
extern float TIMESTEP;

//Hardening parameters
static const float
	LAMBDA = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
	MU = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);

//Various compiler options
#define WIN_SIZE 640
#define WIN_METERS 1
#define SLO_MO 0e6
#define LIMIT_FPS false
#define REALTIME false
#define SUPPORTS_POINT_SMOOTH true
#define SCREENCAST false
#define SCREENCAST_DIR "../screencast/"
#define ENABLE_IMPLICIT false

#endif

