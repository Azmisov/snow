#ifndef SIMCONSTANTS_H
#define	SIMCONSTANTS_H

static const float
	PARTICLE_DIAM = .0072,		//Diameter of each particle; smaller = higher resolution
	FRAMERATE = 1/1000.0,		//Frames per second
	FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC (.95)
	CRIT_COMPRESS = 1-2.4e-2,	//Fracture threshold for compression (1-2.5e-2)
	CRIT_STRETCH = 1+7.5e-3,	//Fracture threshold for stretching (1+7.5e-3)
	HARDENING = 8.0,			//How much plastic deformation strengthens material (10)
	DENSITY = 1.0,				//Density of snow in kg/m^2 (400)
	YOUNGS_MODULUS = 1.4e5,		//Young's modulus (springiness) (1.4e5)
	POISSONS_RATIO = .2,		//Poisson's ratio (transverse/axial strain ratio) (.2)
	IMPLICIT_RATIO = 0,			//Percentage that should be implicit vs explicit for velocity update
	MAX_IMPLICIT_ITERS = 30,	//Maximum iterations for the conjugate residual
	MAX_IMPLICIT_ERR = 1e4,		//Maximum allowed error for conjugate residual
	MIN_IMPLICIT_ERR = 1e-4;	//Minimum allowed error for conjugate residual

//Actual timestep is adaptive, based on grid resolution and max velocity
const float TIMESTEP = 1e-3;

//Various compiler options
#define WIN_SIZE 640
#define WIN_METERS 1
#define SLO_MO 0e6
#define LIMIT_FPS true
#define SUPPORTS_POINT_SMOOTH true
#define SCREENCAST false
#define SCREENCAST_DIR "../screencast/"

#endif

