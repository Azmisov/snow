#ifndef SIMCONSTANTS_H
#define	SIMCONSTANTS_H

static const float
	TIMESTEP = 5e-4,			//Time between each simulation step (5e-4 = implicit, 5e-5 = explicit)
	FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC (.95)
	CRIT_COMPRESS = 1-2.5e-2,	//Fracture threshold for compression (1-2.5e-2)
	CRIT_STRETCH = 1+7.5e-3,	//Fracture threshold for stretching (1+7.5e-3)
	HARDENING = 10,				//How much plastic deformation strengthens material (10)
	DENSITY = 400,				//Density of snow in kg/m^3 (400)
	YOUNGS_MODULUS = 1.4e5,		//Young's modulus (springiness) (1.4e5)
	POISSONS_RATIO = .2,		//Poisson's ratio (transverse/axial strain ratio) (.2)
	IMPLICIT_RATIO = 0,			//Percentage that should be implicit vs explicit for velocity update
	MAX_IMPLICIT_ITERS = 30,	//Maximum iterations for the conjugate residual
	MAX_IMPLICIT_ERR = 1e4,		//Maximum allowed error for conjugate residual
	MIN_IMPLICIT_ERR = 1e-4;	//Minimum allowed error for conjugate residual

//Various compiler options
#define LIMIT_FPS 0
#define FPS 60
#define REDRAW_EVERY LIMIT_FPS ? 1/TIMESTEP/FPS : 1;
#define REALTIME_PLAYBACK 0
#define SUPPORTS_POINT_SMOOTH 0
#define SCREENCAST 0

#endif

