#ifndef SIMCONSTANTS_H
#define	SIMCONSTANTS_H

static const float
	TIMESTEP = 1e-4,			//For implicit update use 5e-4
	FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC
	CRIT_COMPRESS = 1-2.5e-2,	//Fracture threshold for compression
	CRIT_STRETCH = 1+7.5e-3,	//Fracture threshold for stretching
	HARDENING = 10,				//How much plastic deformation strengthens material
	DENSITY = 240,				//Density of snow in kg/m^3 (Stomakhin suggests 400, I'm adjusting for 2D)
	YOUNGS_MODULUS = 1.4e5,		//Young's modulus (springiness)
	POISSONS_RATIO = .2,		//Poisson's ratio (transverse/axial strain ratio)
	IMPLICIT_RATIO = .5,		//Percentage that should be implicit vs explicit for velocity update
	MAX_IMPLICIT_ITERS = 30,	//Maximum iterations for the conjugate residual
	MAX_IMPLICIT_ERR = 1e-6;	//Maximum allowed squared error for conjugate residual

#endif

