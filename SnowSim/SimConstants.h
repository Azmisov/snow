#ifndef SIMCONSTANTS_H
#define	SIMCONSTANTS_H

static const float
	TIMESTEP = .01,
	FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC
	CRIT_COMPRESS = 0.975,		//Fracture threshold for compression
	CRIT_STRETCH = 1.0075,		//Fracture threshold for stretching
	HARDENING = 10;				//How much plastic deformation strengthens material

#endif

