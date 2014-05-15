#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include "Vector2f.h"
#include "Matrix2f.h"

class Particle {
public:
	float volume, mass, density;
	Vector2f position, velocity;
	//Lame parameters
	float lambda, mu;
	//Deformation gradient (elastic and plastic parts)
	Matrix2f def_elastic, def_plastic;
	//Grid interpolation weights
	Vector2f grid_position;
	Vector2f weight_gradient[16];
	float weights[16];
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();
	
	Vector2f stressForce();
};

#endif

