#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include "Vector2f.h"

class Particle {
public:
	float volume, mass;
	Vector2f position, velocity;
	//Deformation gradient (elastic and plastic parts)
	Vector2f grad_elastic, grad_plastic;
	//Grid interpolation weights
	Vector2f grid_position;
	float weights[16];
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();
	
	void setPosition(float x, float y);
};

#endif

