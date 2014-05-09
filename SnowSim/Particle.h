#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include "Vector2f.h"

class Particle {
	
public:
	float volume, density, mass;
	Vector2f position, velocity, gradient;
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();
	
	void setPosition(float x, float y);
};

#endif

