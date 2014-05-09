#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

void Particle::setPosition(float x, float y){
	this->x = x;
	this->y = y;
}

