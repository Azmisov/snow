/* 
 * File:   Particle.cpp
 * Author: inygaard
 * 
 * Created on March 31, 2014, 2:26 PM
 */

#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

void Particle::setPosition(float x, float y, float z){
	this->x = x;
	this->y = y;
	this->z = z;
}

