/* 
 * File:   Particle.h
 * Author: inygaard
 *
 * Created on March 31, 2014, 2:26 PM
 */

#ifndef PARTICLE_H
#define	PARTICLE_H

class Particle {
	
public:
	float x, y, z;
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();
	
	void setPosition(float x, float y, float z);
private:

};

#endif	/* PARTICLE_H */

