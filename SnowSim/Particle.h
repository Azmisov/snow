#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include "Point3f.h"

class Particle {
	
public:
	float x, y, volume, density, mass;
	Point3f velocity, gradient;
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();
	
	void setPosition(float x, float y);
	
	//Generate square of particles (starting at 0,0)
	static std::vector<Particle>* createSquare(int size, float spacing){
		//point cloud definition; simple square
		std::vector<Particle> *cloud = new std::vector<Particle>(size*size);
		for (int x=0, n=0; x<size; x++){
			for (int y=0; y<size; y++){
				(*cloud)[n++].setPosition(x*spacing, y*spacing);
			}
		}
		return cloud;
	}
};

#endif

