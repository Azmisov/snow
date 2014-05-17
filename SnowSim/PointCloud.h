#ifndef OBJECT_H
#define	OBJECT_H

#include <vector>
#include "SimConstants.h"
#include "Particle.h"

class PointCloud {
public:
	int size;
	Particle* particles;
	
	PointCloud();
	PointCloud(const PointCloud& orig);
	virtual ~PointCloud();
	
	//Transform points in cloud
	void scale(Vector2f origin, Vector2f scale);
	void translate(Vector2f off);
	//Update particle data
	void update();
	
	//Get bounding box [xmin, xmax, ymin, ymax]
	void bounds(float bounds[4]);
	
	//Generate square of particles (starting at 0,0 spaced 1 unit apart)
	static PointCloud* createSquare(float mpdim, int ppdim){
		PointCloud *obj = new PointCloud();
		obj->size = ppdim*ppdim;
		obj->particles = new Particle[obj->size];

		float spacing = ppdim == 1 ? 0 : mpdim/(ppdim-1);
		//point cloud definition; simple square
		//we give each particle equal mass (TODO: make mass dependent on volume?)
		float volume = mpdim*mpdim,	//TODO: this is 2D, will need to adjust for 3D
			m = (DENSITY*volume)/obj->size;
		//we also give each particle equal lame parameters (TODO: randomized/customizable lames)
		float lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
			mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
		for (int x=0, n=0; x<ppdim; x++){
			for (int y=0; y<ppdim; y++){
				Particle &p = obj->particles[n++];
				p.position.setPosition(x*spacing, y*spacing);
				p.velocity.setPosition(0);
				p.lambda_s = lambda;
				p.mu_s = mu;
				p.mass = m;
			}
		}
		return obj;
	}
};

#endif

