#ifndef OBJECT_H
#define	OBJECT_H

#include <vector>
#include "SimConstants.h"
#include "Particle.h"

class PointCloud {
public:
	int size;
	std::vector<Particle> particles;

	PointCloud();
	PointCloud(int cloud_size);
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
		PointCloud *obj = new PointCloud(ppdim*ppdim);

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
				obj->particles.push_back(Particle(
					Vector2f(x*spacing, y*spacing),
					Vector2f(6,0), m, lambda+1000*x+1000*y, mu+1000*x+1000*y //
				));
			}
		}
		return obj;
	}
};

#endif

