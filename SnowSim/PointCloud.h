#ifndef OBJECT_H
#define	OBJECT_H

#include <vector>
#include "SimConstants.h"
#include "Particle.h"
#include "Shape.h"

#define VOLUME_EPSILON 1e-5

inline float random_number(float lo, float hi){
	 return lo + rand() / (float) (RAND_MAX/(hi-lo));
}

class PointCloud {
public:
	int size;
	float mass;
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
	
	//Generate square of particles (starting at 0,0)
	static PointCloud* createSquare(float mpdim, int ppdim, Vector2f velocity){
		PointCloud *obj = new PointCloud(ppdim*ppdim);

		float spacing = ppdim == 1 ? 0 : mpdim/(ppdim-1);
		//point cloud definition; simple square
		//we give each particle equal mass (TODO: make mass dependent on volume?)
		float volume = mpdim*mpdim*mpdim;
		obj->mass = DENSITY*volume;
		float m = volume/obj->size;
		//we also give each particle equal lame parameters (TODO: randomized/customizable lames)
		float lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
			mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
		for (int x=0, n=0; x<ppdim; x++){
			for (int y=0; y<ppdim; y++){
				float adjust = (rand() % 6000)*(fabs(x-ppdim/2.0) + fabs(y-ppdim/2.0))/ppdim;
				obj->particles.push_back(Particle(
					Vector2f(x*spacing, y*spacing),
					Vector2f(velocity), m, lambda+adjust, mu+adjust
				));
			}
		}
		return obj;
	}
	
	//Generate particles that fill a set of shapes
	static PointCloud* createShape(std::vector<Shape*>& snow_shapes, int particles, Vector2f velocity){
		//Compute area of all the snow shapes
		float volume = 0;
		int len = snow_shapes.size(), num_shapes = 0;
		for (int i=0; i<len; i++){
			float indi_volume = snow_shapes[i]->volume();
			if (indi_volume > VOLUME_EPSILON && snow_shapes[i]->vertices.size() > 2){
				num_shapes++;
				volume += indi_volume;
			}
		}
		//If there is no volume, we can't really do a snow sim
		if (volume < 1e-5)
			return NULL;
		
		//Lame parameters
		float lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
			mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
		
		//Otherwise, create our object
		PointCloud *obj = new PointCloud(particles);
		//Use volume to estimate mass per particle
		obj->mass = volume*DENSITY;
		float mass = obj->mass/particles;
		//Randomly scatter points
		float bounds[4];
		int shape_num = 0, total_points = 0;
		for (int i=0; i<len; i++){
			float v = snow_shapes[i]->volume();
			if (v > VOLUME_EPSILON && snow_shapes[i]->vertices.size() > 2){
				int points;
				//Points given to each shape is proportional to their area
				if (++shape_num < num_shapes)
					points = v*particles/volume;
				//Last shape gets remainder, so we don't have round-off errors
				else points = particles-total_points;
				total_points += points;
				
				//Estimate the centroid of the shape with the bounds
				snow_shapes[i]->bounds(bounds);
				float cx = (bounds[0]+bounds[1])/2.0,
					cy = (bounds[2]+bounds[3])/2.0,
					cw = bounds[1] - cx,
					ch = bounds[3] - cy;
				
				//Randomly scatter points in the shape until the quota is met
				int points_found = 0;
				while (points_found != points){
					float tx = random_number(bounds[0], bounds[1]),
						ty = random_number(bounds[2], bounds[3]);
					//Check if this point is inside the shape
					if (snow_shapes[i]->contains(tx, ty)){
						//Randomly adjust hardness of snow
						float adjust = (rand()/(float) RAND_MAX)*70;
						//Make snow hardest on the outer edges
						adjust *= ((fabs(tx-cx)/cw) + (fabs(ty-cy)/ch))/2.0;
						obj->particles.push_back(Particle(
							Vector2f(tx, ty), Vector2f(velocity), mass, lambda*adjust, mu*adjust
						));
						points_found++;
					}
				}
			}
		}
		
		return obj;
	}
};

#endif

