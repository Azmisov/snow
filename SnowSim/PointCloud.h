#ifndef OBJECT_H
#define	OBJECT_H

#include <vector>
#include "SimConstants.h"
#include "Particle.h"
#include "Shape.h"

#define AREA_EPSILON 1e-5

inline float random_number(float lo, float hi){
	return lo + rand() / (float) (RAND_MAX/(hi-lo));
}

class PointCloud {
public:
	int size;
	float max_velocity;
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
	
	//Merge two point clouds
	void merge(const PointCloud& other);
	//Get bounding box [xmin, xmax, ymin, ymax]
	void bounds(float bounds[4]);

	//Generate particles that fill a set of shapes
	static PointCloud* createShape(std::vector<Shape*>& snow_shapes, Vector2f velocity){
		//Compute area of all the snow shapes
		float area = 0;
		int len = snow_shapes.size(), num_shapes = 0;
		for (int i=0; i<len; i++){
			float indi_area = snow_shapes[i]->area();
			if (indi_area > AREA_EPSILON && snow_shapes[i]->vertices.size() > 2){
				num_shapes++;
				area += indi_area;
			}
		}
		//If there is no volume, we can't really do a snow sim
		if (area < AREA_EPSILON)
			return NULL;
		
		//Lame parameters
		float lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO)),
			mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
		
		//Otherwise, create our object
		//Calculate particle settings
		float particle_area = PARTICLE_DIAM*PARTICLE_DIAM,
			particle_mass = particle_area*DENSITY;
		int particles = area / particle_area;
		//Randomly scatter points
		PointCloud *obj = new PointCloud(particles);
		float bounds[4];
		int shape_num = 0, total_points = 0;
		for (int i=0; i<len; i++){
			float a = snow_shapes[i]->area();
			if (a > AREA_EPSILON && snow_shapes[i]->vertices.size() > 2){
				int points;
				//Points given to each shape is proportional to their area
				if (++shape_num < num_shapes)
					points = a*particles/area;
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
				srand(5);
				while (points_found != points){
					float tx = random_number(bounds[0], bounds[1]),
						ty = random_number(bounds[2], bounds[3]);
					//Check if this point is inside the shape
					if (snow_shapes[i]->contains(tx, ty)){
						//Randomly adjust hardness of snow
						float adjust = (rand()/(float) RAND_MAX)*10;
						//Make snow hardest on the outer edges
						adjust *= ((fabs(tx-cx)/cw) + (fabs(ty-cy)/ch))/2.0;
						//Add the snow particle
						obj->particles.push_back(Particle(
							Vector2f(tx, ty), Vector2f(velocity), particle_mass, lambda, mu
						));						
						points_found++;
					}
				}
			}
		}
		//Set initial max velocity
		obj->max_velocity = velocity.length_squared();
		
		return obj;
	}
};

#endif

