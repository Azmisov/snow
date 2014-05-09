#ifndef OBJECT_H
#define	OBJECT_H

#include <vector>
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
	//Get bounding box [xmin, xmax, ymin, ymax]
	void bounds(float bounds[4]);
	
	//Generate square of particles (starting at 0,0 spaced 1 unit apart)
	static PointCloud* createSquare(float mass, int dim){
		PointCloud *obj = new PointCloud();
		obj->size = dim*dim;
		obj->particles = new Particle[obj->size];
		
		//point cloud definition; simple square
		//we give each particle equal mass (TODO: make mass dependent on volume?)
		float m = obj->size/mass;
		for (int x=0, n=0; x<dim; x++){
			for (int y=0; y<dim; y++){
				Particle &p = obj->particles[n++];
				p.position.setPosition(x, y);
				p.velocity = 0;
				p.mass = m;
			}
		}
		return obj;
	}
};

#endif

