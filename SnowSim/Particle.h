#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include "Vector2f.h"
#include "Matrix2f.h"

class Particle {
public:
	float volume, mass, density;
	Vector2f position, velocity;
	//Lame parameters
	float lambda, mu;
	//Deformation gradient (elastic and plastic parts)
	Matrix2f def_elastic, def_plastic;
	//Cached SVD's for elastic deformation gradient
	Matrix2f svd_w, svd_v;
	Vector2f svd_e;
	//Grid interpolation weights
	Vector2f grid_position;
	Vector2f weight_gradient[16];
	float weights[16];
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();

	//Compute SVD for deformation gradient
	void updateSVD();
	//Update deformation gradient and Lame parameters
	void updateGradient();
	//Compute stress force; call after updateSVD()
	Matrix2f stressForce() const;
};

#endif

