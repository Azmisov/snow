#ifndef PARTICLE_H
#define	PARTICLE_H

#include <cmath>
#include "SimConstants.h"
#include "Vector2f.h"
#include "Matrix2f.h"

class Particle {
public:
	float volume, mass, density;
	Vector2f position, velocity;
	Matrix2f velocity_gradient;
	//Lame parameters
	float lambda, mu;
	//Deformation gradient (elastic and plastic parts)
	Matrix2f def_elastic, def_plastic;
	//Cached SVD's for elastic deformation gradient
	Matrix2f svd_w, svd_v;
	Vector2f svd_e;
	//Cached determinants
	float det_elastic, det_plastic;
	//Grid interpolation weights
	Vector2f grid_position;
	Vector2f weight_gradient[16];
	float weights[16];
	
	Particle();
	Particle(const Particle& orig);
	virtual ~Particle();

	//Computes various intermediate data to speed up future calculations
	void updateSVD();
	void updateDet();
	//Update deformation gradient
	void gradientUpdate();
	//Compute stress force; call after updateSVD() and updateDet()
	Matrix2f stressForce() const;
};

#endif

