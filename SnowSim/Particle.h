#ifndef PARTICLE_H
#define	PARTICLE_H

#include <cmath>
#include "Vector2f.h"
#include "Matrix2f.h"
#include "SimConstants.h"

class Particle {
public:
	float volume, mass, density;
	Vector2f position, velocity;
	Matrix2f velocity_gradient;
	//Lame parameters (_s denotes starting configuration)
	float lambda, mu;
	//Deformation gradient (elastic and plastic parts)
	Matrix2f def_elastic, def_plastic;
	//Cached SVD's for elastic deformation gradient
	Matrix2f svd_w, svd_v;
	Vector2f svd_e;
	//Cached polar decomposition
	Matrix2f polar_r, polar_s;
	//Grid interpolation weights
	Vector2f grid_position;
	Vector2f weight_gradient[16];
	float weights[16];

	Particle();
	Particle(const Vector2f& pos, const Vector2f& vel, float mass, float lambda, float mu);
	virtual ~Particle();

	//Update position, based on velocity
	void updatePos();
	//Update deformation gradient
	void updateGradient();
	void applyPlasticity();
	//Compute stress tensor
	const Matrix2f energyDerivative();
	
	//Computes stress force delta, for implicit velocity update
	const Vector2f deltaForce(const Vector2f& u, const Vector2f& weight_grad);
};

#endif

