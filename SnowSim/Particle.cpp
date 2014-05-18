#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Vector2f& pos, const Vector2f& vel, float mass, float lame_lambda, float lame_mu){
	position.setPosition(pos);
	velocity.setPosition(vel);
	this->mass = mass;
	lambda_s = lame_lambda;
	mu_s = lame_mu;
	//To start out with, we assume the deformation gradient is zero
	//Or in other words, all particle velocities are the same
	def_elastic.loadIdentity();
	def_plastic.loadIdentity();
	det_elastic = 1;
	det_plastic = 1;
	svd_e.setPosition(1, 1);
	svd_w.setData(0, 1, 1, 0);
	svd_v.setData(0, 1, 1, 0);
}
Particle::~Particle(){}

void Particle::updatePos(){
	//Simple euler integration
	position += TIMESTEP*velocity;
}
void Particle::updateGradient(){
	//So, initially we make all updates elastic
	velocity_gradient *= TIMESTEP;
	velocity_gradient.diag_sum(1);
	def_elastic.setData(velocity_gradient * def_elastic);
	Matrix2f f_all = def_elastic * def_plastic;
	
	//We compute the SVD decomposition
	//The singular values (basically a scale transform) tell us if 
	//the particle has exceeded critical stretch/compression
	def_elastic.svd(&svd_w, &svd_e, &svd_v);
	//Clamp singular values to within elastic region
	for (int i=0; i<2; i++){
		if (svd_e[i] < CRIT_COMPRESS)
			svd_e[i] = CRIT_COMPRESS;
		else if (svd_e[i] > CRIT_STRETCH)
			svd_e[i] = CRIT_STRETCH;
	}
	
	//Recompute elastic and plastic gradient
	//We're basically just putting the SVD back together again
	Matrix2f v_cpy(svd_v), w_cpy(svd_w);
	v_cpy.diag_product_inv(svd_e);
	w_cpy.diag_product(svd_e);
	def_plastic = v_cpy*svd_w.transpose()*f_all;
	def_elastic = w_cpy*svd_v.transpose();
	
	//Update lame parameters to account for hardening
	det_elastic = def_elastic.determinant();
	det_plastic = def_plastic.determinant();
	float scale = exp(HARDENING*(1-det_plastic));
	mu = mu_s*scale;
	lambda = lambda_s*scale;
}
Matrix2f Particle::cauchyStress(){
	/* Stress force on each particle is: -volume*cauchy_stress
		We transfer the force to the FEM grid using the shape function gradient
		cauchy_stress can be computed via: (2u(Fe - Re)*Fe^T + y(Je - 1)Je*I)/J
			I: identity matrix
			Fe: elastic deformation gradient
			Fp: plastic deformation gradient
			Re: rotation term of Fe's polar decomposition
				equivalent to WV* from SVD decomposition
			J: determinant of Fe*Fp
			Je: determinant of Fe
			u/y: Lame parameters
	 
		For an implicit solution, we also need to calculate the forces with unresolved velocities:
			(1/J) * [2u(Fe' - Re') + y(Je' - 1)Je'*Fe'^-T] * Fe^T
		We define Fe' = X*Fe (the same way we do in updateGradient()),
		such that X = (I + TIMESTEP*weight_gradient).

		Note: volume = volume_initial*J, so the J drops out
	*/
	
	Matrix2f temp = def_elastic - svd_w*svd_v.transpose();
	temp *= 2*mu;
	Matrix2f temp2 = temp*def_elastic.transpose();
	temp2.diag_sum(lambda*det_elastic*(det_elastic-1));
	return -volume * temp2;
	
	/* For implicit solution only:
	//First compute temporary elastic deformation gradient for next timestep
	velocity_gradient *= TIMESTEP;
	velocity_gradient.diag_sum(1);
	Matrix2f fep = velocity_gradient * def_elastic,
			fep_inv_trans = fep.inverse();
	fep_inv_trans.transpose();
	float fe_det = fep.determinant();
	//Compute SVD of temporary Fe
	//We'll use the member variables of the class, since they aren't being used anymore
	fep.svd(&svd_w, &svd_e, &svd_v);
	//Now put it all together
	fep -= svd_w*svd_v.transpose();
	fep *= 2*mu;
	fep_inv_trans *= lambda*fe_det*(fe_det-1);
	fep += fep_inv_trans;
	//Final force
	return -volume * fep * def_elastic.transpose();
	//*/
}