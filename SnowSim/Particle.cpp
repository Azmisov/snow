#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

void Particle::updateSVD(){
	def_elastic.svd(&svd_w, &svd_e, &svd_v);
	//Clamp singular values to within elastic region
	for (int i=0; i<2; i++){
		if (svd_e[i] < CRIT_COMPRESS)
			svd_e[i] = CRIT_COMPRESS;
		else if (svd_e[i] > CRIT_STRETCH)
			svd_e[i] = CRIT_STRETCH;
	}
}
void Particle::updateDet(){
	det_elastic = def_elastic.determinant();
	det_plastic = def_plastic.determinant();
}
void Particle::gradientUpdate(){
	//So, initially we make all updates elastic
	Matrix2f dx = TIMESTEP*velocity_gradient;
	dx.diag_sum(1);
	def_elastic.setData(def_elastic * dx);
	Matrix2f f_all = def_elastic * def_plastic;
	
	//We compute the SVD decomposition
	//The singular values (basically a scale transform) tell us if 
	//the particle has exceeded critical stretch/compression
	updateSVD();
	
	//Recompute elastic and plastic gradient
	//We're basically just putting the SVD back together again
	Matrix2f v_cpy(svd_v), w_cpy(svd_w);
	v_cpy.diag_product_inv(svd_e);
	w_cpy.diag_product(svd_e);
	def_plastic = v_cpy*svd_w.transpose()*f_all;
	def_elastic = w_cpy*svd_v.transpose();
	
	//Update lame parameters to account for hardening
	//TODO:
	
}
Matrix2f Particle::stressForce() const{
	/* Stress force on each particle is: -volume*cauchy_stress
		We transfer the force to the FEM grid using the shape function gradient
		cauchy_stress = (2u(Fe - Re)*Fe^T + y(Je - 1)Je*I)/J
			I: identity matrix
			Fe: elastic deformation gradient
			Fp: plastic deformation gradient
			Re: rotation term of Fe's polar decomposition
				equivalent to WV* from SVD decomposition
			J: determinant of Fe*Fp
			Je: determinant of Fe
			u/y: Lame parameters
	*/
	//I move all scalar operations to the left, to reduce Matrix-Scalar operations
	float j = det_elastic*det_plastic,
		pc = lambda*det_elastic*(det_elastic-1);
	Matrix2f fet = def_elastic.transpose(),
			stress = 2*mu*(def_elastic - svd_w*svd_v.transpose())*fet;
	stress.diag_sum(pc);
	//Final force
	return -volume/j*stress;
}

