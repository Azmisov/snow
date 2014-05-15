#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

void Particle::updateSVD(){
	def_elastic.svd(&svd_w, &svd_e, &svd_v);
}
void Particle::updateGradient(){
	
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
	float je = def_elastic.determinant(),
		j = je*def_plastic.determinant(),
		pc = lambda*je*(je-1);
	Matrix2f fet = def_elastic.transpose(),
			stress = 2*mu*(def_elastic - svd_w*svd_v.transpose())*fet;
	stress.diag_sum(pc);
	//Final force
	return -volume/j*stress;
}

