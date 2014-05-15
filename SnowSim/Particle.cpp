#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

Vector2f Particle::stressForce(){
	/* Stress force on each particle is: -volume*cauchy_stress
		We transfer the force to the FEM grid using the shape function gradient
		cauchy_stress = (2u(Fe - Re) + y(Je - 1)Je(Fe^-T))*(Fe^T)/J
			Fe: elastic deformation gradient
			Fp: plastic deformation gradient
			Re: rotation term of Fe's polar decomposition
				equivalent to WV* from SVD decomposition
			J: determinant of Fe*Fp
			Je: determinant of Fe
			u/y: Lame parameters
	*/
}

