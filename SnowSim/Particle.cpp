#include "Particle.h"

Particle::Particle(){}
Particle::Particle(const Particle& orig){}
Particle::~Particle(){}

Vector2f Particle::stressForce(){
	//Stress force on each particle is: -volume*cauchy_stress
	//We transfer the force to the FEM grid using the shape function gradient
	//To calculate cauchy stress:
	// partial = 2u(Fe - Re) + y(Je - 1)JeFe^-T
	// stress = partial*Fe^T/J ???? what is J? why are these matrices?
}

