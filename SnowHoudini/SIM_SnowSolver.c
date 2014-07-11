#include "SIM_SnowSolver.h"
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>

#include <iostream>
#include <stdio.h>
#include <string>

//Houdini hook
void initializeSIM(void *){
	IMPLEMENT_DATAFACTORY(SIM_SnowSolver);
}

//Constructor
SIM_SnowSolver::SIM_SnowSolver(const SIM_DataFactory *factory) : BaseClass(factory){}
SIM_SnowSolver::~SIM_SnowSolver(){}

//Gets node description data
const SIM_DopDescription* SIM_SnowSolver::getDescription(){
	//TODO: how to add a long description for these parameters, when user hovers over them?
	//TODO: have a horizontal rule separating lagrangian and eulerian parameters (like how the Gas Linear Combination subsolver has it)
	//TODO: maybe move some of the particles attributes into SIM data fields instead of having them attached to particles?
	//TODO: make sim constansts parameters???
	
	/* Particle parameters (lagrangian):
		We will just use the default point attribute "P" for particle position;
		
		Also, we can just use temporary matrices for SVD/Polar Decomp/Velocity Grad,
		since they get recomputed at each timestep anyways. If we split this up into
		multiple sub-solver nodes, we'd need to make them particle attributes.
	*/
	static PRM_Name p_field("particles", "Particles");			//particles
	static PRM_Name p_fe("p_fe", "Fe Attr");					//particle elastic deformation gradient
	static PRM_Name p_fp("p_fp", "Fp Attr");					//particle plastic deformation gradient
	static PRM_Name p_mass("p_mass", "Mass Attr");				//particle mass
	static PRM_Name p_vel("p_vel", "Velocity Attr");			//particle velocity
	static PRM_Name p_vol("p_vol", "Volume Attr");				//particle volume
	static PRM_Name p_d("p_d", "Density Attr");					//particle density
	//It may be better to remove these, if recomputing these values is actually faster than caching
	static PRM_Name p_w("p_w", "Weights Attr");					//particle weight (for each node within 2-node radius)
	static PRM_Name p_wg("p_wg", "Weight Gradients Attr");		//particle weight gradient (for each node within 2-node radius)

	//Grid parameters (eulerian):
	static PRM_Name g_mass("g_mass", "Mass Field");				//grid mass
	static PRM_Name g_force("g_force", "Force Field");			//grid forces
	static PRM_Name g_nvel("g_nvel", "New Velocity Field");		//grid velocity (after applying forces)
	static PRM_Name g_ovel("g_ovel", "Old Velocity Field");		//grid velocity (before applying forces)
	static PRM_Name g_active("g_active", "Activated Field");	//boolean field that tells whether there are particles within a radius of 2
	
	//TODO: import collision field

	static PRM_Template theTemplates[] = {
		//particles
		PRM_Template(PRM_STRING, 1, &p_field),
		PRM_Template(PRM_STRING, 1, &p_fe),
		PRM_Template(PRM_STRING, 1, &p_fp),
		PRM_Template(PRM_STRING, 1, &p_mass),
		PRM_Template(PRM_STRING, 1, &p_vel),
		PRM_Template(PRM_STRING, 1, &p_vol),
		PRM_Template(PRM_STRING, 1, &p_d),
		PRM_Template(PRM_STRING, 1, &p_w),
		PRM_Template(PRM_STRING, 1, &p_wg),
		//grid
		PRM_Template(PRM_STRING, 1, &g_mass),
		PRM_Template(PRM_STRING, 1, &g_force),
		PRM_Template(PRM_STRING, 1, &g_nvel),
		PRM_Template(PRM_STRING, 1, &g_ovel),
		PRM_Template(PRM_STRING, 1, &g_active),
		//what is this for ???
		PRM_Template()
	};

	static SIM_DopDescription desc(
		true,					// true, to make this node a DOP
		"hdk_SnowSolver",		// internal name
		"Snow Solver",			// node label
		"Solver",				// data name (for details view)
		classname(),			// type of this dop
		theTemplates			// input parameters
	);
	return &desc;
}

//Do the interpolation calculations
bool SIM_SnowSolver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep){
	/// STEP #0: Retrieve all data objects from Houdini
	//TODO: should we reset grid here, or in a separate node?
	//		we're resizing the grid separately, so we could just do the reset after that...
	

	/// STEP #1: Transfer mass to grid
	
	/// STEP #2: First timestep only - Estimate particle volumes using grid mass
	
	/// STEP #3: Transfer velocity to grid
	
	/// STEP #4: Compute new grid velocities
	
	/// STEP #5: Grid collision resolution
	
	/// STEP #6: Transfer grid velocities to particles and integrate
	
	/// STEP #7: Update particle deformation gradient
	
	/// STEP #8: Particle collision resolution
	
	return true;
}






