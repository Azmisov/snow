#include "SIM_GridInterpolate.h"
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
	IMPLEMENT_DATAFACTORY(SIM_GridInterpolate);
}

//Constructor
SIM_GridInterpolate::SIM_GridInterpolate(const SIM_DataFactory *factory) : BaseClass(factory){}
SIM_GridInterpolate::~SIM_GridInterpolate(){}

//Gets node description data
const SIM_DopDescription* SIM_GridInterpolate::getDescription(){
	//Register input parameters
	static PRM_Template theTemplates[] = {
		PRM_Template(PRM_STRING, 1, &PRM_Name("particles", "Particles")),
		PRM_Template(PRM_STRING, 1, &PRM_Name("inputAttr", "Attribute")),
		PRM_Template(PRM_STRING, 1, &PRM_Name("outputAttr", "Destination Field")),
		PRM_Template()
	};

	return &SIM_DopDescription(
		true,					// true, to make this node a DOP
		"hdk_gridinterpolate",	// internal name
		"Grid Interpolate",		// node label
		"Solver",				// data name (for details view)
		classname(),			// type of this dop
		theTemplates			// input parameters
	);
}

//Do the interpolation calculations
bool SIM_GridInterpolate::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep){
	return true;
}
