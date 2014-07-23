#include "SIM_GridInterpolate.h"
#include <GU/GU_DetailHandle.h>
#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_AttributeRef.h>
#include <GA/GA_Iterator.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_GeometryCopy.h>
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
	static PRM_Name p_field("particles", "Particles");
	static PRM_Name i_field("inputAttr", "Attribute");
	static PRM_Name o_field("outputAttr", "Destination Field");

	static PRM_Template theTemplates[] = {
		PRM_Template(PRM_STRING, 1, &p_field),
		PRM_Template(PRM_STRING, 1, &i_field),
		PRM_Template(PRM_STRING, 1, &o_field),
		PRM_Template()
	};

	static SIM_DopDescription desc(
		true,					// true, to make this node a DOP
		"hdk_gridinterpolate",	// internal name
		"Grid Interpolate",		// node label
		"Solver",				// data name (for details view)
		classname(),			// type of this dop
		theTemplates			// input parameters
	);
	return &desc;
}

//Do the interpolation calculations
bool SIM_GridInterpolate::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep){
	SIM_GeometryCopy* geometry = (SIM_GeometryCopy*)obj->getNamedSubData(INPUT_PARTICLES_NAME);
	if (geometry)
	{
		// Get voxel field
		SIM_ScalarField *field;
		SIM_DataArray field_data;
		getMatchingData(field_data, obj, GRID_FIELD_NAME);
		
		field = SIM_DATA_CAST(field_data[0], SIM_ScalarField);

		// Get geometry attributes
		GU_DetailHandle gdh = geometry->getGeometry().getWriteableCopy();
		const GU_Detail* gdp_in = gdh.readLock(); // Must unlock later
		GU_Detail* gdp_out = gdh.writeLock();

		UT_String input_attr;
		getInputAttr(input_attr);
		GA_ROAttributeRef in_attr_ref = gdp_in->findPointAttribute(input_attr);
		GA_ROHandleT<fpreal> in_attr(in_attr_ref.getAttribute());

		// Not sure how to write to particle attributes yet . . .
		// UT_String output_attr;
		// getOutputAttr(output_attr);
		// GA_RWAttributeRef out_attr_ref = gdp_out->findPointAttribute(output_attr);

		GA_ROAttributeRef pos = gdp_in->findPointAttribute("P"); // I think P always exists
		GA_ROHandleT<UT_Vector3> point_pos(pos.getAttribute());

		if (in_attr.isValid() && point_pos.isValid())
		{
			for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance())
			{
				UT_Vector3 P(point_pos.get(it.getOffset())); // Position of the particle in world space
				fpreal val = in_attr.get(it.getOffset()); // Value of the particle input attribute

				int grid_x, grid_y, grid_z; // Indices to the voxel this particle is in
				UT_Vector3 grid_P; // Position of the center of this voxel in world space
				UT_Vector3 voxel_size = field->getVoxelSize();
				field->posToIndex(P, grid_x, grid_y, grid_z);
				field->indexToPos(grid_x, grid_y, grid_z, grid_P);
			}
		}
		gdh.unlock(gdp_out);
	    gdh.unlock(gdp_in);
	}
	return true;
}
