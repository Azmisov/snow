#include "SIM_SnowSolver.h"
#include <GU/GU_DetailHandle.h>
#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_AttributeRef.h>
#include <GA/GA_Iterator.h>
#include <GA/GA_Types.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_GeometryCopy.h>
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

	SIM_GeometryCopy* geometry = (SIM_GeometryCopy*)obj->getNamedSubData("particles");
	if (!geometry)
	{
		return true;
	}
	
	//Get particle data
	//Do we use the attribute name???
	GU_DetailHandle gdh = geometry->getGeometry().getWriteableCopy();
	const GU_Detail* gdp_in = gdh.readLock(); // Must unlock later
	GU_Detail* gdp_out = gdh.writeLock();

	GA_ROAttributeRef p_ref_position = gdp_in->findPointAttribute("P");
	GA_ROHandleT<UT_Vector3> p_position(p_ref_position.getAttribute());

	GA_ROAttributeRef p_ref_mass = gdp_in->findPointAttribute("mass");
	GA_ROHandleT<float> p_mass(p_ref_mass.getAttribute());

	GA_ROAttributeRef p_ref_volume = gdp_in->findPointAttribute("vol");
	GA_ROHandleT<float> p_volume(p_ref_volume.getAttribute());

	GA_RWAttributeRef p_ref_density = gdp_out->findPointAttribute("density");
	GA_RWHandleF p_density(p_ref_density.getAttribute());

	GA_ROAttributeRef p_ref_vel = gdp_in->findPointAttribute("vel");
	GA_ROHandleT<UT_Vector3> p_vel(p_ref_vel.getAttribute());

	GA_ROAttributeRef p_ref_Fe = gdp_in->findPointAttribute("Fe");
	GA_ROHandleT<UT_Matrix3> p_Fe(p_ref_Fe.getAttribute());

	GA_ROAttributeRef p_ref_Fp = gdp_in->findPointAttribute("Fp");
	GA_ROHandleT<UT_Matrix3> p_Fp(p_ref_Fp.getAttribute());

	
	//Get grid data
	SIM_ScalarField *g_mass_field;
	SIM_DataArray g_mass_data;
	getMatchingData(g_mass_data, obj, "g_mass");	
	g_mass_field = SIM_DATA_CAST(g_mass_data[0], SIM_ScalarField);

	SIM_VectorField *g_force_field;
	SIM_DataArray g_force_data;
	getMatchingData(g_force_data, obj, "g_force");	
	g_force_field = SIM_DATA_CAST(g_force_data[0], SIM_VectorField);

	SIM_VectorField *g_nvel_field;
	SIM_DataArray g_nvel_data;
	getMatchingData(g_nvel_data, obj, "g_nvel");
	g_nvel_field = SIM_DATA_CAST(g_nvel_data[0], SIM_VectorField);

	SIM_VectorField *g_ovel_field;
	SIM_DataArray g_ovel_data;
	getMatchingData(g_ovel_data, obj, "g_ovel");
	g_ovel_field = SIM_DATA_CAST(g_ovel_data[0], SIM_VectorField);

	SIM_ScalarField *g_active_field;
	SIM_DataArray g_active_data;
	getMatchingData(g_active_data, obj, "g_active");	
	g_active_field = SIM_DATA_CAST(g_active_data[0], SIM_ScalarField);

	UT_VoxelArrayF* g_mass = g_mass_field->getField()->fieldNC();

	UT_VoxelArrayF* g_forceX = g_force_field->getField(0)->fieldNC();
	UT_VoxelArrayF* g_forceY = g_force_field->getField(1)->fieldNC();
	UT_VoxelArrayF* g_forceZ = g_force_field->getField(2)->fieldNC();

	UT_VoxelArrayF* g_nvelX = g_nvel_field->getField(0)->fieldNC();
	UT_VoxelArrayF* g_nvelY = g_nvel_field->getField(1)->fieldNC();
	UT_VoxelArrayF* g_nvelZ = g_nvel_field->getField(2)->fieldNC();

	UT_VoxelArrayF* g_ovelX = g_nvel_field->getField(0)->fieldNC();
	UT_VoxelArrayF* g_ovelY = g_nvel_field->getField(1)->fieldNC();
	UT_VoxelArrayF* g_ovelZ = g_nvel_field->getField(2)->fieldNC();

	UT_VoxelArrayF* g_active = g_active_field->getField()->fieldNC();
	

	//TODO: should we reset grid here, or in a separate node?
	//		we're resizing the grid separately, so we could just do the reset after that...
	
	GA_RWHandleT<UT_Vector> p_wh(gdp_out->findPointAttribute("p_w"));
	if (!p_wh.isValid())
	{
		GA_RWAttributeRef p_w;
		// Don't create this attribute from the particles, it must be created here.
		p_w = gdp_out->addFloatTuple(GA_ATTRIB_POINT, "p_w", 64, GA_Defaults(0.0));
		// There is no type for array attributes, but apparently a type is not required.
		// p_w.setTypeInfo(GA_TYPE_VECTOR);
		p_wh = GA_RWHandleT<UT_Vector>(p_w);
	}

	// Not sure how to make this work yet
	GA_RWHandleT<UT_VectorT<UT_Vector3> > p_wgh(gdp_out->findPointAttribute("p_wg"));
	if(!p_wgh.isValid())
	{
		GA_RWAttributeRef p_wg;
		p_wg = gdp_out->addFloatTuple(GA_ATTRIB_POINT, "p_wg", 64, GA_Defaults(0.0)); // UT_Vector3(0.0, 0.0, 0.0) doesn't work . . .
		p_wgh = GA_RWHandleT<UT_VectorT<UT_Vector3> >(p_wg);
	}

	/// STEP #1: Transfer mass to grid 	

	if (p_position.isValid())
	{	
		int num = 0;
		for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
		{
			const UT_Vector3 pos(p_position.get(it.getOffset())); //Particle position pos
			UT_Vector weights(0, 64);
			UT_VectorT<UT_Vector3> weight_gradients(0, 64);
			int p_gridx = 0;
			int p_gridy = 0;
			int p_gridz = 0;
			g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position

			for (int idx=0, z=p_gridz-1, z_end=p_gridz+2; z<=z_end; z++){
				//Z-dimension interpolation
				float z_pos = z-p_gridz,
					wz = SIM_SnowSolver::bspline(z_pos),
					dz = SIM_SnowSolver::bsplineSlope(z_pos);

				for (int y=p_gridy-1, y_end=p_gridy+2; y<=y_end; y++){
					//Y-dimension interpolation
					float y_pos = y-p_gridy,
						wy = SIM_SnowSolver::bspline(y_pos),
						dy = SIM_SnowSolver::bsplineSlope(y_pos);
					for (int x=p_gridx-1, x_end=p_gridx+2; x<=x_end; x++, idx++){
						//X-dimension interpolation
						float x_pos = x-p_gridx,
							wx = SIM_SnowSolver::bspline(x_pos),
							dx = SIM_SnowSolver::bsplineSlope(x_pos);
						
						//Final weight is dyadic product of weights in each dimension
						float weight = wx*wy;
						weights.assign(&weight, idx, idx);

						//Weight gradient is a vector of partial derivatives
						const UT_Vector3 newWeight(UT_Vector3(dx*wx, dy*wy, dz*wz));		 
						weight_gradients.setSubvector3(idx, newWeight);

						//Interpolate mass
					    float node_mass = g_mass->getValue(x,y,z);
						float particle_mass(p_mass.get(it.getOffset()));
						node_mass += weight*particle_mass; 
						g_mass->setValue(x,y,z,node_mass);
					}
				}
			}
			//Causes a crash!
			//p_wh.set(it.getOffset(), weights);
			//p_wgh.set(it.getOffset(), weight_gradients);
		}
	}

	/// STEP #2: First timestep only - Estimate particle volumes using grid mass
	if(time == 0.0){
		
		for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
		{
			const UT_Vector3 pos(p_position.get(it.getOffset())); //Particle position pos

			float density = 0;				
			int p_gridx = 0;
			int p_gridy = 0;
			int p_gridz = 0;
			g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position
		
			for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
				for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
					for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
						float w = .1;//p.weights[idx];
						if (w > BSPLINE_EPSILON){
							density += w * (g_mass->getValue(x,y,z));
						}
					}
				}
			}

			float oldDensity(p_density.get(it.getOffset()));
			oldDensity /= .1*.1*.1;

			//This seems to only apply to the current call
			p_density.set(it.getOffset(),oldDensity);

			//p.density /= .1*.1*.1;			
			//Have division size as parameter!!
			//p.volume = p.mass / p.density;
			//p_mass.set(it.getOffset(),.1); //Particle position pos
		}
	}
	int num = 0;
	for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
	{
		//This prints out 1000 for density on the 1st timestep (calculated above on that timstep), but after goes back to the default of 1
		float oldDensity(p_density.get(it.getOffset()));
		if(num == 100){
			cout << oldDensity << endl;
		}
		num++;
	}

	/// STEP #3: Transfer velocity to grid
	
	/// STEP #4: Compute new grid velocities
	
	/// STEP #5: Grid collision resolution
	
	/// STEP #6: Transfer grid velocities to particles and integrate
	
	/// STEP #7: Update particle deformation gradient
	
	/// STEP #8: Particle collision resolution

	gdh.unlock(gdp_out);
    gdh.unlock(gdp_in);
	
	return true;
}






