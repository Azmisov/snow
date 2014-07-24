#include "SIM_SnowSolver.h"
#include <GU/GU_DetailHandle.h>
#include <GU/GU_Detail.h>
#include <GA/GA_Handle.h>
#include <GA/GA_AttributeRef.h>
#include <GA/GA_Iterator.h>
#include <GA/GA_Types.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_MatrixSolver.h>
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
#include <vector>

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
	static PRM_Name p_mass("p_mass", "Particle Mass");			//particle mass
	static PRM_Name p_vel("p_vel", "Velocity Attr");			//particle velocity
	static PRM_Name p_vol("p_vol", "Volume Attr");				//particle volume
	static PRM_Name p_d("p_d", "Density Attr");					//particle density
	//It may be better to remove these, if recomputing these values is actually faster than caching
	static PRM_Name p_w("p_w", "Weights Attr");					//particle weight (for each node within 2-node radius)
	static PRM_Name p_wg("p_wg", "Weight Gradients Attr");		//particle weight gradient (for each node within 2-node radius)

	//Grid parameters (eulerian):
	static PRM_Name g_div("div_size", "Division Size");			//grid division size
	static PRM_Name g_mass("g_mass", "Mass Field");				//grid mass
	static PRM_Name g_nvel("g_nvel", "New Velocity Field");		//grid velocity (after applying forces)
	static PRM_Name g_ovel("g_ovel", "Old Velocity Field");		//grid velocity (before applying forces)
	static PRM_Name g_active("g_active", "Activated Field");	//boolean field that tells whether there are particles within a radius of 2
	
	//TODO: import collision field

	static PRM_Template theTemplates[] = {
		//particles
		PRM_Template(PRM_STRING, 1, &p_field),
		PRM_Template(PRM_STRING, 1, &p_fe),
		PRM_Template(PRM_STRING, 1, &p_fp),
		PRM_Template(PRM_FLT_J, 1, &p_mass),
		PRM_Template(PRM_STRING, 1, &p_vel),
		PRM_Template(PRM_STRING, 1, &p_vol),
		PRM_Template(PRM_STRING, 1, &p_d),
		PRM_Template(PRM_STRING, 1, &p_w),
		PRM_Template(PRM_STRING, 1, &p_wg),
		//grid
		PRM_Template(PRM_FLT_J, 1, &g_div),
		PRM_Template(PRM_STRING, 1, &g_mass),
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

	// SIM_GeometryCopy* geometry = (SIM_GeometryCopy*)obj->getNamedSubData("particles");
	SIM_Geometry* geometry = (SIM_Geometry*)obj->getNamedSubData("particles");
	if (!geometry)
	{
		return true;
	}
	
	//Get particle data
	//Do we use the attribute name???
	// GU_DetailHandle gdh = geometry->getGeometry().getWriteableCopy();
	GU_DetailHandle gdh = geometry->getOwnGeometry();
	const GU_Detail* gdp_in = gdh.readLock(); // Must unlock later
	GU_Detail* gdp_out = gdh.writeLock();

	GA_ROAttributeRef p_ref_position = gdp_in->findPointAttribute("P");
	GA_ROHandleT<UT_Vector3> p_position(p_ref_position.getAttribute());

	GA_RWAttributeRef p_ref_volume = gdp_out->findPointAttribute("vol");
	GA_RWHandleT<float> p_volume(p_ref_volume.getAttribute());

	GA_RWAttributeRef p_ref_density = gdp_out->findPointAttribute("density");
	GA_RWHandleF p_density(p_ref_density.getAttribute());

	GA_ROAttributeRef p_ref_vel = gdp_in->findPointAttribute("vel");
	GA_ROHandleT<UT_Vector3> p_vel(p_ref_vel.getAttribute());

	GA_ROAttributeRef p_ref_Fe = gdp_in->findPointAttribute("Fe");
	GA_ROHandleT<UT_Matrix3> p_Fe(p_ref_Fe.getAttribute());

	GA_ROAttributeRef p_ref_Fp = gdp_in->findPointAttribute("Fp");
	GA_ROHandleT<UT_Matrix3> p_Fp(p_ref_Fp.getAttribute());

	//EVALUATE PARAMETERS
	float particle_mass = .01;
	float division_size = .1;
	float mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
	float lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO));

	//Posible other parameters:
	//lambda_s (lame_lambda)
	//mu_s (lame_mu)

	float voxelArea = division_size*division_size*division_size;

	//Get grid data
	SIM_ScalarField *g_mass_field;
	SIM_DataArray g_mass_data;
	getMatchingData(g_mass_data, obj, "g_mass");	
	g_mass_field = SIM_DATA_CAST(g_mass_data[0], SIM_ScalarField);

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

	UT_VoxelArrayF* g_nvelX = g_nvel_field->getField(0)->fieldNC();
	UT_VoxelArrayF* g_nvelY = g_nvel_field->getField(1)->fieldNC();
	UT_VoxelArrayF* g_nvelZ = g_nvel_field->getField(2)->fieldNC();

	UT_VoxelArrayF* g_ovelX = g_ovel_field->getField(0)->fieldNC();
	UT_VoxelArrayF* g_ovelY = g_ovel_field->getField(1)->fieldNC();
	UT_VoxelArrayF* g_ovelZ = g_ovel_field->getField(2)->fieldNC();

	UT_VoxelArrayF* g_active = g_active_field->getField()->fieldNC();
	
	UT_Vector3 fieldDims(g_nvel_field->getDivisions());
	//TODO: should we reset grid here, or in a separate node?
	//		we're resizing the grid separately, so we could just do the reset after that...
	
	// GA_RWHandleT<UT_Vector> p_wh(gdp_out->findPointAttribute("p_w"));
	// if (!p_wh.isValid())
	// {
	// 	GA_RWAttributeRef p_w;
	// 	// Don't create this attribute from the particles, it must be created here.
	// 	p_w = gdp_out->addFloatTuple(GA_ATTRIB_POINT, "p_w", 64, GA_Defaults(0.0));
	// 	// There is no type for array attributes, but apparently a type is not required.
	// 	// p_w.setTypeInfo(GA_TYPE_VECTOR);
	// 	p_wh = GA_RWHandleT<UT_Vector>(p_w);

	// 	// DEBUG ==========================================================================
	// 	std::string filepath = "~/snow_solver_test.txt";
	// 	FILE *tfile = fopen(filepath.c_str(), "w");
	// 	fputs("Added tuple attribute\n", tfile);
	// 	fclose(tfile);
	// 	// DEBUG ==========================================================================
	// }

	// Not sure how to make this work yet
	// GA_RWHandleT<UT_VectorT<UT_Vector3> > p_wgh(gdp_out->findPointAttribute("p_wg"));
	// if(!p_wgh.isValid())
	// {
	// 	GA_RWAttributeRef p_wg;
	// 	p_wg = gdp_out->addFloatTuple(GA_ATTRIB_POINT, "p_wg", 64, GA_Defaults(0.0)); // UT_Vector3(0.0, 0.0, 0.0) doesn't work . . .
	// 	p_wgh = GA_RWHandleT<UT_VectorT<UT_Vector3> >(p_wg);
	// }

	int point_count = gdp_out->getPointRange().getEntries();
	int weight_count = 64;
	double p_w[point_count][weight_count];
	
	// UT_Vector3 p_wgh[point_count][64]; // Doesn't work with C99
	
	std::vector<std::vector<UT_Vector3> > p_wgh;
	for(int i=0; i<point_count; i++)
	{
		std::vector<UT_Vector3> empty;
		for(int j=0; j<weight_count; j++)
		{
			empty.push_back(UT_Vector3(0.0, 0.0, 0.0));
		}
		p_wgh.push_back(empty);
	}
	
	/// STEP #1: Transfer mass to grid 	

	if (p_position.isValid())
	{	
		for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
		{
			const UT_Vector3 pos(p_position.get(it.getOffset())); //Particle position pos
			UT_Vector weights(0, 64);
			UT_VectorT<UT_Vector3> weight_gradients(0, 64);
			int p_gridx = 0;
			int p_gridy = 0;
			int p_gridz = 0;
			g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position

			for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
				//Z-dimension interpolation
				float z_pos = p_gridz-z,
					wz = SIM_SnowSolver::bspline(z_pos),
					dz = SIM_SnowSolver::bsplineSlope(z_pos);

				for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
					//Y-dimension interpolation
					float y_pos = p_gridy-y,
						wy = SIM_SnowSolver::bspline(y_pos),
						dy = SIM_SnowSolver::bsplineSlope(y_pos);
					for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
						//X-dimension interpolation
						float x_pos = p_gridx-x,
							wx = SIM_SnowSolver::bspline(x_pos),
							dx = SIM_SnowSolver::bsplineSlope(x_pos);
						
						//Final weight is dyadic product of weights in each dimension
						float weight = wx*wy*wz;
						p_w[it.getOffset()-1][idx] = weight;

						//Weight gradient is a vector of partial derivatives
						const UT_Vector3 newWeight(UT_Vector3(dx*wy*wz, wx*dy*wz, wx*wy*dz));		 
						p_wgh[it.getOffset()-1][idx] = newWeight/division_size;//Always use cubed voxels??

						//Interpolate mass
					    float node_mass = g_mass->getValue(x,y,z);
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

			float density(p_density.get(it.getOffset()));				
			int p_gridx = 0;
			int p_gridy = 0;
			int p_gridz = 0;
			g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position
		
			for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
				for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
					for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
						float w = p_w[it.getOffset()-1][idx];
						if (w > BSPLINE_EPSILON){
							density += w * (g_mass->getValue(x,y,z));
						}
					}
				}
			}
			
			density /= voxelArea;
			p_density.set(it.getOffset(),density);		
			p_volume.set(it.getOffset(), particle_mass / density);
		}
	}

	/// STEP #3: Transfer velocity to grid
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
	{
		const UT_Vector3 pos(p_position.get(it.getOffset())); //Particle position pos
		const UT_Vector3 vel(p_vel.get(it.getOffset())); //Particle velocity vel
		
		int p_gridx = 0;
		int p_gridy = 0;
		int p_gridz = 0;
		g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position
	
		for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
			for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
				for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
					float w = p_w[it.getOffset()-1][idx];
					if (w > BSPLINE_EPSILON){
						float massVal = w*particle_mass;
						float nodex_vel = g_ovelX->getValue(x,y,z) + (vel[0]*massVal);
					    float nodey_vel = g_ovelY->getValue(x,y,z) + (vel[1]*massVal);
					    float nodez_vel = g_ovelZ->getValue(x,y,z) + (vel[2]*massVal);
						g_ovelX->setValue(x,y,z,nodex_vel);
						g_ovelY->setValue(x,y,z,nodey_vel);
						g_ovelZ->setValue(x,y,z,nodez_vel);			
						g_active->setValue(x,y,z,1.0);			
					}
				}
			}
		}
	}

	for(int iX=0; iX < fieldDims[0]; iX++){
		for(int iY=0; iY < fieldDims[1]; iY++){
			for(int iZ=0; iZ < fieldDims[2]; iZ++){
				float node_active = g_active->getValue(iX,iY,iZ);
				if(node_active == 1.0){
					float node_mass = g_mass->getValue(iX,iY,iZ);
					g_ovelX->setValue(iX,iY,iZ,(g_ovelX->getValue(iX,iY,iZ)/node_mass));
					g_ovelY->setValue(iX,iY,iZ,(g_ovelY->getValue(iX,iY,iZ)/node_mass));
					g_ovelZ->setValue(iX,iY,iZ,(g_ovelZ->getValue(iX,iY,iZ)/node_mass));
				}
			}
		}
	}
	//Collision detection!! ! ! ! !

	/// STEP #4: Compute new grid velocities

	/*
	int num = 0;
	UT_VectorF svd_e(1,3);//Plasticity stuff!?!?!
	UT_MatrixF svd_w(1,3,1,3);//Plasticity stuff!?!?!
	UT_MatrixF svd_v(1,3,1,3);//Plasticity stuff!?!?!
	UT_Matrix3 energy;
	UT_MatrixSolverF solver;
	svd_w.makeIdentity();
	svd_v.makeIdentity();
	svd_e(0) = 1;
	svd_e(1) = 1;
	svd_e(2) = 1;
	svd_w(1,2) = .5;
	svd_w(2,3) = 1.7;
	svd_w(3,1) = 3.8;
	
	if(!solver.SVDDecomp(svd_w,svd_e,svd_v, 100)){
		cout << "NOOOO" << endl;
	}
	cout << svd_e << endl;
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()) //Iterate through particles
	{
		//Apply plasticity here!
		const UT_Vector3 pos(p_position.get(it.getOffset())); //Particle position pos
		const UT_Vector3 vel(p_vel.get(it.getOffset())); //Particle velocity vel
		//const UT_Matrix3 energy(p_energy.get(it.getOffset()));
		const float volume(p_volume.get(it.getOffset()));
		
		int p_gridx = 0;
		int p_gridy = 0;
		int p_gridz = 0;
		g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz); //Get grid position

		//Calculating energy!!!!
		const UT_Matrix3 def_plastic(p_Fp.get(it.getOffset()));
		const UT_Matrix3 def_elastic(p_Fe.get(it.getOffset()));
		
		//SVD		
		svd_w.setSubmatrix3(1, 1, def_elastic);
		if(!solver.SVDDecomp(svd_w,svd_e,svd_v)){
			svd_w.makeIdentity();
			svd_v.makeIdentity();
			svd_e(0) = 1;
			svd_e(1) = 1;
			svd_e(2) = 1;
			cout << "Warning - check step 4" << endl;
		}
		for (int i=0; i<3; i++){
			if (svd_e(i) < CRIT_COMPRESS)
				svd_e(i) = CRIT_COMPRESS;
			else if (svd_e(i) > CRIT_STRETCH)
				svd_e(i) = CRIT_STRETCH;
		}
		
		
		

		///////////
		float harden = exp(HARDENING*(1-def_plastic.determinant()));
		float Je = (svd_e[0])*(svd_e[1])*(svd_e[2]);
	
		UT_Matrix3 svd_mult = (svd_w*svd_v);
		svd_mult.transpose();
		def_elastic.transpose();
		temp = 2*mu*(def_elastic - (svd_mult*def_elastic));

		float addToDiag = lambda*Je*(Je-1);
		temp[0] += addToDiag;
		temp[4] += addToDiag;
		temp[8] += addToDiag;
		UT_Matrix3 energy = volume * harden * temp;
		///
		
		for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
			for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
				for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
					float w = p_w[it.getOffset()-1][idx];
					if (w > BSPLINE_EPSILON){
						UT_Vector3 ngrad = p_wgh[it.getOffset()-1][idx];
						g_nvelX->setValue(x,y,z,g_nvelX->getValue(x,y,z) + ngrad.dot(energy[0]));
						g_nvelY->setValue(x,y,z,g_nvelY->getValue(x,y,z) + ngrad.dot(energy[1]));
						g_nvelZ->setValue(x,y,z,g_nvelZ->getValue(x,y,z) + ngrad.dot(energy[2]));						
					}
				}
			}
		}
		num++;
	}
	
	int node = 0;
	for(int iX=0; iX < fieldDims[0]; iX++){
		for(int iY=0; iY < fieldDims[1]; iY++){
			for(int iZ=0; iZ < fieldDims[2]; iZ++){
				
				float node_active = g_active->getValue(iX,iY,iZ);
				if(node_active){
					float nodex_ovel = g_ovelX->getValue(iX,iY,iZ);
					float nodey_ovel = g_ovelY->getValue(iX,iY,iZ);
					float nodez_ovel = g_ovelZ->getValue(iX,iY,iZ);
					float nodex_nvel = g_nvelX->getValue(iX,iY,iZ);
					float nodey_nvel = g_nvelY->getValue(iX,iY,iZ);
					float nodez_nvel = g_nvelZ->getValue(iX,iY,iZ);
					float node_mass = g_mass->getValue(iX,iY,iZ);

					nodex_ovel += timestep*(GRAVITY[0] - (nodex_nvel/node_mass)); 
					nodey_ovel += timestep*(GRAVITY[1] - (nodey_nvel/node_mass)); 
					nodez_ovel += timestep*(GRAVITY[2] - (nodez_nvel/node_mass)); 
					
					g_nvelX->setValue(iX,iY,iZ,nodex_ovel);
					g_nvelY->setValue(iX,iY,iZ,nodey_ovel);
					g_nvelZ->setValue(iX,iY,iZ,nodez_ovel);
				}
				node++;
			}
		}
	}
	*/
	/// STEP #5: Grid collision resolution
	
	/// STEP #6: Transfer grid velocities to particles and integrate
	
	/// STEP #7: Update particle deformation gradient
	
	/// STEP #8: Particle collision resolution

	//gdh.unlock(gdp_out);
    //gdh.unlock(gdp_in);
	
	return true;
}






