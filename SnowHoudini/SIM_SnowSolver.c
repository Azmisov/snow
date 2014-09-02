#include "SIM_SnowSolver.h"
#include "Eigen/Dense"

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
#include <SIM/SIM_Engine.h>
#include <OP/OP_Node.h>
#include <OP/OP_Context.h>
#include <CH/CH_Manager.h>
#include <boost/array.hpp>

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <ctime>

typedef Eigen::Matrix<freal,3,3> eigen_matrix3;
typedef Eigen::Matrix<freal,3,1> eigen_vector3;

inline bool computeSDFNormal(const UT_VoxelArrayF *g_col, int iX, int iY, int iZ, vector3 &norm);

//Houdini hook
void initializeSIM(void *){
	IMPLEMENT_DATAFACTORY(SIM_SnowSolver);
}

//Constructor
SIM_SnowSolver::SIM_SnowSolver(const SIM_DataFactory *factory) : BaseClass(factory){}
SIM_SnowSolver::~SIM_SnowSolver(){}

//Gets node description data
const SIM_DopDescription* SIM_SnowSolver::getDescription(){
	//TODO: maybe move some of the particles attributes into SIM data fields instead of having them attached to particles?
	
	/* Particle parameters (lagrangian):
		We will just use the default point attribute "P" for particle position;
		
		Also, we can just use temporary matrices for SVD/Polar Decomp/Velocity Grad,
		since they get recomputed at each timestep anyways. If we split this up into
		multiple sub-solver nodes, we'd need to make them particle attributes.
	*/
	static PRM_Name p_field("particles", "Particles");			//particles
	static PRM_Name p_fe("p_fe", "Fe Attr");					//particle elastic deformation gradient
	static PRM_Name p_fp("p_fp", "Fp Attr");					//particle plastic deformation gradient
	static PRM_Name p_vel("p_vel", "Velocity Attr");			//particle velocity
	static PRM_Name p_vol("p_vol", "Volume Attr");				//particle volume
	static PRM_Name p_d("p_d", "Density Attr");					//particle density
	//It may be better to remove these, if recomputing these values is actually faster than caching
	static PRM_Name p_w("p_w", "Weights Attr");					//particle weight (for each node within 2-node radius)
	static PRM_Name p_wg("p_wg", "Weight Gradients Attr");		//particle weight gradient (for each node within 2-node radius)

	//Grid parameters (eulerian):
	static PRM_Name g_mass("g_mass", "Mass Field");				//grid mass
	static PRM_Name g_nvel("g_nvel", "New Velocity Field");		//grid velocity (after applying forces)
	static PRM_Name g_ovel("g_ovel", "Old Velocity Field");		//grid velocity (before applying forces)
	static PRM_Name g_active("g_active", "Activated Field");	//boolean field that tells whether there are particles within a radius of 2
	static PRM_Name g_density("g_density", "Density Field");	//grid density
	static PRM_Name g_col("g_col", "Collision Field"); 			// grid collision
	static PRM_Name g_colVel("g_colVel", "Collision Velocity Field"); 			// grid collision velocity

	static PRM_Template theTemplates[] = {
		//particles
		PRM_Template(PRM_STRING, 1, &p_field),
		PRM_Template(PRM_STRING, 1, &p_fe),
		PRM_Template(PRM_STRING, 1, &p_fp),
		PRM_Template(PRM_STRING, 1, &p_vel),
		PRM_Template(PRM_STRING, 1, &p_vol),
		PRM_Template(PRM_STRING, 1, &p_d),
		PRM_Template(PRM_STRING, 1, &p_w),
		PRM_Template(PRM_STRING, 1, &p_wg),
		//grid
		PRM_Template(PRM_STRING, 1, &g_mass),
		PRM_Template(PRM_STRING, 1, &g_nvel),
		PRM_Template(PRM_STRING, 1, &g_ovel),
		PRM_Template(PRM_STRING, 1, &g_active),
		PRM_Template(PRM_STRING, 1, &g_density),
		PRM_Template(PRM_STRING, 1, &g_col),
		PRM_Template(PRM_STRING, 1, &g_colVel),
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
bool SIM_SnowSolver::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time framerate){
	printf("Solving %.3f (%.2d), 00%%", (double) time, engine.getSimulationFrame(time));
	clock_t timer = clock();

	/// STEP #0: Retrieve all data objects from Houdini

	OP_Context context(CHgetEvalTime());
	OP_Node* solver_node = engine.getOwner()->findNode("/obj/DOP/mpm_solver/solver");
	OP_Node* collision_node = engine.getOwner()->findNode("/obj/DOP/mpm_solver/collision_mask");
	OP_Node* mpm_node = engine.getOwner()->findNode("/obj/DOP/mpm_solver");
	//Scalar params
	freal particle_mass = mpm_node->evalFloat("p_mass",0,time),
		YOUNGS_MODULUS = mpm_node->evalFloat("youngs_modulus",0,time),
		POISSONS_RATIO = mpm_node->evalFloat("poissons_ratio",0,time),
		CRIT_COMPRESS = mpm_node->evalFloat("crit_comp",0,time),
		CRIT_STRETCH = mpm_node->evalFloat("crit_stretch",0,time),
		FLIP_PERCENT = mpm_node->evalFloat("flip_percent",0,time),
		HARDENING = mpm_node->evalFloat("hardening",0,time),
		CFL = mpm_node->evalFloat("cfl",0,time),
		COF = mpm_node->evalFloat("cof",0,time),
		division_size = mpm_node->evalFloat("div_size",0,time),
		MAX_TIMESTEP = mpm_node->evalFloat("max_timestep",0,time);	
	//Vector param
	vector3 GRAVITY(
		mpm_node->evalFloat("gravity",0,time),
		mpm_node->evalFloat("gravity",1,time),
		mpm_node->evalFloat("gravity",2,time)
	);
	vector3 bbox_min_limit(
		mpm_node->evalFloat("bbox_min",0,time),
		mpm_node->evalFloat("bbox_min",1,time),
		mpm_node->evalFloat("bbox_min",2,time)
	);
	vector3 bbox_max_limit(
		mpm_node->evalFloat("bbox_max",0,time),
		mpm_node->evalFloat("bbox_max",1,time),
		mpm_node->evalFloat("bbox_max",2,time)
	);

	UT_String s_p, s_vol, s_den, s_vel, s_fe, s_fp;
	solver_node->evalString(s_p, "particles",0,time);
	solver_node->evalString(s_vol,"p_vol",0,time);
	solver_node->evalString(s_den,"p_d",0,time);
	solver_node->evalString(s_vel,"p_vel",0,time);
	solver_node->evalString(s_fe,"p_fe",0,time);
	solver_node->evalString(s_fp,"p_fp",0,time);

	SIM_Geometry* geometry = (SIM_Geometry*) obj->getNamedSubData(s_p);
	if (!geometry) return true;
	
	//Get particle data
	//Do we use the attribute name???
	// GU_DetailHandle gdh = geometry->getGeometry().getWriteableCopy();
	GU_DetailHandle gdh = geometry->getOwnGeometry();
	const GU_Detail* gdp_in = gdh.readLock(); // Must unlock later
	GU_Detail* gdp_out = gdh.writeLock();

	GA_RWAttributeRef p_ref_position = gdp_out->findPointAttribute("P");
	GA_RWHandleT<vector3> p_position(p_ref_position.getAttribute());

	GA_RWAttributeRef p_ref_volume = gdp_out->findPointAttribute(s_vol);
	GA_RWHandleT<freal> p_volume(p_ref_volume.getAttribute());

	GA_RWAttributeRef p_ref_density = gdp_out->findPointAttribute(s_den);
	GA_RWHandleT<freal> p_density(p_ref_density.getAttribute());

	GA_RWAttributeRef p_ref_vel = gdp_out->findPointAttribute(s_vel);
	GA_RWHandleT<vector3> p_vel(p_ref_vel.getAttribute());

	GA_RWAttributeRef p_ref_Fe = gdp_out->findPointAttribute(s_fe);
	GA_RWHandleT<matrix3> p_Fe(p_ref_Fe.getAttribute());

	GA_RWAttributeRef p_ref_Fp = gdp_out->findPointAttribute(s_fp);
	GA_RWHandleT<matrix3> p_Fp(p_ref_Fp.getAttribute());

	//EVALUATE PARAMETERS
	freal mu = YOUNGS_MODULUS/(2+2*POISSONS_RATIO);
	freal lambda = YOUNGS_MODULUS*POISSONS_RATIO/((1+POISSONS_RATIO)*(1-2*POISSONS_RATIO));
	freal timestep;

	//Get grid data
	SIM_ScalarField *g_mass_field;
	SIM_DataArray g_mass_data;
	getMatchingData(g_mass_data, obj, "g_mass");	
	g_mass_field = SIM_DATA_CAST(g_mass_data(0), SIM_ScalarField);

	SIM_VectorField *g_nvel_field;
	SIM_DataArray g_nvel_data;
	getMatchingData(g_nvel_data, obj, "g_nvel");
	g_nvel_field = SIM_DATA_CAST(g_nvel_data(0), SIM_VectorField);

	SIM_VectorField *g_ovel_field;
	SIM_DataArray g_ovel_data;
	getMatchingData(g_ovel_data, obj, "g_ovel");
	g_ovel_field = SIM_DATA_CAST(g_ovel_data(0), SIM_VectorField);

	SIM_ScalarField *g_active_field;
	SIM_DataArray g_active_data;
	getMatchingData(g_active_data, obj, "g_active");	
	g_active_field = SIM_DATA_CAST(g_active_data(0), SIM_ScalarField);

	SIM_ScalarField *g_density_field;
	SIM_DataArray g_density_data;
	getMatchingData(g_density_data, obj, "g_density");	
	g_density_field = SIM_DATA_CAST(g_density_data(0), SIM_ScalarField);

	SIM_ScalarField *g_col_field;
	SIM_DataArray g_col_data;
	getMatchingData(g_col_data, obj, "g_col");	
	g_col_field = SIM_DATA_CAST(g_col_data(0), SIM_ScalarField);

	SIM_VectorField *g_colVel_field;
	SIM_DataArray g_colVel_data;
	getMatchingData(g_colVel_data, obj, "g_colVel");	
	g_colVel_field = SIM_DATA_CAST(g_colVel_data(0), SIM_VectorField);
	
	UT_VoxelArrayF *g_nvelX, *g_nvelY, *g_nvelZ, *g_ovelX, *g_ovelY, *g_ovelZ, *g_active, *g_mass, *g_density,
			*g_col, *g_colVelX, *g_colVelY, *g_colVelZ;

	int point_count = gdp_out->getPointRange().getEntries();
	std::vector<boost::array<freal,64> > p_w(point_count);
	std::vector<boost::array<vector3,64> > p_wgh(point_count);

	/*
	//Find initial maximum velocity, for adaptive timestep
	//Also find initial bounding box of particles
	vector3 bbox_min, bbox_max;
	bool bbox_reset = true;
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()){
		int pid = it.getOffset();
		//Adjust bbox
		vector3 pos(p_position.get(pid));
		if (bbox_reset){
			bbox_min = pos;
			bbox_max = pos;
			bbox_reset = false;
		}
		else{
			for (int i=0; i<3; i++){
				if (bbox_min[i] > pos[i])
					bbox_min[i] = pos[i];
				else if (bbox_max[i] < pos[i])
					bbox_max[i] = pos[i];
			}
		}
	}
	if (bbox_reset)
		cout << "This should never happen!!!!!!!" << endl;
	


	//Get new grid dimensions
	grid_origin = (bbox_max+bbox_min)/2.0;	//grid center (to be adjusted to grid origin later)
	grid_size = bbox_max-bbox_min; //grid size (to be adjusted to grid scale later)
	//this gives us 2+ voxel padding, approximately; there should be two border layers
	//for the bspline weight function; and an additional border layer to compute collision normals
	grid_size += division_size*7;
	//*/

	//Grid dimensions, for resizing
	vector3 grid_origin, grid_divs;
	freal voxelArea = division_size*division_size*division_size;

	g_mass = g_mass_field->getField()->fieldNC();
	g_nvelX = g_nvel_field->getField(0)->fieldNC();
	g_nvelY = g_nvel_field->getField(1)->fieldNC();
	g_nvelZ = g_nvel_field->getField(2)->fieldNC();
	g_ovelX = g_ovel_field->getField(0)->fieldNC();
	g_ovelY = g_ovel_field->getField(1)->fieldNC();
	g_ovelZ = g_ovel_field->getField(2)->fieldNC();
	
	g_colVelX = g_colVel_field->getField(0)->fieldNC();
	g_colVelY = g_colVel_field->getField(1)->fieldNC();
	g_colVelZ = g_colVel_field->getField(2)->fieldNC();
	g_col = g_col_field->getField()->fieldNC();
	g_active = g_active_field->getField()->fieldNC();
	grid_divs = g_mass_field->getDivisions();

	//TODO: only resize grids if the bounding box changes?
#if 0
	if (0){ //!mapping_density){
		//Substep progress
		/*
		int progress = (int) (100*adaptive_time/framerate);
		if (progress < 100){
			printf("\b\b\b%02i%%", progress);
			fflush(stdout);
		}//*/

		//Resize grid
		g_mass_field->resizeKeepData(grid_size, grid_origin, false);
		g_nvel_field->resizeKeepData(grid_size, grid_origin, false);
		g_ovel_field->resizeKeepData(grid_size, grid_origin, false);
		g_active_field->resizeKeepData(grid_size, grid_origin, false);
		
		g_col_field->resizeKeepData(grid_size, grid_origin, true);
		g_colVel_field->resizeKeepData(grid_size, grid_origin, true);
		bbox_reset = true;
	
		//Pointers may be invalid after resize
		g_mass = g_mass_field->getField()->fieldNC();
		g_nvelX = g_nvel_field->getField(0)->fieldNC();
		g_nvelY = g_nvel_field->getField(1)->fieldNC();
		g_nvelZ = g_nvel_field->getField(2)->fieldNC();
		g_ovelX = g_ovel_field->getField(0)->fieldNC();
		g_ovelY = g_ovel_field->getField(1)->fieldNC();
		g_ovelZ = g_ovel_field->getField(2)->fieldNC();
		
		g_colVelX = g_colVel_field->getField(0)->fieldNC();
		g_colVelY = g_colVel_field->getField(1)->fieldNC();
		g_colVelZ = g_colVel_field->getField(2)->fieldNC();
		g_col = g_col_field->getField()->fieldNC();
		g_active = g_active_field->getField()->fieldNC();
		grid_divs = g_mass_field->getDivisions();

		//Reset grid
		for(int iX=0; iX < grid_divs[0]; iX++){
			for(int iY=0; iY < grid_divs[1]; iY++){
				for(int iZ=0; iZ < grid_divs[2]; iZ++){
					g_mass->setValue(iX,iY,iZ,0);
					g_active->setValue(iX,iY,iZ,0);
					g_ovelX->setValue(iX,iY,iZ,0);
					g_ovelY->setValue(iX,iY,iZ,0);
					g_ovelZ->setValue(iX,iY,iZ,0);
					g_nvelX->setValue(iX,iY,iZ,0);
					g_nvelY->setValue(iX,iY,iZ,0);
					g_nvelZ->setValue(iX,iY,iZ,0);
				}
			}
		}
	}
	else{
		/*
		//Resize density field
		g_density_field->resizeKeepData(grid_size, grid_origin, false);
		g_density = g_density_field->getField()->fieldNC();
		grid_divs = g_density_field->getDivisions();
		for(int iX=0; iX < grid_divs[0]; iX++){
			for(int iY=0; iY < grid_divs[1]; iY++){
				for(int iZ=0; iZ < grid_divs[2]; iZ++){
					g_density->setValue(iX,iY,iZ,0);
				}
			}
		}
		resized_field = g_density_field;//*/
	}
#endif

	//Get world-to-grid conversion ratios
	//Particle's grid position can be found via (pos - grid_origin)/grid_cellsize
	grid_origin = g_mass_field->getCenter() - g_mass_field->getSize()/2.0;

	/// STEP #1: Transfer mass to grid
	
	if (p_position.isValid()){

		//Iterate through particles
		for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()){
			int pid = it.getOffset();
							
			//Get grid position
			vector3 gpos = (p_position.get(pid) - grid_origin)/division_size;
			int p_gridx = 0, p_gridy = 0, p_gridz = 0;
			g_mass_field->posToIndex(p_position.get(pid),p_gridx,p_gridy,p_gridz);
			freal particle_density = p_density.get(pid);
			//Compute weights and transfer mass
			for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
				//Z-dimension interpolation
				freal z_pos = gpos[2]-z,
					wz = SIM_SnowSolver::bspline(z_pos),
					dz = SIM_SnowSolver::bsplineSlope(z_pos);
				for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
					//Y-dimension interpolation
					freal y_pos = gpos[1]-y,
						wy = SIM_SnowSolver::bspline(y_pos),
						dy = SIM_SnowSolver::bsplineSlope(y_pos);
					for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
						//X-dimension interpolation
						freal x_pos = gpos[0]-x,
							wx = SIM_SnowSolver::bspline(x_pos),
							dx = SIM_SnowSolver::bsplineSlope(x_pos);
						
						//Final weight is dyadic product of weights in each dimension
						freal weight = wx*wy*wz;
						p_w[pid-1][idx] = weight;

						if (1){  // !mapping_density){
							//Weight gradient is a vector of partial derivatives
							p_wgh[pid-1][idx] = vector3(dx*wy*wz, wx*dy*wz, wx*wy*dz)/division_size;

							//Interpolate mass
							freal node_mass = g_mass->getValue(x,y,z);
							node_mass += weight*particle_mass;
							g_mass->setValue(x,y,z,node_mass);
						}
						else{
							freal density = particle_density*weight+g_density->getValue(x,y,z);
							g_density->setValue(x,y,z,density);
						}
					}
				}
			}
		}
	}

	//Exit early
	//if (mapping_density) break;
	
	/// STEP #2: First timestep only - Estimate particle volumes using grid mass

	/*
	if (time == 0.0){
		//Iterate through particles
		for (GA_Iterator it(gdp_out->getPointRange()); !it.atEnd(); it.advance()){
			int pid = it.getOffset();
			freal density = 0;

			//Get grid position
			int p_gridx = 0, p_gridy = 0, p_gridz = 0;
			g_nvel_field->posToIndex(0,p_position.get(pid),p_gridx,p_gridy,p_gridz);
			//Transfer grid density (within radius) to particles
			for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
				for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
					for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
						freal w = p_w[pid-1][idx];
						if (w > EPSILON){
							//Transfer density
							density += w * g_mass->getValue(x,y,z);
						}
					}
				}
			}
			
			density /= voxelArea;
			p_density.set(pid,density);
			p_volume.set(pid, particle_mass/density);
		}
	}
	//*/
	
	/// STEP #3: Transfer velocity to grid

	//This must happen after transferring mass, to conserve momentum
	//Iterate through particles and transfer
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()){
		int pid = it.getOffset();
		vector3 vel_fac = p_vel.get(pid)*particle_mass;

		//Get grid position
		int p_gridx = 0, p_gridy = 0, p_gridz = 0;
		g_nvel_field->posToIndex(0,p_position.get(pid),p_gridx,p_gridy,p_gridz);

		//Transfer to grid nodes within radius
		for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
			for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
				for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
					freal w = p_w[pid-1][idx];
					if (w > EPSILON){
						freal nodex_vel = g_ovelX->getValue(x,y,z) + vel_fac[0]*w;
						freal nodey_vel = g_ovelY->getValue(x,y,z) + vel_fac[1]*w;
						freal nodez_vel = g_ovelZ->getValue(x,y,z) + vel_fac[2]*w;
						g_ovelX->setValue(x,y,z,nodex_vel);
						g_ovelY->setValue(x,y,z,nodey_vel);
						g_ovelZ->setValue(x,y,z,nodez_vel);			
						g_active->setValue(x,y,z,1.0);			
					}
				}
			}
		}
	}

	//Division is slow; we only want to do divide by mass once, for each active node
	for(int iX=0; iX < grid_divs[0]; iX++){
		for(int iY=0; iY < grid_divs[1]; iY++){
			for(int iZ=0; iZ < grid_divs[2]; iZ++){
				//Only check nodes that have mass
				if (g_active->getValue(iX,iY,iZ)){
					freal node_mass = 1/(g_mass->getValue(iX,iY,iZ));
					g_ovelX->setValue(iX,iY,iZ,(g_ovelX->getValue(iX,iY,iZ)*node_mass));
					g_ovelY->setValue(iX,iY,iZ,(g_ovelY->getValue(iX,iY,iZ)*node_mass));
					g_ovelZ->setValue(iX,iY,iZ,(g_ovelZ->getValue(iX,iY,iZ)*node_mass));
				}
			}
		}
	}
	
	/// STEP #4: Compute new grid velocities

	//Temporary variables for plasticity and force calculation
	//We need one set of variables for each thread that will be running
	eigen_matrix3 def_elastic, def_plastic, energy, svd_u, svd_v;
	Eigen::JacobiSVD<eigen_matrix3, Eigen::NoQRPreconditioner> svd;
	eigen_vector3 svd_e;
	matrix3  HDK_def_plastic, HDK_def_elastic, HDK_energy;
	freal* data_dp = HDK_def_plastic.data();
	freal* data_de = HDK_def_elastic.data();
	freal* data_energy = HDK_energy.data();
	//Map Eigen matrices to HDK matrices
	Eigen::Map<eigen_matrix3> data_dp_map(data_dp);
	Eigen::Map<eigen_matrix3> data_de_map(data_de);
	Eigen::Map<eigen_matrix3> data_energy_map(data_energy);	

	//Compute force at each particle and transfer to Eulerian grid
	//We use "nvel" to hold the grid force, since that variable is not in use
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()){
		int pid = it.getOffset();
		
		//Apply plasticity to deformation gradient, before computing forces
		//We need to use the Eigen lib to do the SVD; transfer houdini matrices to Eigen matrices
		HDK_def_plastic = p_Fp.get(pid);
		HDK_def_elastic = p_Fe.get(pid);
		def_plastic = Eigen::Map<eigen_matrix3>(data_dp);
		def_elastic = Eigen::Map<eigen_matrix3>(data_de);
		
		//Compute singular value decomposition (uev*)
		svd.compute(def_elastic, Eigen::ComputeFullV | Eigen::ComputeFullU);
		svd_e = svd.singularValues();
		svd_u = svd.matrixU();
		svd_v = svd.matrixV();
		//Clamp singular values
		for (int i=0; i<3; i++){
			if (svd_e[i] < CRIT_COMPRESS) 
				svd_e[i] = CRIT_COMPRESS;
			else if (svd_e[i] > CRIT_STRETCH)
				svd_e[i] = CRIT_STRETCH;
		}
		//Put SVD back together for new elastic and plastic gradients
		def_plastic = svd_v * svd_e.asDiagonal().inverse() * svd_u.transpose() * def_elastic * def_plastic;
		svd_v.transposeInPlace();
		def_elastic = svd_u * svd_e.asDiagonal() * svd_v;
		
		//Now compute the energy partial derivative (which we use to get force at each grid node)
		energy = 2*mu*(def_elastic - svd_u*svd_v)*def_elastic.transpose();
		//Je is the determinant of def_elastic (equivalent to svd_e.prod())
		freal Je = svd_e.prod(),
			contour = lambda*Je*(Je-1),
			jp = def_plastic.determinant(),
			particle_vol = p_volume.get(pid);
		for (int i=0; i<3; i++)
			energy(i,i) += contour;
		energy *=  particle_vol * exp(HARDENING*(1-jp));
		//p_density.set(pid,1/jp);
		
		//Transfer Eigen matrices back to HDK
		data_dp_map = def_plastic;
		data_de_map = def_elastic;
		data_energy_map = energy;
		
		p_Fp.set(pid,HDK_def_plastic);
		p_Fe.set(pid,HDK_def_elastic);
		
		//Transfer energy to surrounding grid nodes
		int p_gridx = 0, p_gridy = 0, p_gridz = 0;
		g_nvel_field->posToIndex(0,p_position.get(pid),p_gridx,p_gridy,p_gridz);
		for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
			for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
				for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
					freal w = p_w[pid-1][idx];
					if (w > EPSILON){
						vector3 ngrad = p_wgh[pid-1][idx];
						g_nvelX->setValue(x,y,z,g_nvelX->getValue(x,y,z) + ngrad.dot(HDK_energy[0]));
						g_nvelY->setValue(x,y,z,g_nvelY->getValue(x,y,z) + ngrad.dot(HDK_energy[1]));
						g_nvelZ->setValue(x,y,z,g_nvelZ->getValue(x,y,z) + ngrad.dot(HDK_energy[2]));						
					}
				}
			}
		}
	}

	//Use new forces to solve for new velocities
	for(int iX=0; iX < grid_divs[0]; iX++){
		for(int iY=0; iY < grid_divs[1]; iY++){
			for(int iZ=0; iZ < grid_divs[2]; iZ++){
				//Only compute for active nodes
				if (g_active->getValue(iX,iY,iZ)){
					freal nodex_ovel = g_ovelX->getValue(iX,iY,iZ);
					freal nodey_ovel = g_ovelY->getValue(iX,iY,iZ);
					freal nodez_ovel = g_ovelZ->getValue(iX,iY,iZ);
					freal forcex = g_nvelX->getValue(iX,iY,iZ);
					freal forcey = g_nvelY->getValue(iX,iY,iZ);
					freal forcez = g_nvelZ->getValue(iX,iY,iZ);
					freal node_mass = 1/(g_mass->getValue(iX,iY,iZ));

					nodex_ovel += timestep*(GRAVITY[0] - forcex*node_mass);
					nodey_ovel += timestep*(GRAVITY[1] - forcey*node_mass);
					nodez_ovel += timestep*(GRAVITY[2] - forcez*node_mass);
					
					g_nvelX->setValue(iX,iY,iZ,nodex_ovel);
					g_nvelY->setValue(iX,iY,iZ,nodey_ovel);
					g_nvelZ->setValue(iX,iY,iZ,nodez_ovel);
				}
			}
		}
	}

	/// STEP #5: Grid collision resolution

	vector3 sdf_normal;
	for(int iX=1; iX < grid_divs[0]-1; iX++){
		for(int iY=1; iY < grid_divs[1]-1; iY++){
			for(int iZ=1; iZ < grid_divs[2]-1; iZ++){
				if (g_active->getValue(iX,iY,iZ)){
					//Compute surface normal at this point; (gradient of SDF)
					if (!computeSDFNormal(g_col, iX, iY, iZ, sdf_normal))
						continue;
					//Collider velocity
					vector3 vco(
						g_colVelX->getValue(iX,iY,iZ),
						g_colVelY->getValue(iX,iY,iZ),
						g_colVelZ->getValue(iX,iY,iZ)
					);
					//Grid velocity
					vector3 v(
						g_nvelX->getValue(iX,iY,iZ),
						g_nvelY->getValue(iX,iY,iZ),
						g_nvelZ->getValue(iX,iY,iZ)
					);
					//Skip if bodies are separating
					vector3 vrel = v - vco;
					freal vn = vrel.dot(sdf_normal);
					if (vn >= 0) continue;
					//Resolve collisions; also add velocity of collision object to snow velocity
					//Sticks to surface (too slow to overcome static friction)
					vector3 vt = vrel - (sdf_normal*vn);
					freal stick = vn*COF, vt_norm = vt.length();
					if (vt_norm <= -stick)
						vt = vco;
					//Dynamic friction
					else vt += stick*vt/vt_norm + vco;
					g_nvelX->setValue(iX,iY,iZ,vt[0]);	
					g_nvelY->setValue(iX,iY,iZ,vt[1]);
					g_nvelZ->setValue(iX,iY,iZ,vt[2]);
				}
			}
		}
	}

	/// STEP #6: Transfer grid velocities to particles and integrate
	/// STEP #7: Particle collision resolution

	vector3 pic, flip, col_vel;
	matrix3 vel_grad;
	//Iterate through particles
	for (GA_Iterator it(gdp_in->getPointRange()); !it.atEnd(); it.advance()){
		int pid = it.getOffset();
		//Particle position
		vector3 pos(p_position.get(pid));
		
		//Reset velocity
		pic[0] = 0.0;
		pic[1] = 0.0;
		pic[2] = 0.0;
		flip = p_vel.get(pid);
		vel_grad.zero();
		freal density = 0;

		 //Get grid position
		int p_gridx = 0, p_gridy = 0, p_gridz = 0;
		g_nvel_field->posToIndex(0,pos,p_gridx,p_gridy,p_gridz);
		for (int idx=0, z=p_gridz-1, z_end=z+3; z<=z_end; z++){
			for (int y=p_gridy-1, y_end=y+3; y<=y_end; y++){
				for (int x=p_gridx-1, x_end=x+3; x<=x_end; x++, idx++){
					freal w = p_w[pid-1][idx];
					if (w > EPSILON){
						const vector3 node_wg = p_wgh[pid-1][idx];
						const vector3 node_nvel(
							g_nvelX->getValue(x,y,z),
							g_nvelY->getValue(x,y,z),
							g_nvelZ->getValue(x,y,z)
						);

						//Transfer velocities
						pic += node_nvel*w;	
						flip[0] += (node_nvel[0] - g_ovelX->getValue(x,y,z))*w;	
						flip[1] += (node_nvel[1]- g_ovelY->getValue(x,y,z))*w;	
						flip[2] += (node_nvel[2] - g_ovelZ->getValue(x,y,z))*w;
						//Transfer density
						density += w * g_mass->getValue(x,y,z);
						//Transfer veloctiy gradient
						vel_grad.outerproductUpdate(1.0, node_nvel, node_wg);
					}
				}
			}
		}

		//Finalize velocity update
		vector3 vel = flip*FLIP_PERCENT + pic*(1-FLIP_PERCENT);
		
		//Interpolate surrounding nodes' SDF info to the particle (trilinear interpolation)
		freal col_sdf = g_col_field->getValue(pos);
		sdf_normal = g_col_field->getGradient(pos);
		col_vel = g_colVel_field->getValue(pos);

		//Resolve particle collisions
		if (col_sdf > 0){
			//Skip if bodies are separating
			vector3 vrel = vel - col_vel;
			freal vn = vrel.dot(sdf_normal);
			if (vn < 0){
				//Resolve and add velocity of collision object to snow velocity
				//Sticks to surface (too slow to overcome static friction)
				vel = vrel - (sdf_normal*vn);
				freal stick = vn*COF, vel_norm = vel.length();
				if (vel_norm <= -stick)
					vel = col_vel;
				//Dynamic friction
				else vel += stick*vel/vel_norm + col_vel;
			}
		}

		SKIP_PCOLLIDE:

		//Finalize density update
		density /= voxelArea;
		p_density.set(pid,density);

		//Update particle position
		pos += timestep*vel;
		//Limit particle position
		int mask = 0;
		for (int i=0; i<3; i++){
			if (pos[i] > bbox_max_limit[i]){
				pos[i] = bbox_max_limit[i];
				vel[i] = 0.0;
				mask |= 1 << i;
			}
			else if (pos[i] < bbox_min_limit[i]){
				pos[i] = bbox_min_limit[i];
				vel[i] = 0.0;
				mask |= 1 << i;
			}
		}
		//Slow particle down at bounds (not really necessary...)
		if (mask){
			for (int i=0; i<3; i++){
				if (mask & 0x1)
					vel[i] *= .7;
				mask >>= 1;
			}
		}
		p_vel.set(pid,vel);
		p_position.set(pid,pos);

		//Update particle deformation gradient
		//Note: plasticity is computed on the next timestep...
		vel_grad *= timestep;
		vel_grad(0,0) += 1;
		vel_grad(1,1) += 1;
		vel_grad(2,2) += 1;
		
		p_Fe.set(pid, vel_grad*p_Fe.get(pid));
		/*
		//Update bounding box (for grid resizing)
		if (bbox_reset){
			bbox_reset = false;
			bbox_min = pos;
			bbox_max = pos;
		}
		else{
			for (int i=0; i<3; i++){
				if (bbox_min[i] > pos[i])
					bbox_min[i] = pos[i];
				else if (bbox_max[i] < pos[i])
					bbox_max[i] = pos[i];
			}
		}//*/
	}
	

	gdh.unlock(gdp_out);
    gdh.unlock(gdp_in);
	/*
    //Output simulation time
    timer = clock() - timer;
    float secs = ((float) timer)/CLOCKS_PER_SEC;
    if (secs > 60){
		int mins = secs/60;
		secs -= mins*60;
		printf("\b\b\b%dm%ds\n",mins,(int) secs);
	}
    else printf ("\b\b\b%.2fs\n",secs);*/
	
	return true;
}

inline bool computeSDFNormal(const UT_VoxelArrayF *g_col, int iX, int iY, int iZ, vector3 &norm){
	//Make sure this is a border cell?????
	//*
	if (g_col->getValue(iX, iY, iZ) <= 0)
		return false;
	//*/
	freal sdf_xb = g_col->getValue(iX-1,iY,iZ),
		sdf_xa = g_col->getValue(iX+1,iY,iZ),
		sdf_yb = g_col->getValue(iX,iY-1,iZ),
		sdf_ya = g_col->getValue(iX,iY+1,iZ),
		sdf_zb = g_col->getValue(iX,iY,iZ-1),
		sdf_za = g_col->getValue(iX,iY,iZ+1);
	/*
	if (!(sdf_xb < 1 || sdf_xa < 1 || sdf_yb < 1 || sdf_ya < 1 || sdf_zb < 1 || sdf_za < 1) ||
		!(sdf_xb || sdf_xa || sdf_yb || sdf_ya || sdf_zb || sdf_za))
		return false;
	//*/
	norm[0] = sdf_xb - sdf_xa;
	norm[1] = sdf_yb - sdf_ya;
	norm[2] = sdf_zb - sdf_za;
	norm.normalize();
	return true;
}



