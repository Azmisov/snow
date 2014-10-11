#ifndef __SIM_SNOWSOLVER_H__
#define __SIM_SNOWSOLVER_H__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
#include <math.h>

#define MPM_PARTICLES "particles"
#define MPM_P_FE "p_fe"
#define MPM_P_FP "p_fp"
#define MPM_P_VEL "p_vel"
#define MPM_P_VOL "p_vol"
#define MPM_P_D "p_d"
#define MPM_P_W "p_w"
#define MPM_P_WG "p_wg"

#define MPM_G_MASS "g_mass"
#define MPM_G_NVEL "g_nvel"
#define MPM_G_OVEL "g_ovel"
#define MPM_G_ACTIVE "g_active"
#define MPM_G_DENSITY "g_density"
#define MPM_G_COL "g_col"
#define MPM_G_COLVEL "g_colVel"
#define MPM_G_EXTFORCE "g_extForce"

#define MPM_P_MASS "p_mass"
#define MPM_YOUNGS_MODULUS "youngs_modulus"
#define MPM_POISSONS_RATIO "poissons_ratio"
#define MPM_CRIT_COMP "crit_comp"
#define MPM_CRIT_STRETCH "crit_stretch"
#define MPM_FLIP_PERCENT "flip_percent"
#define MPM_HARDENING "hardening"
#define MPM_CFL "cfl"
#define MPM_COF "cof"
#define MPM_DIV_SIZE "div_size"
#define MPM_MAX_TIMESTEP "max_timestep"
#define MPM_MAX_VEL "max_vel"

#define MPM_GRAVITY "gravity"
#define MPM_BBOX_MIN "bbox_min"
#define MPM_BBOX_MAX "bbox_max"

typedef double freal;
typedef UT_Vector3T<freal> vector3;
typedef UT_Matrix3T<freal> matrix3;

static const int BSPLINE_RADIUS = 2;			//Radius of B-spline interpolation function
static const freal EPSILON = 1e-10;

class SIM_SnowSolver : public GAS_SubSolver{
public:
	GET_DATA_FUNC_S(MPM_PARTICLES, Particles);
	GET_DATA_FUNC_S(MPM_P_FE, PFe);
	GET_DATA_FUNC_S(MPM_P_FP, PFp);
	GET_DATA_FUNC_S(MPM_P_VEL, PVel);
	GET_DATA_FUNC_S(MPM_P_VOL, PVol);
	GET_DATA_FUNC_S(MPM_P_D, PD);
	GET_DATA_FUNC_S(MPM_P_W, PW);
	GET_DATA_FUNC_S(MPM_P_WG, PWg);

	GET_DATA_FUNC_S(MPM_G_MASS, GMass);
	GET_DATA_FUNC_S(MPM_G_NVEL, GNvel);
	GET_DATA_FUNC_S(MPM_G_OVEL, GOvel);
	GET_DATA_FUNC_S(MPM_G_ACTIVE, GActive);
	GET_DATA_FUNC_S(MPM_G_DENSITY, GDensity);
	GET_DATA_FUNC_S(MPM_G_COL, GCol);
	GET_DATA_FUNC_S(MPM_G_COLVEL, GColVel);
	GET_DATA_FUNC_S(MPM_G_EXTFORCE, GExtForce);

    GETSET_DATA_FUNCS_F(MPM_P_MASS, PMass);
    GETSET_DATA_FUNCS_F(MPM_YOUNGS_MODULUS, YoungsModulus);
	GETSET_DATA_FUNCS_F(MPM_POISSONS_RATIO, PoissonsRatio);
	GETSET_DATA_FUNCS_F(MPM_CRIT_COMP, CritComp);
	GETSET_DATA_FUNCS_F(MPM_CRIT_STRETCH, CritStretch);
	GETSET_DATA_FUNCS_F(MPM_FLIP_PERCENT, FlipPercent);
	GETSET_DATA_FUNCS_F(MPM_HARDENING, Hardening);
	GETSET_DATA_FUNCS_F(MPM_CFL, Cfl);
	GETSET_DATA_FUNCS_F(MPM_COF, Cof);
	GETSET_DATA_FUNCS_F(MPM_DIV_SIZE, DivSize);
	GETSET_DATA_FUNCS_F(MPM_MAX_VEL, MaxVel);
	GETSET_DATA_FUNCS_F(MPM_MAX_TIMESTEP, MaxTimestep);

	GET_DATA_FUNC_V3(MPM_GRAVITY, Gravity);
	GET_DATA_FUNC_V3(MPM_BBOX_MIN, BboxMin);
	GET_DATA_FUNC_V3(MPM_BBOX_MAX, BboxMax);
protected:
	//Constructor
    explicit SIM_SnowSolver(const SIM_DataFactory *factory);
    virtual ~SIM_SnowSolver();

	//Runs our sub-solver for a single object
	virtual bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep);
    
private:
    //Allows this solver to be a DataFactory?
	//Not sure why it needs to be a DataFactory...
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(
		SIM_SnowSolver,				//class
		GAS_SubSolver,				//super class
		"Snow Solver",				//description
		getDescription()			//dop parameters
	);
    //Description of our sub-solver
    static const SIM_DopDescription *getDescription();
    
	static freal bspline(freal x){
		x = fabs(x);
		freal w;
		if (x < 1)
			w = x*x*(x/2 - 1) + 2/3.0;
		else if (x < 2)
			w = x*(x*(-x/6 + 1) - 2) + 4/3.0;
		else return 0;
		//Clamp between 0 and 1... if needed
		if (w < EPSILON) return 0;
		return w;
	}
	static freal bsplineSlope(freal x){
		freal abs_x = fabs(x);
		if (abs_x < 1)
			return 1.5*x*abs_x - 2*x;
		else if (x < 2)
			return -x*abs_x/2 + 2*x - 2*x/abs_x;
		else return 0;
		//Clamp between -2/3 and 2/3... if needed
	}
};

#endif
