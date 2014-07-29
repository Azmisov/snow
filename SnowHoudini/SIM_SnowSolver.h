#ifndef __SIM_SNOWSOLVER_H__
#define __SIM_SNOWSOLVER_H__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
#include <math.h>

static const float 
				BSPLINE_EPSILON = 1e-4,
				TWO_THIRDS = 2/3.0,
				YOUNGS_MODULUS = 2.0,		//Young's modulus (springiness) (1.4e5)
				POISSONS_RATIO = .3,		//Poisson's ratio (transverse/axial strain ratio) (.2)
				CRIT_COMPRESS = 1-2.4e-2,	//Fracture threshold for compression (1-2.5e-2)
				CRIT_STRETCH = 1+7.5e-3,	//Fracture threshold for stretching (1+7.5e-3)
				FLIP_PERCENT = .95,			//Weight to give FLIP update over PIC (.95)
				HARDENING = 5.0;		//How much plastic deformation strengthens material (10)


static const UT_Vector3 GRAVITY(0.0,-1.0,0.0);

class SIM_SnowSolver : public GAS_SubSolver{
public:
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
	static float bspline(float x){
		x = fabs(x);
		float w;
		if (x < 1)
			w = x*x*(x/2 - 1) + 2/3.0;
		else if (x < 2)
			w = x*(x*(-x/6 + 1) - 2) + 4/3.0;
		else return 0;
		//Clamp between 0 and 1... if needed
		if (w < BSPLINE_EPSILON) return 0;
		return w;
	}

	static float bsplineSlope(float x){
		float abs_x = fabs(x);
		if (abs_x < 1)
			return 1.5*x*abs_x - 2*x;
		else if (x < 2)
			return -x*abs_x/2 + 2*x - 2*x/abs_x;
		else return 0;
		//Clamp between -2/3 and 2/3... if needed
	}
};

#endif
