#ifndef __SIM_SNOWSOLVER_H__
#define __SIM_SNOWSOLVER_H__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
#include <math.h>
 
const float BSPLINE_EPSILON = 1e-4;
const float TWO_THIRDS = 2/3.0;

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
			w = x*x*(x/2 - 1) + TWO_THIRDS;
		else if (x < 2)
			w = x*(x*(-x/6 + 1) - 2) + 2*TWO_THIRDS;
		else return 0;
		//Clamp between 0 and 1... if needed
		if (w < BSPLINE_EPSILON) return 0;
		return w;
	}

	static float bsplineSlope(float x){
		float abs_x = fabs(x), w;
		if (abs_x < 1)
			return 1.5*x*abs_x - 2*x;
		else if (x < 2)
			return -x*abs_x/2 + 2*x - 2*x/abs_x;
		else return 0;
		//Clamp between -2/3 and 2/3... if needed
	}
};

#endif
