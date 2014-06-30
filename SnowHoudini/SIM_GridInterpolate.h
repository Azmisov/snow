#ifndef __SIM_GRIDINTERPOLATE_H__
#define __SIM_GRIDINTERPOLATE_H__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
 
class SIM_GridInterpolate : public GAS_SubSolver{
public:
	//These create setter/getter functions to access input params
	GET_DATA_FUNC_S("particles", InputParticles);
	GET_DATA_FUNC_S("inputAttr", InputAttr);
	GET_DATA_FUNC_S("ouputAttr", OutputAttr);

protected:
	//Constructor
    explicit SIM_GridInterpolate(const SIM_DataFactory *factory);
    virtual ~SIM_GridInterpolate();

	//Runs our sub-solver for a single object
	virtual bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep);
    
private:
    //Allows this solver to be a DataFactory?
	//Not sure why it needs to be a DataFactory...
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(
		SIM_GridInterpolate,		//class
		GAS_SubSolver,				//super class
		"Grid Interpolate",			//description
		getDescription()			//dop parameters
	);
    //Description of our sub-solver
    static const SIM_DopDescription *getDescription();
};

#endif
