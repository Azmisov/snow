#ifndef __SIM_SNOWSOLVER_H__
#define __SIM_SNOWSOLVER_H__

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>
 
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
};

#endif
