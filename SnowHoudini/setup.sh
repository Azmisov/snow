if [ -z "${HFS}" ]
then
    export HFS=/opt/hfs13.0.343
fi

pushd ${HFS}
source ./houdini_setup
popd

#compile each houdini node
hcustom SIM_CalculateVelocity.C
hcustom SIM_GridInterpolate.c
hcustom SIM_SnowSolver.c
