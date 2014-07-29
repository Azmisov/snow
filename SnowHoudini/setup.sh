#!/bin/bash

HFS=/opt/hfs*
if [ -z "${HFS}" ]
then
	export HFS=$HFS
fi

pushd ${HFS}
source ./houdini_setup
popd

#compile each houdini node
#hcustom SIM_CalculateVelocity.C
#hcustom SIM_GridInterpolate.c
hcustom SIM_SnowSolver.c
