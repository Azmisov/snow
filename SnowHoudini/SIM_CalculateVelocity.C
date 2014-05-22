/*
 * Copyright (c) 2013
 *  Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 */

#include "SIM_CalculateVelocity.h"
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>

#include <iostream>
#include <stdio.h>
#include <string>

///
/// This is the hook that Houdini grabs from the dll to link in
/// this.  As such, it merely has to implement the data factory
/// for this node.
///
void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(SIM_CalculateVelocity);
}

/// Standard constructor, note that BaseClass was crated by the
/// DECLARE_DATAFACTORY and provides an easy way to chain through
/// the class hierarchy.
SIM_CalculateVelocity::SIM_CalculateVelocity(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

SIM_CalculateVelocity::~SIM_CalculateVelocity()
{
}

/// Used to automatically populate the node which will represent
/// this data type.
const SIM_DopDescription *
SIM_CalculateVelocity::getDopDescription()
{
    static PRM_Name theVelFieldName(GAS_NAME_FIELDVELOCITY, "Velocity Field");

    static PRM_Template      theTemplates[] = {
    PRM_Template(PRM_STRING, 1, &theVelFieldName),
    PRM_Template()
    };

    static SIM_DopDescription    theDopDescription(
        true,       // Should we make a DOP?
        "hdk_calcvelocity",    // Internal name of the DOP.
        "Calculate Velocity",       // Label of the DOP
        "Solver",       // Default data name
        classname(),    // The type of this DOP, usually the class.
        theTemplates);  // Template list for generating the DOP

    return &theDopDescription;
}

bool
SIM_CalculateVelocity::solveGasSubclass(SIM_Engine &engine,
            SIM_Object *obj,
            SIM_Time time,
            SIM_Time timestep)
{
    SIM_VectorField *vel_field;

    SIM_DataArray vel_data;
    getMatchingData(vel_data, obj, GAS_NAME_FIELDVELOCITY);
    
    vel_field = SIM_DATA_CAST(vel_data[0], SIM_VectorField);

    UT_Vector3 divisions = vel_field->getDivisions();
    int XDIV = divisions[0];
    int YDIV = divisions[1];
    int ZDIV = divisions[2];

    UT_VoxelArrayF* velX = vel_field->getField(0)->fieldNC();
    UT_VoxelArrayF* velY = vel_field->getField(1)->fieldNC();
    UT_VoxelArrayF* velZ = vel_field->getField(2)->fieldNC();
    
    for(int x=1; x<XDIV-1; x++)
    {
        for(int y=1; y<YDIV-1; y++)
        {
            for(int z=1; z<ZDIV-1; z++)
            {

                double x_vel = velX->getValue(x,y,z);
                double y_vel = velY->getValue(x,y,z);
                double z_vel = velZ->getValue(x,y,z);
                //TODO: calculate stuff
                velX->setValue(x,y,z, x_vel);
                velY->setValue(x,y,z, y_vel);
                velZ->setValue(x,y,z, z_vel);
            }
        }
    }

    return true;
}

void
SIM_CalculateVelocity::calculateVelPartial(SIM_RawField *dst, const SIM_RawField *src, const UT_JobInfo &info)
{

}