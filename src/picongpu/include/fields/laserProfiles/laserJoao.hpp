/**
 * Copyright 2013-2016 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch, Stefan Tietze
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */



#pragma once

#include "pmacc_types.hpp"
#include "simulation_defines.hpp"

namespace picongpu
{
/** not focusing wavepaket with spacial gaussian envelope
 *
 *  no phase shifts, just spacial envelope
 *  including correction to laser formular derived from vector potential, so the integration
 *  along propagation direction gives 0
 *  this is important for few-cycle laser pulses
 */
namespace laserJoao
{

/** Compute the
 *
 */
HINLINE float3_X laserLongitudinal(uint32_t currentStep, float_X& phase)
{
    float_X envelope = float_X(AMPLITUDE);
    float3_X elong(float3_X::create(0.0));

    // a symmetric pulse will be initialized at position z=0 for
    // a time of RAMP_INIT * PULSE_LENGTH + LASER_NOFOCUS_CONSTANT = INIT_TIME.
    // we shift the complete pulse for the half of this time to start with
    // the front of the laser pulse.
    const float_64 mue = 0.5 * INIT_TIME;

    const float_64 runTime = DELTA_T*currentStep - mue;
    const float_64 f = SPEED_OF_LIGHT / WAVE_LENGTH;

    const float_64 w = 2.0 * PI * f;

    const float_64 tau = PULSE_LENGTH * sqrt(2.0);

    const float_64 exponent = (runTime / PULSE_LENGTH / sqrt(2.0));
    envelope *= math::exp(-0.5 * exponent * exponent);

    phase += float_X(w * runTime) + LASER_PHASE;

    if( Polarisation == LINEAR_X )
    {
        elong.x() = float_X(envelope * math::sin(phase));
    }
    else if( Polarisation == LINEAR_Z )
    {
        elong.z() = float_X(envelope * math::sin(phase));
    }
    else if( Polarisation == CIRCULAR )
    {
        elong.x() = float_X(envelope / sqrt(2.0) * math::sin(phase) );
        elong.z() = float_X(envelope / sqrt(2.0) * math::cos(phase) );
    }

    return elong;
}

/**
 *
 * @param elong
 * @param phase
 * @param posX
 * @param posZ
 * @return
 */
HDINLINE float3_X laserTransversal(float3_X elong, const float_X, const float_X posX, const float_X posZ)
{

    const float_64 exp_x = float_64(posX) * float_64(posX) / (W0_X * W0_X);
    const float_64 exp_z = posZ * posZ / (W0_Z * W0_Z);

    const float_64 argCosh = Q_COS * posX / (PI * W0_X* W0_X);
    const float_64 cosh =  0.5 * (math::exp( argCosh ) + math::exp( -argCosh ));
    const float_64 ExpSigma = math::exp( float_64(-1.0)*( Q_COS * Q_COS / (4 * PI * PI * W0_X * W0_X)));

    //    printf("joao laser: %f  <  %f \n", posX, 0.5 * jetWidth);
    if ( (math::abs( posX )  <  0.5 * jetWidth) ||  ( (math::abs( posX ) > 20.* W0_X)  || (math::abs( posZ ) > 20.* W0_Z) ))
    {
        return elong * 0.0;
    }
    else
    {
        const float3_64 ret = precisionCast<float_64>(elong)/(p + h * ExpSigma ) * math::exp( float_64(-1.0)*(exp_x + exp_z)) * (p + h * ExpSigma * cosh);
        //printf("JoaoLaser: E=%f  x=%f  z=%f \n", math::abs(ret), posX, posZ);
        return precisionCast<float_X>(ret);
    }
}

}
}



