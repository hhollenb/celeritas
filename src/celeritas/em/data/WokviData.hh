//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/data/WokviData.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Macros.hh"
#include "corecel/Types.hh"
#include "corecel/data/Collection.hh"
#include "celeritas/Quantities.hh"
#include "celeritas/Types.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Particle and action ids used by WokviModel
 */
struct WokviIds
{
    ActionId action;
    ParticleId electron;
    ParticleId positron;

    explicit CELER_FUNCTION operator bool() const
    {
        return action && electron && positron;
    }
};

//---------------------------------------------------------------------------//
/*!
 * Per-element data used by the WokviModel
 */
struct WokviElementData
{
    // Matrix of Mott coefficients
    real_type mott_coeff[5][6];
};

//---------------------------------------------------------------------------//
/*!
 * Supported models of nuclear form factors
 */
enum class NuclearFormFactorType
{
    None,
    Flat,
    Exponential,
    Gaussian
};

//---------------------------------------------------------------------------//
/*!
 * Constant shared data used by the WokviModel
 */
template<Ownership W, MemSpace M>
struct WokviData
{
    using Mass = units::MevMass;
    using MomentumSq = units::MevMomentumSq;

    template<class T>
    using ElementItems = celeritas::Collection<T, W, M, ElementId>;

    // Ids
    WokviIds ids;

    // Per element form factors
    ElementItems<WokviElementData> elem_data;

    // Mass of the electron
    Mass electron_mass;

    // Squared screening radius for incident electrons
    MomentumSq screen_r_sq_elec;

    // Nuclear form factor momentum scale
    MomentumSq form_momentum_scale;

    // Model for the form factor to use
    NuclearFormFactorType form_factor_type;

    // Check if the data is initialized
    explicit CELER_FUNCTION operator bool() const
    {
        return ids && !elem_data.empty();
    }

    // Copy initialize from an existing WokviData
    template<Ownership W2, MemSpace M2>
    WokviData& operator=(WokviData<W2, M2> const& other)
    {
        CELER_EXPECT(other);
        ids = other.ids;
        elem_data = other.elem_data;
        electron_mass = other.electron_mass;
        form_factor_type = other.form_factor_type;
        screen_r_sq_elec = other.screen_r_sq_elec;
        form_momentum_scale = other.form_momentum_scale;
        return *this;
    }
};

using WokviDeviceRef = DeviceCRef<WokviData>;
using WokviHostRef = HostCRef<WokviData>;
using WokviRef = NativeCRef<WokviData>;

//---------------------------------------------------------------------------//
}  // namespace celeritas
