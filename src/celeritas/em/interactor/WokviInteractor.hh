//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/interactor/WokviInteractor.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Macros.hh"
#include "corecel/data/StackAllocator.hh"
#include "corecel/math/ArrayUtils.hh"
#include "celeritas/Constants.hh"
#include "celeritas/Quantities.hh"
#include "celeritas/Types.hh"
#include "celeritas/em/data/WokviData.hh"
#include "celeritas/em/distribution/WokviDistribution.hh"
#include "celeritas/em/interactor/detail/WokviStateHelper.hh"
#include "celeritas/mat/ElementView.hh"
#include "celeritas/mat/MaterialView.hh"
#include "celeritas/phys/Interaction.hh"
#include "celeritas/phys/ParticleTrackView.hh"
#include "celeritas/phys/Secondary.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Applies the Wentzel OK and VI single Coulomb scattering model.
 */
class WokviInteractor
{
  public:
    //!@{
    //! \name Type aliases
    using Energy = units::MevEnergy;
    //!@}

  public:
    // Construct with shared and state data
    inline CELER_FUNCTION WokviInteractor(WokviRef const& shared,
                                          ParticleTrackView const& particle,
                                          Real3 const& inc_direction,
                                          MaterialView const& material,
                                          ElementComponentId const& elcomp_id,
                                          StackAllocator<Secondary>& allocate);

    template<class Engine>
    inline CELER_FUNCTION Interaction operator()(Engine& rng);

  private:
    //// DATA ////

    // Constant shared data
    WokviRef const& data_;

    // Incident direction
    Real3 const& inc_direction_;

    // Allocator for secondary tracks
    StackAllocator<Secondary>& allocate_;

    real_type inc_energy_;
    real_type inc_mass_;
    IsotopeView target_;

    //// HELPER FUNCTIONS ////

    // Calculates the recoil energy for the given scattering direction
    inline CELER_FUNCTION real_type
    calc_recoil_energy(Real3 const& new_direction) const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from shared and state data
 */
CELER_FUNCTION
WokviInteractor::WokviInteractor(WokviRef const& shared,
                                 ParticleTrackView const& particle,
                                 Real3 const& inc_direction,
                                 MaterialView const& material,
                                 ElementComponentId const& elcomp_id,
                                 StackAllocator<Secondary>& allocate)
    : state_(particle,
             material,
             elcomp_id,
             /*cut_energy*/ units::MevEnergy{100.0},
             shared)
    , data_(shared)
    , inc_direction_(inc_direction)
    , allocate_(allocate)
    , inc_energy_(particle.energy())
    , inc_mass_(particle.mass())
    , target_(material...)
{
}

//---------------------------------------------------------------------------//
/*!
 * Sample the Coulomb scattering of the incident particle.
 */
template<class Engine>
CELER_FUNCTION Interaction WokviInteractor::operator()(Engine& rng)
{
    // Distribution model governing the scattering
    WokviDistribution distrib(data_);

    // Incident particle scatters
    Interaction result;

    // Sample the new direction
    const Real3 new_direction = distrib(rng);
    result.direction = rotate(inc_direction_, new_direction);

    // Calculate recoil and final energies
    real_type recoil_energy = calc_recoil_energy(new_direction);
    real_type final_energy = inc_energy_ - recoil_energy;
    if (final_energy < 0)
    {
        recoil_energy = inc_energy_;
        final_energy = 0;
    }
    result.energy = Energy{final_energy};

    // TODO: For high enough recoil energies, ions are produced

    result.energy_deposition = Energy{recoil_energy};

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculates the recoil energy for the given scattering direction calculated
 * by WokviDistribution.
 */
CELER_FUNCTION real_type
WokviInteractor::calc_recoil_energy(Real3 const& new_direction) const
{
    const real_type cos_theta = new_direction[2];
    const real_type inc_mom_sq = inc_energy_ * (inc_energy_ + 2 * inc_mass_);
    const real_type target_mass = value_as<Mass>(target_.nuclear_mass());
    return inc_mom_sq * (1-cos_theta)
           / (target_mass + (inc_mass_ + inc_energy_) * (1-cos_theta));
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
