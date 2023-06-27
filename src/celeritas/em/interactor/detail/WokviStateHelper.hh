//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/interactor/detail/BremFinalStateHelper.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/math/Algorithms.hh"
#include "celeritas/Constants.hh"
#include "celeritas/Quantities.hh"
#include "celeritas/Types.hh"
#include "celeritas/em/data/WokviData.hh"
#include "celeritas/mat/ElementView.hh"
#include "celeritas/mat/MaterialView.hh"
#include "celeritas/phys/ParticleTrackView.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * A utility class for containing common quantities needed by different parts
 * of the Wokvi Model.
 */
struct WokviStateHelper
{
  public:
    //!@{
    //! \name Type aliases
    using Energy = units::MevEnergy;
    using Mass = units::MevMass;
    using MomentumSq = units::MevMomentumSq;
    using Charge = units::ElementaryCharge;
    //!@}

  public:
    // Construct from data
    inline CELER_FUNCTION WokviStateHelper(ParticleTrackView const& particle,
                                           MaterialView const& material,
                                           ElementComponentId const& elcomp_id,
                                           Energy cut_energy,
                                           WokviRef const& data);

    // Incident particle quantities
    const real_type inc_energy;
    const real_type inc_mass;
    const real_type inc_mom_sq;
    const real_type inv_beta_sq;

    // Target element quantities
    const ElementView element;
    WokviElementData const& element_data;
    inline CELER_FUNCTION real_type target_Z() const
    {
        return element.atomic_number().get();
    }
    // TODO: AmuMass or MevMass?
    inline CELER_FUNCTION real_type target_mass() const
    {
        return value_as<units::AmuMass>(element.atomic_mass());
    }

    // Kinetic factor in cross sections
    const real_type kinetic_factor;

    // (Twice the) Screening coefficient (2*A in [PRM 8.51])
    inline CELER_FUNCTION real_type screening_coefficient() const
    {
        return screen_z_;
    }

    // Upper bound of cos(theta) for scattering off electrons
    inline CELER_FUNCTION real_type max_electron_cos_t() const
    {
        return cos_t_max_elec_;
    }

    // Mott factor for electron/positron scattering
    inline CELER_FUNCTION real_type mott_factor() const
    {
        return mott_factor_;
    }

    // cos(theta) bounds for nuclear cross section
    inline CELER_FUNCTION real_type cos_t_min_nuc() const;
    inline CELER_FUNCTION real_type cos_t_max_nuc() const;

    // cos(theta) bounds for electron cross section
    inline CELER_FUNCTION real_type cos_t_min_elec() const;
    inline CELER_FUNCTION real_type cos_t_max_elec() const;

    // Helper function for a common used expression
    //      w = 1 + 2A - cos_t
    inline CELER_FUNCTION real_type w_term(real_type cos_t) const;

  private:
    // Constat model data
    WokviRef const& data_;

    // Calculated quantities
    real_type screen_z_;
    real_type cos_t_max_elec_;
    real_type mott_factor_;

    // Computes the value of the screening coefficient
    inline CELER_FUNCTION real_type compute_screening_coefficient() const;

    // Scaling behavior of the screening coefficient
    // Helper function for compute_screening_coefficient
    inline CELER_FUNCTION real_type
    screening_scaling_function(real_type factor) const;

    // Computes the maximum cos(theta) for scattering off electrons
    inline CELER_FUNCTION real_type
    compute_max_electron_cos_t(Energy cut_energy) const;
};

//---------------------------------------------------------------------------//
// UTILITY FUNCTIONS
//---------------------------------------------------------------------------//
namespace
{
/*!
 * Calculate the momentum squared from the kinetic energy and mass.
 */
CELER_FUNCTION real_type momentum_from_kinetic_energy(real_type t,
                                                      real_type mass)
{
    return t * (t + 2.0 * mass);
}
}  // namespace

//---------------------------------------------------------------------------//
/*!
 * Construct the state from track and shared data.
 */
CELER_FUNCTION
WokviStateHelper::WokviStateHelper(ParticleTrackView const& particle,
                                   MaterialView const& material,
                                   ElementComponentId const& elcomp_id,
                                   Energy cut_energy,
                                   WokviRef const& data)
    : inc_energy(value_as<Energy>(particle.energy()))
    , inc_mass(value_as<Mass>(particle.mass()))
    , inc_mom_sq(value_as<MomentumSq>(particle.momentum_sq()))
    , inv_beta_sq(1.0 + ipow<2>(inc_mass) / inc_mom_sq)
    , element(material.make_element_view(elcomp_id))
    , element_data(data.elem_data[material.element_id(elcomp_id)])
    , kinetic_factor(value_as<CoeffQuantity>(data.coeff)
                     * element.atomic_number().get()
                     * ipow<2>(value_as<Charge>(particle.charge()))
                     * inv_beta_sq / inc_mom_sq)
    , data_(data)
{
    // Order of initialization doesn't matter
    screen_z_ = compute_screening_coefficient();
    cos_t_max_elec_ = compute_max_electron_cos_t(cut_energy);

    // Mott factor for incident electron / positrons
    // TODO: Reference?
    mott_factor_ = 1.0 + 2.0e-4 * ipow<2>(target_Z());
}

//---------------------------------------------------------------------------//
/*!
 * Cosine of the minimum theta for nuclear scattering
 */
CELER_FUNCTION real_type WokviStateHelper::cos_t_min_nuc() const
{
    return 1.0;
}

//---------------------------------------------------------------------------//
/*!
 * Cosine of the maximum theta for nuclear scattering
 */
CELER_FUNCTION real_type WokviStateHelper::cos_t_max_nuc() const
{
    return -1.0;
}

//---------------------------------------------------------------------------//
/*!
 * Cosine of the minimum theta for electron scattering
 */
CELER_FUNCTION real_type WokviStateHelper::cos_t_min_elec() const
{
    return max(cos_t_min_nuc(), cos_t_max_elec_);
}

//---------------------------------------------------------------------------//
/*!
 * Cosine of the maximum theta for electron scattering
 */
CELER_FUNCTION real_type WokviStateHelper::cos_t_max_elec() const
{
    return max(cos_t_max_nuc(), cos_t_max_elec_);
}

//---------------------------------------------------------------------------//
/*!
 * Computes and saves the cosine of the maximum theta for electron scattering.
 *
 * TODO: Reference?
 */
CELER_FUNCTION real_type
WokviStateHelper::compute_max_electron_cos_t(Energy cut_energy) const
{
    const real_type t_max = 0.5 * inc_energy;
    const real_type t = min(value_as<Energy>(cut_energy), t_max);
    const real_type t1 = inc_energy - t;
    if (t1 > 0.0)
    {
        const real_type mom1_sq = momentum_from_kinetic_energy(
            t, value_as<Mass>(data_.electron_mass));
        const real_type mom2_sq = momentum_from_kinetic_energy(
            t1, value_as<Mass>(data_.electron_mass));
        const real_type ctm = (inc_mom_sq + mom2_sq - mom1_sq) * 0.5
                              / sqrt(inc_mom_sq * mom2_sq);
        return clamp(ctm, 0.0, 1.0);
    }

    // Default value
    return 1.0;
}

//---------------------------------------------------------------------------//
/*!
 * Helper function that calculates part of the Moliere and Bethe screening
 * factor, specifically the factor in the square brackets in
 * [Fern 32] or in [PRM 8.22]
 */
CELER_FUNCTION real_type
WokviStateHelper::screening_scaling_function(real_type factor) const
{
    // TODO: Why min?
    return min(target_Z() * 1.13,
               1.13
                   + 3.76 * ipow<2>(target_Z()) * inv_beta_sq
                         * constants::alpha_fine_structure * factor);
}

//---------------------------------------------------------------------------//
/*!
 * Computes the screening coefficient 2A as described in [PRM 8.22]
 */
CELER_FUNCTION real_type WokviStateHelper::compute_screening_coefficient() const
{
    // TODO: Reference for just proton correction?
    real_type correction = 1.0;
    if (element.atomic_number().get() > 1)
    {
        const real_type tau = inc_energy / inc_mass;
        correction = screening_scaling_function(
            sqrt(tau / (tau + ipow<2>(element.cbrt_z()))));
    }

    return correction * value_as<MomentumSq>(element_data.screen_r_sq_elec)
           / inc_mom_sq;
}

//---------------------------------------------------------------------------//
/*!
 * A convenience function for computing the term
 *      w = 1.0 + 2A - cos(t)
 * which occurs regularly in the cross section calculations.
 */
CELER_FUNCTION real_type WokviStateHelper::w_term(real_type cos_t) const
{
    return 1.0 + screen_z_ - cos_t;
}

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace celeritas
