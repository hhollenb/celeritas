//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/distribution/MollerEnergyDistribution.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Macros.hh"
#include "corecel/Types.hh"
#include "corecel/math/Algorithms.hh"
#include "celeritas/em/interactor/detail/WokviStateHelper.hh"
#include "celeritas/em/xs/MottXsCalculator.hh"
#include "celeritas/em/xs/WokviXsCalculator.hh"
#include "celeritas/random/distribution/UniformRealDistribution.hh"

/*
 * References:
 *      Fern - Fernandez-Varea 1992
 *      PRM - Geant4 Physics Reference Manual Release 11.1
 */
namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Helper class for \c WokviInteractor .
 *
 * Calculates the total cross section of the Wentzel OK and VI model, and
 * samples the scattering direction.
 */
class WokviDistribution
{
  public:
    using MomentumSq = units::MevMomentumSq;

  public:
    // Construct with state and date from WokviInteractor
    inline CELER_FUNCTION
    WokviDistribution(detail::WokviStateHelper const& state,
                      WokviRef const& data);

    // Sample the scattering direction
    template<class Engine>
    inline CELER_FUNCTION Real3 operator()(Engine& rng) const;

    // The total cross section
    inline CELER_FUNCTION real_type cross_section() const;

  private:
    //// DATA ////

    // Precomputed variables for the interaction
    detail::WokviStateHelper const& state_;

    // Shared WokviModel data
    WokviRef const& data_;

    // Nuclear form factor
    const real_type form_factor_A_;

    // Ratio of electron to total cross section
    real_type elec_ratio_;

    // Sum of electron and nuclear cross sections
    real_type total_cross_section_;

    inline CELER_FUNCTION real_type calculate_form_factor(real_type formf,
                                                          real_type z1) const;
    inline CELER_FUNCTION real_type flat_form_factor(real_type x) const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct with state and data from WokviInteractor
 *
 * TODO: Reference for factors?
 */
CELER_FUNCTION
WokviDistribution::WokviDistribution(detail::WokviStateHelper const& state,
                                     WokviRef const& data)
    : state_(state)
    , data_(data)
    , form_factor_A_(state.form_factor_A())  // TODO:
                                             // Reference?
{
    // Calculate cross sections
    const WokviXsCalculator xsec(state);

    const real_type nuc_xsec = xsec.nuclear_xsec();
    const real_type elec_xsec = xsec.electron_xsec();

    total_cross_section_ = nuc_xsec + elec_xsec;
    if (total_cross_section_ > 0.0)
    {
        elec_ratio_ = elec_xsec / total_cross_section_;
    }
}

//---------------------------------------------------------------------------//
/*!
 * The total (nuclear + electron) cross section.
 */
CELER_FUNCTION real_type WokviDistribution::cross_section() const
{
    return total_cross_section_;
}

//---------------------------------------------------------------------------//
/*!
 * Samples the final direction of the interaction. This direction is in
 * the frame where the incident particle's momentum is oriented along the
 * z-axis, so it's final direction in the lab frame will need to be rotated.
 */
template<class Engine>
CELER_FUNCTION Real3 WokviDistribution::operator()(Engine& rng) const
{
    UniformRealDistribution<real_type> sample;

    // Parameters for scattering of a nucleus
    real_type form_factor = form_factor_A_;
    real_type cos_t1 = state_.cos_t_min_nuc();
    real_type cos_t2 = state_.cos_t_max_nuc();

    // Randomly choose if scattered off of electrons instead
    if (elec_ratio_ > 0.0 && sample(rng) <= elec_ratio_)
    {
        form_factor = 0.0;
        cos_t1 = state_.cos_t_min_elec();
        cos_t2 = state_.cos_t_max_elec();
    }

    // Check angular bounds are valid
    if (cos_t1 < cos_t2)
    {
        return {0.0, 0.0, 1.0};
    }

    // Sample scattering angle [Fern 92] where z1 = 2*mu = 1 - cos(t)
    const real_type w1 = state_.w_term(cos_t1);
    const real_type w2 = state_.w_term(cos_t2);
    const real_type z1 = w1 * w2 / (w1 + sample(rng) * (w2 - w1))
                         - state_.screening_coefficient();

    // Calculate rejection
    // TODO: Reference?
    MottXsCalculator mott_xsec(state_);
    const real_type fm = calculate_form_factor(form_factor, z1);
    const real_type g_rej = mott_xsec(sqrt(z1)) * ipow<2>(fm);

    if (sample(rng) > g_rej)
    {
        return {0.0, 0.0, 1.0};
    }

    // Calculate scattered vector assuming azimuthal angle is isotropic
    const real_type cos_t = clamp(1.0 - z1, -1.0, 1.0);
    const real_type sin_t = sqrt((1.0 - cos_t) * (1.0 + cos_t));
    const real_type phi = 2.0 * celeritas::constants::pi * sample(rng);
    return {sin_t * cos(phi), sin_t * sin(phi), cos_t};
}

//---------------------------------------------------------------------------//
/*!
 * Calculates the form factor based on the form factor model.
 * TODO: Reference?
 */
CELER_FUNCTION real_type
WokviDistribution::calculate_form_factor(real_type formf, real_type z1) const
{
    switch (data_.form_factor_type)
    {
        case NuclearFormFactorType::Flat: {
            // In units MeV
            const real_type ccoef = 0.00508;
            const real_type x = sqrt(2.0 * state_.inc_mom_sq * z1) * ccoef
                                * 2.0;
            return flat_form_factor(x)
                   * flat_form_factor(
                       x * 0.6 * fastpow(state_.target_mass(), 1.0 / 3.0));
        }
        break;
        case NuclearFormFactorType::Exponential:
            return 1.0 / ipow<2>(1.0 + formf * z1);
        case NuclearFormFactorType::Gaussian:
            return exp(-2.0 * formf * z1);
        default:
            return 1.0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Flat form factor
 * TODO: Reference?
 */
CELER_FUNCTION real_type WokviDistribution::flat_form_factor(real_type x) const
{
    return 3.0 * (sin(x) - x * cos(x)) / ipow<3>(x);
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
