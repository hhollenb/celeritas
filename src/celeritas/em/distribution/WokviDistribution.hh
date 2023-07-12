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
#include "celeritas/em/data/WokviData.hh"
#include "celeritas/em/xs/MottXsCalculator.hh"
#include "celeritas/em/xs/WokviXsCalculator.hh"
#include "celeritas/mat/IsotopeView.hh"
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
    using Energy = units::MevEnergy;
    using Mass = units::MevMass;

  public:
    // Construct with state and date from WokviInteractor
    inline CELER_FUNCTION
    WokviDistribution(real_type inc_energy,
                      real_type inc_mass,
                      IsotopeView const& target,
                      WokviElementData const& element_data,
                      WokviRef const& data);

    // Sample the scattering direction
    template<class Engine>
    inline CELER_FUNCTION Real3 operator()(Engine& rng) const;

  private:
    //// DATA ////

    // Shared WokviModel data
    WokviRef const& data_;

    // Target element
    IsotopeView const& target_;
    WokviElementData const& element_data_;

    // Incident particle data
    real_type const inc_energy_;
    real_type const inc_mass_;


    inline CELER_FUNCTION real_type calculate_form_factor(real_type formf,
                                                          real_type cos_t) const;
    inline CELER_FUNCTION real_type flat_form_factor(real_type x) const;

    inline CELER_FUNCTION real_type compute_screening_coefficient() const;
    inline CELER_FUNCTION real_type compute_max_electron_cos_t() const;

    inline CELER_FUNCTION int target_Z() const;
    inline CELER_FUNCTION real_type inc_mom_sq() const;
};

namespace
{
CELER_FUNCTION real_type calc_mom_sq(real_type energy, real_type mass)
{
    return energy * (energy + 2 * mass);
}
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
WokviDistribution::WokviDistribution(real_type inc_energy,
                                     real_type inc_mass,
                                     IsotopeView const& target,
                                     WokviElementData const& element_data,
                                     WokviRef const& data)
    : data_(data)
    , target_(target)
    , element_data_(element_data)
    , inc_energy_(inc_energy)
    , inc_mass_(inc_mass)
{}

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

    const real_type screen_coeff = compute_screening_coefficient();
    const real_type cos_t_max_elec = compute_max_electron_cos_t();

    const real_type scale
        = (target_Z() == 1) ? 1 / 3.097e-6
                            : value_as<MomentumSq>(data_.form_momentum_scale);
    const real_type form_factor_A
        = inc_mom_sq()
          * fastpow((real_type)target_.atomic_mass_number().get(), 2 * 0.27)
          / scale;

    // Parameters for scattering of a nucleus
    real_type form_factor = form_factor_A;
    real_type cos_t1 =  1;
    real_type cos_t2 = -1;

    // Randomly choose if scattered off of electrons instead
    const WokviXsCalculator xsec(target_Z(), screen_coeff, cos_t_max_elec);
    const real_type elec_ratio = xsec();
    if (sample(rng) < elec_ratio)
    {
        form_factor = 0;
        // TODO: Can simplify this logic with
        //       -1 <= cos_t_max_elec <= 1
        cos_t1 = max(cos_t1, cos_t_max_elec);
        cos_t2 = max(cos_t2, cos_t_max_elec);
    }

    // Sample scattering angle [Fern 92] where cos(theta) = 1 + 2*mu
    // For incident electrons / positrons, theta_min = 0 always
    const real_type w1 = 1 - cos_t1 + 2 * screen_coeff;
    const real_type w2 = 1 - cos_t2 + 2 * screen_coeff;
    const real_type cos_theta = clamp(
        1 + 2 * screen_coeff - w1 * w2 / (w1 + sample(rng) * (w2 - w1)),
        {-1},
        {1});

    // Calculate rejection
    // TODO: Reference?
    MottXsCalculator mott_xsec(element_data_, inc_energy_, inc_mass_);
    const real_type fm = calculate_form_factor(form_factor, cos_theta);
    const real_type g_rej = mott_xsec(cos_theta) * ipow<2>(fm);

    if (sample(rng) > g_rej)
    {
        return {0, 0, 1};
    }

    // Calculate scattered vector assuming azimuthal angle is isotropic
    const real_type sin_theta = sqrt((1 - cos_theta) * (1 + cos_theta));
    const real_type phi = 2 * celeritas::constants::pi * sample(rng);
    return {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
}

//---------------------------------------------------------------------------//
/*!
 * Calculates the form factor based on the form factor model.
 * TODO: Reference?
 */
CELER_FUNCTION real_type
WokviDistribution::calculate_form_factor(real_type formf, real_type cos_t) const
{
    switch (data_.form_factor_type)
    {
        case NuclearFormFactorType::Flat: {
            // In units MeV
            const real_type ccoef = 0.00508;
            const real_type x = sqrt(2 * inc_mom_sq() * (1-cos_t)) * ccoef * 2;
            return flat_form_factor(x)
                   * flat_form_factor(
                       x * 0.6
                       * fastpow(value_as<Mass>(target_.nuclear_mass()),
                                 static_cast<real_type>(1) / 3));
        }
        break;
        case NuclearFormFactorType::Exponential:
            return 1 / ipow<2>(1 + formf * (1-cos_t));
        case NuclearFormFactorType::Gaussian:
            return exp(-2 * formf * (1-cos_t));
        default:
            return 1;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Flat form factor
 * TODO: Reference?
 */
CELER_FUNCTION real_type WokviDistribution::flat_form_factor(real_type x) const
{
    return 3 * (sin(x) - x * cos(x)) / ipow<3>(x);
}

//---------------------------------------------------------------------------//
CELER_FUNCTION int WokviDistribution::target_Z() const
{
    return target_.atomic_number().get();
}

CELER_FUNCTION real_type WokviDistribution::inc_mom_sq() const
{
    return calc_mom_sq(inc_energy_, inc_mass_);
}

//---------------------------------------------------------------------------//
CELER_FUNCTION real_type WokviDistribution::compute_screening_coefficient() const
{
    // TODO: Reference for just proton correction?
    real_type correction = 1;
    const real_type sq_cbrt_z = fastpow(static_cast<real_type>(target_Z()),
                                        static_cast<real_type>(2) / 3);
    if (target_Z() > 1)
    {
        const real_type tau = inc_energy_ / inc_mass_;
        const real_type factor = sqrt(tau / (tau + sq_cbrt_z));
        const real_type inv_beta_sq = 1 + ipow<2>(inc_mass_) / inc_mom_sq();
        correction = min(
                         target_Z() * 1.13,
                         1.13 + 3.76 * ipow<2>(target_Z()) * inv_beta_sq * constants::alpha_fine_structure * factor);
    }

    return correction * value_as<MomentumSq>(data_.screen_r_sq_elec)
           * sq_cbrt_z / inc_mom_sq();
}

//---------------------------------------------------------------------------//
CELER_FUNCTION real_type WokviDistribution::compute_max_electron_cos_t() const
{
    // TODO: Cut Energy?
    const Energy cut_energy;
    const real_type t_max = 0.5 * inc_energy_;
    const real_type t = min(value_as<Energy>(cut_energy), t_max);
    const real_type t1 = inc_energy_ - t;
    if (t1 > 0)
    {
        const real_type mom1_sq = calc_mom_sq(t, value_as<Mass>(data_.electron_mass));
        const real_type mom2_sq = calc_mom_sq(t1, value_as<Mass>(data_.electron_mass));
        const real_type ctm = (inc_mom_sq() + mom2_sq - mom1_sq) * 0.5
                              / sqrt(inc_mom_sq() * mom2_sq);
        return clamp(ctm, real_type{0}, real_type{1});
    }

    // Default value
    return 1;
}


//---------------------------------------------------------------------------//
}  // namespace celeritas
