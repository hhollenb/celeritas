//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/xs/WokviXsCalculator.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Macros.hh"
#include "corecel/Types.hh"
#include "celeritas/em/interactor/detail/WokviStateHelper.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Calculates the ratio of Mott cross section to the Rutherford cross section.
 */
class MottXsCalculator
{
  public:
    // Construct with state data
    inline CELER_FUNCTION
    MottXsCalculator(detail::WokviStateHelper const& state);

    // Ratio of Mott and Rutherford cross sections
    inline CELER_FUNCTION real_type operator()(real_type fcos_t) const;

  private:
    detail::WokviStateHelper const& state_;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct with state data
 */
CELER_FUNCTION
MottXsCalculator::MottXsCalculator(detail::WokviStateHelper const& state)
    : state_(state)
{
}

//---------------------------------------------------------------------------//
/*!
 * Compute the ratio of Mott to Rutherford cross sections.
 *
 * The parameter fcos_t is equivalent to
 *      sqrt(1 - cos(theta))
 * where theta is the scattered angle in the z-aligned momentum frame.
 *
 * For 1 <= Z <= 92, an interpolated expression is used [PRM 8.48].
 */
CELER_FUNCTION
real_type MottXsCalculator::operator()(real_type fcos_t) const
{
    real_type ratio = 0;
    // Mean velocity of electrons between ~KeV and 900 MeV
    const real_type beta_shift = 0.7181228;
    const real_type beta0 = sqrt(1.0 / state_.inv_beta_sq) - beta_shift;

    // Construct [0,5] powers of beta0
    real_type b[6];
    b[0] = 1.0;
    for (int i = 1; i < 6; i++)
    {
        b[i] = beta0 * b[i - 1];
    }

    // Compute the ratio, summing over powers of fcos_t
    real_type f0 = 1.0;
    for (int j = 0; j <= 4; j++) {
        // Calculate the a_j coefficient
        real_type a = 0.0;
        for (int k = 0; k < 6; k++) {
            a += state_.element_data.mott_coeff[j][k] * b[k];
        }
        // Sum in power series of fcos_t
        ratio += a * f0;
        f0 *= fcos_t;
    }

    return ratio;
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
