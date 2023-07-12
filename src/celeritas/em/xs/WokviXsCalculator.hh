//----------------------------------*-C++-*----------------------------------//
// Copyright 2023-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/xs/WokviXsCalculator.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Macros.hh"
#include "corecel/Types.hh"
#include "corecel/math/Algorithms.hh"
#include "celeritas/phys/AtomicNumber.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Calculate the Coulomb scattering cross section using the Wentzel OK and
 * VI model.
 */
class WokviXsCalculator
{
  public:
    // Construct the calculator from the given values
    inline CELER_FUNCTION WokviXsCalculator(int target_z,
                                            real_type screening_coefficient,
                                            real_type cos_t_max_elec);

    // The ratio of electron to total cross section for Coulomb scattering
    inline CELER_FUNCTION real_type operator()() const;

  private:
    int const target_z_;
    real_type const screening_coefficient_;
    real_type const cos_t_max_elec_;

    // Cross sections
    inline CELER_FUNCTION real_type nuclear_xsec() const;
    inline CELER_FUNCTION real_type electron_xsec() const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct with state data
 */
CELER_FUNCTION
WokviXsCalculator::WokviXsCalculator(int target_z,
                                     real_type screening_coefficient,
                                     real_type cos_t_max_elec)
    : target_z_(target_z)
    , screening_coefficient_(screening_coefficient)
    , cos_t_max_elec_(cos_t_max_elec)
{
    CELER_EXPECT(target_z_ > 0);
    CELER_EXPECT(screening_coefficient > 0);
    CELER_EXPECT(cos_t_max_elec >= -1 && cos_t_max_elec <= 1);
}


//---------------------------------------------------------------------------//
/*!
 * Ratio of electron cross section to the total (nuclear + electron)
 * cross section.
 */
CELER_FUNCTION real_type WokviXsCalculator::operator()() const
{
    const real_type nuc_xsec = nuclear_xsec();
    const real_type elec_xsec = electron_xsec();

    return elec_xsec / (nuc_xsec + elec_xsec);
}

//---------------------------------------------------------------------------//
/*!
 * Reduced integrated nuclear cross section from theta=0 to pi.
 * Since this is only used in the electric ratio, mutual factors with the
 * electron cross section are dropped.
 */
CELER_FUNCTION real_type WokviXsCalculator::nuclear_xsec() const
{
    return target_z_ / (1 + screening_coefficient_);
}

//---------------------------------------------------------------------------//
/*!
 * Integrated electron cross section
 */
CELER_FUNCTION real_type WokviXsCalculator::electron_xsec() const
{
    return (1 - cos_t_max_elec_)
           / (1 - cos_t_max_elec_ + 2 * screening_coefficient_);
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
