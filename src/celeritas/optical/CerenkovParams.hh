//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/CerenkovParams.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Types.hh"
#include "corecel/data/CollectionMirror.hh"
#include "corecel/data/ParamsDataInterface.hh"

#include "CerenkovData.hh"
#include "DataBuilder.hh"

namespace celeritas
{
class OpticalPropertyParams;

//---------------------------------------------------------------------------//
/*!
 * Build and manage Cerenkov data.
 */
class CerenkovParams final : public MirroredParamDataInterface<CerenkovData>
{
  public:
    //!@{
    //! \name Type aliases
    using SPConstProperties = std::shared_ptr<OpticalPropertyParams const>;
    //!@}

  public:
    // Construct with optical property data
    explicit CerenkovParams(SPConstProperties properties);
};

//---------------------------------------------------------------------------//
}  // namespace celeritas
