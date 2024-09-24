//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/model/RayleighModel.hh
//---------------------------------------------------------------------------//
#pragma once

#include "../ImportedModelAdapter.hh"
#include "../Model.hh"

namespace celeritas
{
namespace optical
{
//---------------------------------------------------------------------------//
/*!
 * Set up and launch the optical Rayleigh model interaction.
 */
class RayleighModel : public Model
{
  public:
    // Construct with imported data
    RayleighModel(ActionId id, ImportedModelAdapter imported);

    // Build the mean free paths for this model
    void build_mfps(detail::MfpBuilder) const override final;

    // Execute the model with host data
    void step(CoreParams const&, CoreStateHost&) const override final;

    // Execute the model with device data
    void step(CoreParams const&, CoreStateDevice&) const override final;

  private:
    ImportedModelAdapter imported_;
};

//---------------------------------------------------------------------------//
}  // namespace optical
}  // namespace celeritas
