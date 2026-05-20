//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file accel/gen/CherenkovOffload.hh
//---------------------------------------------------------------------------//
#pragma once

#include <optional>
#include <string>
#include <unordered_set>
#include <vector>
#include <G4Cerenkov.hh>
#include <G4LogicalVolume.hh>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * A replacement for Geant4's \c G4Cerenkov process which constructs \c
 * GeneratorDistributionData from a \c PostStepDoIt call.
 *
 * This process should have stacking photons set to false so that photons are
 * not initialized in Geant4.
 */
class CherenkovOffload : public G4Cerenkov
{
  public:
    using AllowedVolNames = std::vector<std::string>;
    using AllowedVols = std::unordered_set<G4LogicalVolume const*>;

    // Construct with optional volume whitelist
    CherenkovOffload(std::optional<AllowedVolNames> names = std::nullopt);

    // Prepare physics table for particle and enforce photon stacking
    void PreparePhysicsTable(G4ParticleDefinition const&) override;

    // Create a generator distribution for the given track and step
    G4VParticleChange*
    PostStepDoIt(G4Track const& aTrack, G4Step const& aStep) override;

  private:
    std::optional<AllowedVolNames> allowed_names_;
    std::optional<AllowedVols> allowed_vols_;
};

//---------------------------------------------------------------------------//
}  // namespace celeritas
