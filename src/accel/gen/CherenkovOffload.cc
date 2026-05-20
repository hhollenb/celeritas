//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file accel/gen/CherenkovOffload.cc
//---------------------------------------------------------------------------//
#include "CherenkovOffload.hh"

#include <G4LogicalVolumeStore.hh>

#include "corecel/io/Logger.hh"
#include "celeritas/g4/GeantOffloadUtils.hh"
#include "celeritas/optical/gen/GeneratorData.hh"
#include "accel/LocalOpticalGenOffload.hh"
#include "accel/detail/IntegrationSingleton.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with optional volume whitelist.
 */
CherenkovOffload::CherenkovOffload(std::optional<AllowedVolNames> names)
    : allowed_names_(std::move(names)), allowed_vols_(std::nullopt)
{
}

//---------------------------------------------------------------------------//
/*!
 * Prepare physics table for particle and enforce photon stacking.
 *
 * Defers physics table preparation to \c G4Cerenkov, but also enforces that
 * stacking photons is false afterwards.
 */
void CherenkovOffload::PreparePhysicsTable(G4ParticleDefinition const& particle)
{
    G4Cerenkov::PreparePhysicsTable(particle);

    // Enforce don't stack photons
    if (this->GetStackPhotons())
    {
        CELER_LOG(warning) << "CherenkovOffload requires stacking photons set "
                              "to false since it sends optical photon tracks "
                              "directly to Celeritas.";
        this->SetStackPhotons(false);
    }

    // Make list of allowed volumes
    if (allowed_names_)
    {
        auto const* volume_store = G4LogicalVolumeStore::GetInstance();
        CELER_ASSERT(volume_store);
        CELER_ASSERT(volume_store->IsMapValid());

        allowed_vols_.emplace();
        for (auto const& name : *allowed_names_)
        {
            auto const* vol = volume_store->GetVolume(name, false);
            CELER_VALIDATE(vol,
                           << "could not find Geant4 logical volume with name "
                              "\""
                           << name << "\" for CherenkovOffload whitelist");
            allowed_vols_->insert(vol);
        }

        CELER_ENSURE(allowed_vols_->size() == allowed_names_->size());
    }
}

//---------------------------------------------------------------------------//
/*!
 * Create a generator distribution for the given track and step.
 *
 * Stacking photons should be disabled so that photons are not duplicated in
 * Geant4. After calling the \c G4Cerenkov::PostStepDoIt this function creates
 * a \c GeneratorDistributionData and pushes it to the local offload, which
 * should be \c LocalOpticalGenOffload.
 */
G4VParticleChange*
CherenkovOffload::PostStepDoIt(G4Track const& aTrack, G4Step const& aStep)
{
    CELER_EXPECT(!this->GetStackPhotons());

    auto* result = G4Cerenkov::PostStepDoIt(aTrack, aStep);

    // If whitelist exists, check if step is inside it
    if (allowed_vols_)
    {
        auto const* vol
            = aStep.GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
        if (allowed_vols_->count(vol) == 0)
        {
            return result;
        }
    }

    if (this->GetNumPhotons() > 0)
    {
        auto data = distribution_from_step(aStep);
        data.type = GeneratorType::cherenkov;
        data.num_photons = static_cast<size_type>(this->GetNumPhotons());

        // Push generator distribution for this step to offload
        auto& local = detail::IntegrationSingleton::instance().local_offload();
        auto* gen_offload = dynamic_cast<LocalOpticalGenOffload*>(&local);

        CELER_VALIDATE(gen_offload,
                       << "LocalOpticalGenOffload required for "
                          "CherenkovOffload");

        CELER_LOG_LOCAL(debug)
            << "Offloading " << data.num_photons << " Cherenkov photons";

        gen_offload->Push(data);
    }

    return result;
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
