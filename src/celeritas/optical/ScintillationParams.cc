//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/ScintillationParams.cc
//---------------------------------------------------------------------------//
#include "ScintillationParams.hh"

#include <algorithm>
#include <numeric>

#include "corecel/cont/Range.hh"
#include "corecel/data/CollectionBuilder.hh"
#include "corecel/math/SoftEqual.hh"
#include "celeritas/Types.hh"
#include "celeritas/grid/GenericGridBuilder.hh"
#include "celeritas/io/ImportData.hh"
#include "celeritas/phys/ParticleParams.hh"

#include "DataBuilder.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with imported data.
 */
std::shared_ptr<ScintillationParams>
ScintillationParams::from_import(ImportData const& data,
                                 SPConstParticles particle_params)
{
    CELER_EXPECT(!data.optical.empty());

    if (!std::any_of(
            data.optical.begin(), data.optical.end(), [](auto const& iter) {
                return static_cast<bool>(iter.second.scintillation);
            }))
    {
        // No scintillation data present
        return nullptr;
    }

    auto const& num_optmats = data.optical.size();

    Input input;
    input.resolution_scale.resize(num_optmats);
    if (data.optical_params.scintillation_by_particle)
    {
        // Collect ScintillationParticleIds
        input.pid_to_scintpid.resize(data.particles.size());
        ScintillationParticleId scintpid{0};
        for (auto const& [matid, iom] : data.optical)
        {
            auto const& iomsp = iom.scintillation.particles;
            for (auto const& [pdg, ipss] : iomsp)
            {
                if (auto const pid = particle_params->find(PDGNumber{pdg}))
                {
                    // Add new ScintillationParticleId
                    input.pid_to_scintpid[pid.get()] = scintpid++;
                }
            }
        }
        // Resize particle- and material-dependent spectra
        auto const num_scint_particles = scintpid.get();
        input.particles.resize(num_scint_particles * num_optmats);
    }
    else
    {
        // Resize material-only spectra
        input.materials.resize(num_optmats);
    }

    size_type optmatidx{0};
    for (auto const& [matid, iom] : data.optical)
    {
        input.resolution_scale[optmatidx] = iom.scintillation.resolution_scale;

        if (!data.optical_params.scintillation_by_particle)
        {
            // Material spectrum
            auto const& iomsm = iom.scintillation.material;
            ImportMaterialScintSpectrum mat_spec;
            mat_spec.yield_per_energy = iomsm.yield_per_energy;
            mat_spec.components = iomsm.components;
            input.materials[optmatidx] = std::move(mat_spec);
        }
        else
        {
            // Particle and material spectrum
            auto const& iomsp = iom.scintillation.particles;

            for (auto const& [pdg, ipss] : iomsp)
            {
                if (auto const pid = particle_params->find(PDGNumber{pdg}))
                {
                    auto scintpid = input.pid_to_scintpid[pid.get()];
                    CELER_ASSERT(scintpid);
                    ImportParticleScintSpectrum part_spec;
                    part_spec.yield_vector = ipss.yield_vector;
                    part_spec.components = ipss.components;
                    input.particles[num_optmats * scintpid.get() + optmatidx]
                        = std::move(part_spec);
                }
            }
        }
        optmatidx++;
    }

    return std::make_shared<ScintillationParams>(std::move(input));
}

//---------------------------------------------------------------------------//
/*!
 * Construct with scintillation input data.
 */
template<>
class DataBuilder<ScintillationData> : public GridDataBuilder<ScintillationData>
{
  public:
    using Self = DataBuilder<ScintillationData>;
    using ScintillationComponentId = OpaqueId<ScintillationComponent>;

  public:
    DataBuilder(HostVal<ScintillationData>* host_data)
        : GridDataBuilder<ScintillationData>(host_data)
    {
    }

    void operator()(ScintillationParams::Input const& input)
    {
        CELER_EXPECT(input);
        CELER_VALIDATE(input.particles.empty() != input.materials.empty(),
                       << "invalid input data. Material spectra and particle "
                          "spectra are mutually exclusive. Please store "
                          "either "
                          "material or particle spectra, but not both.");

        CELER_VALIDATE(
            std::all_of(input.resolution_scale.begin(),
                        input.resolution_scale.end(),
                        [](auto const& val) -> bool { return val >= 0; }),
            << "invalid resolution_scale for scintillation (should be "
               "nonnegative)");

        auto build_per_material = ListBuilder{
            range(OpticalMaterialId{input.resolution_scale.size()})};
        build_per_material(&host_data_->resolution_scale,
                           make_span(input.resolution_scale));
        CELER_ENSURE(!host_data_->resolution_scale.empty());

        if (input.particles.empty())
        {
            // Store material scintillation data
            build_per_material(&host_data_->materials,
                               make_span(input.materials),
                               *this,
                               &Self::import_material);

            CELER_VALIDATE(
                input.materials.size() == input.resolution_scale.size(),
                << "material and resolution scales do not match");
            CELER_ENSURE(host_data_->materials.size()
                         == host_data_->resolution_scale.size());
        }
        else
        {
            // Store particle data
            CELER_VALIDATE(!input.pid_to_scintpid.empty(),
                           << "missing particle ID to scintillation particle "
                              "ID "
                              "mapping");

            // Store particle ids
            host_data_->num_scint_particles = std::count_if(
                input.pid_to_scintpid.begin(),
                input.pid_to_scintpid.end(),
                [](auto const& id) { return static_cast<bool>(id); });
            CELER_ENSURE(host_data_->num_scint_particles > 0);

            ListBuilder{range(ParticleId{input.pid_to_scintpid.size()})}(
                &host_data_->pid_to_scintpid, make_span(input.pid_to_scintpid));

            // Store particle spectra
            CELER_EXPECT(!input.particles.empty());

            ListBuilder{range(ParticleScintSpectrumId{input.particles.size()})}(
                &host_data_->particles,
                make_span(input.particles),
                *this,
                &Self::import_particles);

            CELER_ENSURE(host_data_->particles.size()
                         == host_data_->num_scint_particles
                                * host_data_->resolution_scale.size());
        }
    }

  private:
    static ScintillationComponent
    import_component(ImportScintComponent const& input_comp,
                     real_type total_yield)
    {
        CELER_VALIDATE(input_comp.lambda_mean > 0,
                       << "invalid lambda_mean=" << input_comp.lambda_mean
                       << " for scintillation component (should be positive)");
        CELER_VALIDATE(input_comp.lambda_sigma > 0,
                       << "invalid lambda_sigma=" << input_comp.lambda_sigma
                       << " for scintillation component"
                       << " (should be positive)");
        CELER_VALIDATE(input_comp.rise_time >= 0,
                       << "invalid rise_time=" << input_comp.rise_time
                       << " for scintillation component"
                       << " (should be nonnegative)");
        CELER_VALIDATE(input_comp.fall_time > 0,
                       << "invalid fall_time=" << input_comp.fall_time
                       << " for scintillation component"
                       << " (should be positive)");
        CELER_VALIDATE(input_comp.yield_per_energy > 0,
                       << "invalid yield=" << input_comp.yield_per_energy
                       << " for scintillation component"
                       << " (should be positive)");

        ScintillationComponent comp;
        comp.lambda_mean = input_comp.lambda_mean;
        comp.lambda_sigma = input_comp.lambda_sigma;
        comp.rise_time = input_comp.rise_time;
        comp.fall_time = input_comp.fall_time;
        comp.yield_frac = input_comp.yield_per_energy / total_yield;
        return comp;
    }

    Range<ScintillationComponentId>
    import_components(std::vector<ImportScintComponent> const& comps) const
    {
        real_type total_yield = std::transform_reduce(
            comps.begin(), comps.end(), 0, std::plus{}, [](auto const& comp) {
                return comp.yield_per_energy;
            });

        return ListBuilder{range(ScintillationComponentId{comps.size()})}(
            &host_data_->components,
            make_span(comps),
            [total_yield](auto const& comp) {
                return Self::import_component(comp, total_yield);
            });
    }

    MaterialScintillationSpectrum
    import_material(ImportMaterialScintSpectrum const& mat) const
    {
        // Check validity of input scintillation data
        CELER_ASSERT(mat);
        CELER_VALIDATE(mat.yield_per_energy > 0,
                       << "invalid yield=" << mat.yield_per_energy
                       << " for scintillation (should be positive)");

        // Material-only data
        MaterialScintillationSpectrum mat_spec;
        mat_spec.yield_per_energy = mat.yield_per_energy;
        mat_spec.components = import_components(mat.components);
        return mat_spec;
    }

    ParticleScintillationSpectrum
    import_particles(ImportParticleScintSpectrum const& spec)
    {
        CELER_VALIDATE(spec.yield_vector,
                       << "particle yield vector is not assigned "
                          "correctly");
        ParticleScintillationSpectrum part_spec;
        part_spec.yield_vector = build_grid_(spec.yield_vector);
        part_spec.components = import_components(spec.components);
        return part_spec;
    };
};

ScintillationParams::ScintillationParams(Input const& input)
    : MirroredParamDataInterface<ScintillationData>(input)
{
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
