//----------------------------------*-C++-*----------------------------------//
// Copyright 2023-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/model/WokviModel.cc
//---------------------------------------------------------------------------//
#include "WokviModel.hh"

#include "celeritas_config.h"
#include "corecel/sys/ScopedMem.hh"
#include "celeritas/em/data/WokviData.hh"
#include "celeritas/em/executor/WokviExecutor.hh"
#include "celeritas/em/model/detail/MottInterpolatedCoefficients.hh"
#include "celeritas/global/ActionLauncher.hh"
#include "celeritas/global/CoreParams.hh"
#include "celeritas/global/TrackExecutor.hh"
#include "celeritas/io/ImportProcess.hh"
#include "celeritas/mat/MaterialParams.hh"
#include "celeritas/phys/InteractionApplier.hh"
#include "celeritas/phys/PDGNumber.hh"
#include "celeritas/phys/ParticleParams.hh"
#include "celeritas/phys/ParticleView.hh"

namespace celeritas
{

//---------------------------------------------------------------------------//
WokviModel::WokviModel(ActionId id,
                       ParticleParams const& particles,
                       MaterialParams const& materials,
                       SPConstImported data)
    : imported_(data,
                particles,
                ImportProcessClass::coulomb_scat,  // TODO: Check the
                                                   // ImportProcessClass tags
                ImportModelClass::e_coulomb_scattering,  // TODO: Check the
                                                         // ImportModelClass
                                                         // tags
                {pdg::electron(), pdg::positron()})
{
    CELER_EXPECT(id);

    ScopedMem record_mem("WokviModel.construct");
    HostVal<WokviData> host_data;

    // This is where the data is built and transfered to the device
    host_data.ids.action = id;
    host_data.ids.electron = particles.find(pdg::electron());
    host_data.ids.positron = particles.find(pdg::positron());

    CELER_VALIDATE(host_data.ids,
                   << "missing IDs (required for " << this->description()
                   << ")");

    // Electron mass
    host_data.electron_mass = particles.get(host_data.ids.electron).mass();

    // TODO: Select form factor
    host_data.form_factor_type = NuclearFormFactorType::Exponential;

    // Prefactor of the screen R squared
    // const real_type fact = 1.0;  // TODO:
    // G4EmParameters::Instance()->ScreeningFactor();

    // Thomas-Fermi constant C_TF
    const real_type ctf
        = fastpow(3 * constants::pi / 4, static_cast<real_type>(2) / 3) / 2;
    std::cout
        << "MODEL: \n"
        << "\tCTF = " << ctf << "\n"
        << "\tElectron Mass = "
        << native_value_to<units::MevMass>(constants::electron_mass).value()
        << "\n"
        << "\talpha2 = " << ipow<2>(constants::alpha_fine_structure) << "\n"
        << "\tCeleritas R^2: "
        << native_value_to<units::MevMomentumSq>(
               ipow<2>(constants::hbar_planck / (2 * ctf * constants::a0_bohr)))
               .value()
        << "\n"
        << "\tMatching R^2: "
        << 0.25
               * ipow<2>(
                   constants::alpha_fine_structure
                   * native_value_to<units::MevMass>(constants::electron_mass)
                         .value()
                   / ctf)
        << "\n";
    host_data.screen_r_sq_elec = native_value_to<units::MevMomentumSq>(
        ipow<2>(constants::hbar_planck / (2 * ctf * constants::a0_bohr)));

    // This is the inverse of Geant's constn
    // need to multiply by 2 to match Geant's magic number
    std::cout
        << "\n\n"
        << "\tForm Momentum Scale: "
        << native_value_to<units::MevMomentumSq>(
               12
               / ipow<2>(1.27 * (1e-15 * units::meter) / constants::hbar_planck)
               / 2)
               .value()
        << "\n";
    host_data.form_momentum_scale = native_value_to<units::MevMomentumSq>(
        12.0 / ipow<2>(1.27 * (1e-15 * units::meter) / constants::hbar_planck)
        / 2.0);

    build_data(host_data, materials);

    data_ = CollectionMirror<WokviData>{std::move(host_data)};

    CELER_ENSURE(this->data_);
}

auto WokviModel::applicability() const -> SetApplicability
{
    Applicability electron_applic;
    electron_applic.particle = this->host_ref().ids.electron;
    // TODO: construct actual energy range
    electron_applic.lower = zero_quantity();
    electron_applic.upper = max_quantity();

    Applicability positron_applic;
    positron_applic.particle = this->host_ref().ids.positron;
    positron_applic.lower = zero_quantity();
    positron_applic.upper = max_quantity();

    return {electron_applic, positron_applic};
}

auto WokviModel::micro_xs(Applicability applic) const -> MicroXsBuilders
{
    return imported_.micro_xs(std::move(applic));
}

void WokviModel::execute(CoreParams const& params, CoreStateHost& state) const
{
    auto execute = make_action_track_executor(
        params.ptr<MemSpace::native>(),
        state.ptr(),
        this->action_id(),
        InteractionApplier{WokviExecutor{this->host_ref()}});
    return launch_action(*this, params, state, execute);
}

#if !CELER_USE_DEVICE
void WokviModel::execute(CoreParams const&, CoreStateDevice&) const
{
    CELER_NOT_CONFIGURED("CUDA OR HIP");
}
#endif

ActionId WokviModel::action_id() const
{
    return this->host_ref().ids.action;
}

void WokviModel::build_data(HostVal<WokviData>& host_data,
                            MaterialParams const& materials)
{
    // Build element data (?)
    unsigned int const num_elements = materials.num_elements();
    auto elem_data = make_builder(&host_data.elem_data);
    elem_data.reserve(num_elements);

    // Thomas-Fermi screening radii
    // Formfactors from A.V. Butkevich et al., NIM A 488 (2002) 282
    // TODO: Reference? Check math, put into helper function
    // Numerical factor in the form factor (units (MeV/c)^-2)
    // const real_type constn = 6.937e-6; // magic number?

    for (auto el_id : range(ElementId{num_elements}))
    {
        ElementView const& element = materials.get(el_id);
        const AtomicNumber z = element.atomic_number();

        WokviElementData z_data;

        // Load Mott coefficients
        // Currently only support up to Z=92 (Uranium) as taken from Geant4
        int const index = (z.get() <= 92) ? z.get() : 0;
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                z_data.mott_coeff[i][j]
                    = detail::interpolated_mott_coeffs[index][i][j];
            }
        }

        elem_data.push_back(z_data);
    }
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
