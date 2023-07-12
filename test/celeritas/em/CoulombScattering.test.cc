//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/em/CoulombScattering.test.cc
//---------------------------------------------------------------------------//
#include "celeritas/Quantities.hh"
#include "celeritas/Units.hh"
#include "celeritas/em/interactor/WokviInteractor.hh"
#include "celeritas/em/model/WokviModel.hh"
#include "celeritas/em/model/detail/MottInterpolatedCoefficients.hh"
#include "celeritas/mat/MaterialTrackView.hh"
#include "celeritas/mat/MaterialView.hh"
#include "celeritas/phys/InteractionIO.hh"
#include "celeritas/phys/InteractorHostTestBase.hh"

#include "celeritas_test.hh"

namespace celeritas
{
namespace test
{
//---------------------------------------------------------------------------//

class CoulombScatteringTest : public InteractorHostTestBase
{
  protected:
    void SetUp() override
    {
        // Set up shared material data
        // TODO: Use multiple elements to test elements are picked correctly
        MaterialParams::Input mat_inp;
        mat_inp.elements
            = {{AtomicNumber{29}, units::AmuMass{63.546}, {}, "Cu"}};
        mat_inp.materials = {
            {0.141 * constants::na_avogadro,
             293.0,
             MatterState::solid,
             {{ElementId{0}, 1.0}},
             "Cu"},
        };
        this->set_material_params(mat_inp);

        // Create mock import data
        {
            ImportProcess ip_electron = this->make_import_process(
                pdg::electron(),
                pdg::electron(),  // TODO: no secondary?
                ImportProcessClass::coulomb_scat,
                {ImportModelClass::e_coulomb_scattering});
            ImportProcess ip_positron = ip_electron;
            ip_positron.particle_pdg = pdg::positron().get();
            this->set_imported_processes(
                {std::move(ip_electron), std::move(ip_positron)});
        }

        model_ = std::make_shared<WokviModel>(ActionId{0},
                                              *this->particle_params(),
                                              *this->material_params(),
                                              this->imported_processes());

        // Set cutoffs
        CutoffParams::Input input;
        CutoffParams::MaterialCutoffs material_cutoffs;
        material_cutoffs.push_back({MevEnergy{0.02064384}, 0.07});
        input.materials = this->material_params();
        input.particles = this->particle_params();
        input.cutoffs.insert({pdg::gamma(), material_cutoffs});
        this->set_cutoff_params(input);

        // Set incident particle to be an electron at 10 MeV
        this->set_inc_particle(pdg::electron(), MevEnergy{10.0});
        this->set_inc_direction({0, 0, 1});
        this->set_material("Cu");
    }

  protected:
    std::shared_ptr<WokviModel> model_;
};

TEST_F(CoulombScatteringTest, wokvi_data)
{
    WokviHostRef const& data = model_->host_ref();

    EXPECT_SOFT_EQ(constants::electron_mass,
                   native_value_from(data.electron_mass));
    EXPECT_EQ(NuclearFormFactorType::Exponential, data.form_factor_type);
    EXPECT_SOFT_EQ(1.0 / 6.937e-6, data.form_momentum_scale.value());
    EXPECT_SOFT_EQ(8.87463e-6, data.screen_r_sq_elec.value());

    // Check element data is filled in correctly
    unsigned int const num_elements = this->material_params()->num_elements();
    for (auto el_id : range(ElementId(num_elements)))
    {
        int const z = this->material_params()->get(el_id).atomic_number().get();
        int const mott_index = (z <= 92) ? z : 0;

        WokviElementData const& element_data = data.elem_data[el_id];
        for (auto i : range(5))
        {
            for (auto j : range(6))
            {
                EXPECT_EQ(detail::interpolated_mott_coeffs[mott_index][i][j],
                          element_data.mott_coeff[i][j]);
            }
        }
    }
}

TEST_F(CoulombScatteringTest, wokvi_xs)
{
    WokviXsCalculator xsec(29, 1.73, -0.6);
    EXPECT_SOFT_EQ(0.0289065, xsec());
}

TEST_F(CoulombScatteringTest, mott_xs)
{
    WokviHostRef const& data = model_->host_ref();

    real_type inc_energy
        = value_as<units::MevEnergy>(particle_track().energy());
    real_type inc_mass = value_as<units::MevMass>(particle_track().mass());
    WokviElementData const& element_data = data.elem_data[ElementId(0)];

    cout << "Energy: " << inc_energy << "\n"
         << "Mass: " << inc_mass << "\n";

    MottXsCalculator xsec(element_data, inc_energy, inc_mass);
    EXPECT_SOFT_EQ(0.837583, xsec(0.21));

    static double const cos_ts[] = {0.21, 0.5, 0.9, 0, -0.1, -0.6, -0.7};
    static double const expected_xsecs[]
        = {0.837776, 0.986799, 1.09074, 0.712007, 0.64827, 0.302596, 0.229266};
    std::vector<double> xsecs;
    for (double cos_t : cos_ts)
    {
        xsecs.push_back(xsec(cos_t));
    }

    EXPECT_VEC_SOFT_EQ(xsecs, expected_xsecs);
}

// TEST_F(CoulombScatteringTest, helper_state)
// {
//     detail::WokviStateHelper state(this->particle_track(),
//                                    this->material_track().make_material_view(),
//                                    ElementComponentId{0},
//                                    units::MevEnergy{10},
//                                    model_->host_ref());
//
//     // Check the incident particle quantities are correct
//     EXPECT_SOFT_EQ(state.inc_energy, 10.0);
//     EXPECT_SOFT_EQ(state.inc_mass, 0.51099894609999996);
//     EXPECT_SOFT_EQ(state.inc_mom_sq, 110.22);
//     EXPECT_SOFT_EQ(state.inv_beta_sq, 1.0023691);
//
//     // Check the element data is correct
//     EXPECT_EQ(state.element.atomic_number().get(), 29);
//     EXPECT_SOFT_EQ(state.target_Z(), 29);
//     EXPECT_SOFT_EQ(state.target_mass(), 63.546);
//     EXPECT_SOFT_EQ(state.element_data.mott_coeff[0][1], 3.24569e-5);
//
//     EXPECT_SOFT_EQ(state.kinetic_factor,
//                    value_as<WokviRef::CoeffQuantity>(model_->host_ref().coeff)
//                        * 29.0 * 1.0023691 / 110.22);
//     EXPECT_SOFT_EQ(state.mott_factor(), 1.1682);
//     EXPECT_SOFT_EQ(state.screening_coefficient(), 2.0);
//
//     // Check angle bounds
//     EXPECT_SOFT_EQ(state.cos_t_min_nuc(), 1.0);
//     EXPECT_SOFT_EQ(state.cos_t_max_nuc(), -1.0);
//     EXPECT_SOFT_EQ(state.max_electron_cos_t(), -1.0);
//     EXPECT_SOFT_EQ(state.cos_t_min_elec(), 1.0);
//     EXPECT_SOFT_EQ(state.cos_t_max_elec(), -1.0);
// }

//---------------------------------------------------------------------------//
}  // namespace test
}  // namespace celeritas
