//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/OpticalModelBuilder.hh
//---------------------------------------------------------------------------//
#pragma once

#include <functional>
#include <memory>
#include <set>
#include <unordered_map>

#include "corecel/cont/Range.hh"
#include "celeritas/Types.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
class ImportData;
class ImportedOpticalMaterials;
class OpticalModel;

//---------------------------------------------------------------------------//
/*!
 * Enum of types of optical models that may be imported.
 */
enum class ImportOpticalModelClass
{
    other,
    // Optical
    absorption,
    rayleigh,
    wavelength_shifting,
    size_
};

//---------------------------------------------------------------------------//
/*!
 * Options used for building optical models.
 */
struct OpticalModelBuilderOptions
{
};

//---------------------------------------------------------------------------//
/*!
 * Construct Celeritas optical models from imported data.
 *
 * This is a factory class which constructs optical models from both user
 * models and built-in models, based on the provided user build map. If no
 * map is provided, then built-in optical models are used by default.
 */
class OpticalModelBuilder
{
  public:
    //!@{
    //! \name Type aliases
    using IOMC = ImportOpticalModelClass;
    using Options = OpticalModelBuilderOptions;
    using SPModel = std::shared_ptr<OpticalModel>;
    using SPConstImported = std::shared_ptr<ImportedOpticalMaterials const>;
    using ActionIdIter = RangeIter<ActionId>;
    //!@}

    //! Input argument for user-provided optical model construction
    struct UserBuildInput
    {
        SPConstImported imported;
    };

    //!@{
    //! \name User builder type aliases
    using UserBuildFunction
        = std::function<SPModel(ActionIdIter, UserBuildInput const&)>;
    using UserBuildMap = std::unordered_map<IOMC, UserBuildFunction>;
    //!@}

  public:
    //! Get an ordered set of all available optical models
    static std::set<IOMC> get_all_model_classes();

    //! Construct from imported and shared data
    OpticalModelBuilder(ImportData const& data,
                        UserBuildMap user_build,
                        Options options);

    //! Construct from imported and shared data without custom user builders
    OpticalModelBuilder(ImportData const& data,
                        Options options);

    //! Create an optical model from the data
    SPModel operator()(IOMC iomc, ActionIdIter start_id) const;

  private:
    UserBuildInput input_;
    UserBuildMap user_build_map_;

    SPConstImported const imported() const { return input_.imported; }

    //!@{
    //! \name Helper functions to build specific models
    SPModel build_absorption(ActionId) const;
    SPModel build_rayleigh(ActionId) const;
    // SPModel build_wls() const;
    //!@}
};

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//

//! Get the string form of the ImportOpticalModelClass enumeration
char const* to_cstring(ImportOpticalModelClass value);

//---------------------------------------------------------------------------//
}  // namespace celeritas
