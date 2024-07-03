//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/CerenkovParams.cc
//---------------------------------------------------------------------------//
#include "CerenkovParams.hh"

#include <utility>
#include <vector>

#include "corecel/cont/Range.hh"
#include "corecel/data/CollectionBuilder.hh"
#include "corecel/data/DedupeCollectionBuilder.hh"
#include "corecel/math/Algorithms.hh"
#include "celeritas/Quantities.hh"
#include "celeritas/Types.hh"
#include "celeritas/grid/GenericGridInserter.hh"

#include "OpticalPropertyParams.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct with optical property data.
 */
template<>
class DataBuilder<CerenkovData> : public GridDataBuilder<CerenkovData>
{
  public:
    DataBuilder(HostVal<CerenkovData>* host_data)
        : GridDataBuilder<CerenkovData>(host_data)
    {
    }

    void operator()(CerenkovParams::SPConstProperties properties)
    {
        CELER_EXPECT(properties);

        auto const& host_ref = properties->host_ref();
        auto opt_materials
            = range(OpticalMaterialId{host_ref.refractive_index.size()});

        ListBuilder{opt_materials}(
            &host_data_->angle_integral,
            host_ref.refractive_index[opt_materials],
            [&](auto const& ri_grid) -> GenericGridData {
                if (!ri_grid)
                {
                    // No refractive index data stored for this material
                    return {};
                }

                return import_angle_integral(host_ref.reals[ri_grid.value],
                                             host_ref.reals[ri_grid.grid]);
            });

        CELER_ASSERT(host_data_->angle_integral.size()
                     == host_ref.refractive_index.size());
    }

  private:
    GenericGridData
    import_angle_integral(Span<real_type const> refractive_index,
                          Span<real_type const> energy)
    {
        std::vector<real_type> integral(energy.size());
        for (size_type i = 1; i < energy.size(); ++i)
        {
            integral[i] = integral[i - 1]
                          + real_type(0.5) * (energy[i] - energy[i - 1])
                                * (1 / ipow<2>(refractive_index[i - 1])
                                   + 1 / ipow<2>(refractive_index[i]));
        }

        return build_grid_(energy, make_span(integral));
    }
};

CerenkovParams::CerenkovParams(SPConstProperties properties)
    : MirroredParamDataInterface<CerenkovData>(properties)
{
    CELER_ENSURE(data_ || properties->host_ref().refractive_index.empty());
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
