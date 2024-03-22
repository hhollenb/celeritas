//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/optical/ImportedOpticalProcessAdapter.hh
//---------------------------------------------------------------------------//
#pragma once

#include "celeritas/io/ImportProcess.hh"
#include "celeritas/io/ImportPhysicsTable.hh"

namespace celeritas
{

class IOPAContextException : public RichContextException
{
  public:
    IOPAContextException(ImportProcessClass ipc, MaterialId mid);

    //! This class type
    char const* type() const final { return "ImportOpticalProcessAdapterContext"; }

    //! Save context to a JSON object
    void output(JsonPimpl*) const final {}

    //! Get an explanatory message
    char const* what() const noexcept final { return what_.c_str(); }

  private:
    std::string what_;
};

struct ImportOpticalProcess
{
    const ImportProcessType process_type{ImportProcessType::optical};
    ImportProcessClass process_class{ImportProcessClass::size_};
    ImportPhysicsTable lambda_table;

    explicit operator bool() const
    {
        return process_type == ImportProcessType::optical
            && process_class != ImportProcessClass::size_
            && lambda_table.table_type == ImportTableType::lambda
            && lambda_table;
};

class ImportedOpticalProcesses
{
  public:
    //!@{
    //! \name Type aliases
    using ImportedOpticalProcessId = OpaqueId<ImportOpticalProcess>;
    using key_type = ImportProcessClass;
    //!@}

  public:
    // Construct with imported data
    static std::shared_ptr<ImportedOpticalProcesses>
    from_import(ImportData const& data);

    // Construct with imported tables
    explicit ImportedOpticalProcesses(std::vector<ImportOpticalProcess> io);

    // Return the process ID for the given process class
    ImportOpticalProcessId find(key_type) const;

    // Get the table for the given process ID
    inline ImportOpticalProcess const& get(ImportOpticalProcessId id) const;

    // Number of imported optical processes
    inline ImportOpticalProcessId::size_type size() const;

  private:
    std::vector<ImportOpticalProcess> processes_;
    std::map<key_type, ImportOpticalProcessId> ids_;
};

//---------------------------------------------------------------------------//
/*!
 *
 */
class ImportedOpticalProcessAdapter
{
  public:
    //!@{
    //! \name Type aliases
    using SPConstImported = std::shared_ptr<ImportedOpticalProcesses const>;
    using GridBuilder = OpticalProcess::GridBuilder;
    using StepLimitBuilder = std::unique_ptr<GridBuilder const>;
    //!@}

  public:
    //! Construct from shared table data
    ImportedOpticalProcessAdapter(SPConstImported imported,
                                  ImportProcessClass process_class);

    //! Construct step limits for the process
    std::vector<OpticalValueGridId> step_limits(GenericGridInserter&, MaterialParams const&) const;

    //! Get the lambda table for the process
    inline ImportPhysicsTable const& get_lambda() const;

    //! Access the imported process
    inline ImportOpticalProcess const& process() const;

  private:
    SPConstImported imported_;
    ImportProcessClass process_class_;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
ImportOpticalProcess const& ImportOpticalProcesses::get(ImportOpticalProcessId id) const
{
    CELER_EXPECT(id < this->size());
    return processes_[id.get()];
}

ImportOpticalProcessId::size_type ImportOpticalProcesses::size() const
{
    return processes_.size();
}

ImportPhysicsTable const& ImportOpticalProcessAdapter::get_lambda() const
{
    return process().lambda_table;
}

ImportOpticalProcess const& ImportOpticalProcessAdapter::process() const
{
    return imported_->get(imported_->find(process_class));
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
