#pragma once

#include "corecel/data/CollectionMirror.hh"
#include "corecel/data/ParamsDataInterface.hh"
#include "celeritas/grid/GenericGridBuilder.hh"

namespace celeritas
{

template<class IdType>
class ListBuilder
{
  public:
    template<class T>
    using Items = Collection<T, Ownership::value, MemSpace::host, IdType>;

    template<class ValueType, class SourceType, class BuilderType>
    using BuilderMethod = ValueType (BuilderType::*)(SourceType const&);

    template<class ValueType, class SourceType, class BuilderType>
    using ConstBuilderMethod
        = ValueType (BuilderType::*)(SourceType const&) const;

  public:
    ListBuilder(Range<IdType> ids) : ids_(ids) {}

    template<class ValueType, class SourceType, class BuilderType>
    Range<IdType>
    operator()(Items<ValueType>* data,
               Span<SourceType const> source,
               BuilderType const& builder,
               ConstBuilderMethod<ValueType, SourceType, BuilderType> method
               = &BuilderType::operator())
    {
        auto build_data = CollectionBuilder(data);
        IdType start_id{data->size()};
        IdType end_id = start_id;
        for (auto id : ids_)
        {
            end_id = build_data.push_back(
                ((&builder)->*method)(source[id.get()]));
        }
        return {start_id, end_id};
    }

    template<class ValueType, class SourceType, class BuilderType>
    Range<IdType>
    operator()(Items<ValueType>* data,
               Span<SourceType const> source,
               BuilderType& builder,
               BuilderMethod<ValueType, SourceType, BuilderType> method
               = &BuilderType::operator())
    {
        auto build_data = CollectionBuilder(data);
        IdType start_id{data->size()};
        IdType end_id = start_id;
        for (auto id : ids_)
        {
            end_id = build_data.push_back(
                ((&builder)->*method)(source[id.get()]));
        }
        return {start_id, end_id};
    }

    template<class ValueType>
    Range<IdType>
    operator()(Items<ValueType>* data, Span<ValueType const> source)
    {
        return (*this)(
            data, source, [](auto const& x) { return ValueType{x}; });
    }

  private:
    Range<IdType> ids_;
};

template<template<Ownership, MemSpace> class P>
class DataBuilder;

template<template<Ownership, MemSpace> class P>
class GridDataBuilder
{
  public:
    GridDataBuilder(HostVal<P>* host_data)
        : host_data_(host_data), build_grid_(&host_data_->reals)
    {
    }

  protected:
    GenericGridBuilder build_grid_;
    HostVal<P>* host_data_;
};

template<template<Ownership, MemSpace> class P>
class MirroredParamDataInterface : public ParamsDataInterface<P>
{
  public:
    using HostRef = typename ParamsDataInterface<P>::HostRef;
    using DeviceRef = typename ParamsDataInterface<P>::DeviceRef;

    HostRef const& host_ref() const override final { return data_.host_ref(); }
    DeviceRef const& device_ref() const override final
    {
        return data_.device_ref();
    }

  protected:
    template<class... Args>
    MirroredParamDataInterface(Args const&... args)
    {
        HostVal<P> host_data;

        DataBuilder<P>{&host_data}(args...);

        data_ = CollectionMirror<P>{std::move(host_data)};
        CELER_ENSURE(data_);
    }

    CollectionMirror<P> data_;
};

}  // namespace celeritas
