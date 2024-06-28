//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file corecel/cont/Registry.hh
//---------------------------------------------------------------------------//
#pragma once

namespace celeritas
{
// Forward declarations
template <class Value>            class RegistryAdapter;
template <class Value, class Key> class MappedRegistryAdapter;

//---------------------------------------------------------------------------//
/*!
 * Type-safe registry for uniquely assigning opaque IDs to a set of values.
 */
template <class Value>
class Registry
{
  public:
    //!@{
    //! \name Type aliases
    using ValueId = OpaqueId<Value>;
    //!@}

  public:
    //! Create a registry from an existing set of values
    explicit Registry(std::vector<Value> values);

    //! Retrieve a value from the registry via an ID
    inline Value const& get(ValueId id) const;

    //! Number of values in the repository
    inline ValueId::size_type size() const;

    //! Get a range of value IDs for this registry
    inline Range<ValueId> id_range() const;

  protected:
    std::vector<Value> values_;
};

//---------------------------------------------------------------------------//
/*!
 * Type-safe registry that further associates a unique key with the opaque
 * IDs of a set of values.
 *
 * Keys are usually enums (e.g. ImportedProcessClass), and each value should
 * have a unique key already associated with it.
 */
template <class Value, class Key>
class MappedRegistry : public Registry<Value>
{
  public:
    //!@{
    //! \name Type aliases
    using MapAdapter = MappedValueAdapater<Value, Key>;
    //!@}

  public:
    //! Create a key-mapped registry for the specified values and key map
    template <class KeyMap>
    explicit MappedRegistry(std::vector<Value> values, KeyMap const& keymap);

    //! Return the value ID for the given key
    ValueId find(Key key) const;

  protected:
    std::map<Key, ValueId> id_map_;
};

//---------------------------------------------------------------------------//
/*!
 * An adapter for a specific ID in a registry.
 *
 * This lightweight type acts as a view for the value associated with this ID,
 * without referring to an explicit memory location. 
 */
template <class Value>
class RegistryAdapter
{
  public:
    //!@{
    //! \name Type aliases
    using SPConstRegistry = std::shared_ptr<Registry<Value> const>;
    using ValueId = OpaqueId<Value>;
    //!@}

  public:
    //! Construct adapter directly from a registry and an ID
    RegistryAdapter(SPConstRegistry registry, ValueId value_id);

    //! Get the value associated with the ID
    inline Value const& get() const;

    //! Get the value associated with the ID
    Value const& operator*() const;

    //! Get the value associated with the ID as a pointer
    Value const* operator->() const;

  private:
    SPConstRegistry registry_;
    ValueId value_id_;
};

//---------------------------------------------------------------------------//
/*!
 * An adapter for a specific key in a mapped registry.
 *
 * This lightweight type acts as a view for the value associated with this key,
 * without referring to an explicit memory location.
 */
template <class Value, class Key>
class MappedRegistryAdapter
{
  public:
    //!@{
    //! \name Type aliases
    using SPConstRegistry = std::shared_ptr<MappedRegistry<Value, Key> const>;
    //!@}

  public:
    //! Construct adapter directly from a registry and an ID
    MappedRegistryAdapter(SPConstRegistry registry, Key key);

    //! Get the value associated with the key
    inline Value const& get() const;

    //! Get the value associated with the key
    Value const& operator*() const;

    //! Get the value associated with the key as a pointer
    Value const* operator->() const;

  private:
    SPConstRegistry registry_;
    Key key;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct a registry from the given values.
 */
template <class Value>
Registry<Value>::Registry(std::vector<Value> values)
    : values_(values)
{}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value assoicated with the given ID.
 */
template <class Value>
auto Registry<Value>::get(ValueId id) const -> Value const&
{
    CELER_EXPECT(id < this->size());
    return values_[id.get()];
}

//---------------------------------------------------------------------------//
/*!
 * Return the number of values in the registry.
 */
template <class Value>
auto Registry<Value>::size() const -> ValueId::size_type
{
    return values_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Return a range of IDs for the values in this registry.
 */
template <class Value>
auto Registry<Value>::id_range() const -> Range<ValueId>
{
    return range(ValueId{values_.size()});
}

//---------------------------------------------------------------------------//
/*!
 * Construct a mapped registry from the specified values.
 *
 * The specified key map should return a unique key associated with a provided
 * value.
 */
template <class Value, class Key>
template <class KeyMap>
MappedRegistry<Value, Key>::MappedRegistry<KeyMap>(std::vector<Value> values, KeyMap const& keymap)
    : Registry<Value>(values), id_map_()
{
    for (ValueId id : id_range())
    {
        Key key = keymap(values_[id.get()]);

        auto insertion = id_map_.insert({key, id});
        CELER_VALIDATE(insertion.second, << "encountered duplicate key when constructing registry.");
    }

    CELER_ENSURE(values_.size() == id_map_.size());
}

//---------------------------------------------------------------------------//
/*!
 * Construct a registry adapter from the given registry and ID.
 */
template <class Value>
RegistryAdapter<Value>::RegistryAdapter(SPConstRegistry registry, ValueId value_id)
    : registry_(registry), value_id_(valud_id)
{}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's ID.
 */
template <class Value>
auto RegistryAdapter<Value>::get() const -> Value const&
{
    return registry_->get(value_id_);
}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's ID.
 */
template <class Value>
auto RegistryAdapter<Value>::operator*() const -> Value const&
{
    return get();
}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's ID as a pointer.
 */
template <class Value>
auto RegistryAdapter<Value>::operator->() const -> Value const*
{
    return &(get());
}

//---------------------------------------------------------------------------//
/*!
 * Construct a mapped registry adapter from the given registry and key.
 */
template <class Value, class Key>
auto MappedRegistryAdapter::MappedRegistryAdapter(SPConstRegistry registry, Key key)
    : registry_(registry), key_(key)
{}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's key.
 */
template <class Value, class Key>
auto MappedRegistryAdapter<Value, Key>::get() const -> Value const&
{
    return registry_->get(key);
}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's key.
 */
template <class Value, class Key>
auto MappedRegistryAdapter<Value, Key>::operator*() const -> Value const&
{
    return get();
}

//---------------------------------------------------------------------------//
/*!
 * Retrieve the value associated with the adaptor's key as a pointer.
 */
template <class Value, class Key>
auto MappedRegistryAdapter<Value, Key>::operator->() const -> Value const*
{
    return &(get());
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
