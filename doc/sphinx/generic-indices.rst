Index manipulation
==================

Defined in header ``<phypp/utility/generic.hpp>``.


mult_ids
--------

.. code-block:: c++

    template<std::size_t D>
    vec1u mult_ids(const std::array<uint_t,D>& dims, uint_t i); // [1]

    template<std::size_t D>
    vec2u mult_ids(const std::array<uint_t,D>& dims, vec1u i); // [2]

    template<std::size_t D, typename T>
    vec1u mult_ids(const vec<D,T>& v, uint_t i); // [3]

    template<std::size_t D, typename T>
    vec2u mult_ids(const vec<D,T>& v, vec1u i); // [4]


flat_id
-------

.. code-block:: c++

    template<std::size_t D, typename ... Args>
    uint_t flat_id(const std::array<uint_t,D>& dims, Args&& ... args); // [1]

    template<std::size_t D, typename TI>
    uint_t flat_id(const std::array<uint_t,D>& dims, const vec<1,TI>& ids); // [2]

    template<std::size_t D, typename T, typename ... Args>
    uint_t flat_id(const vec<D,T>& v, Args&& ... args); // [3]

    template<std::size_t D, typename T, typename TI>
    uint_t flat_id(const vec<D,T>& v, const vec<1,TI>& ids); // [4]


increment_index_list
--------------------

.. code-block:: c++

    void increment_index_list(vec1u& ids, const uint_t& n); // [1]

    void increment_index_list(vec1u& ids, const vec1u& n); // [2]
