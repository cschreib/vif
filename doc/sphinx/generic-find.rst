Finding elements
================

Defined in header ``<phypp/utility/generic.hpp>``.

where, where_first, where_last
------------------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u where(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type>
    uint_t where_first(const vec<Dim,Type>& v); // [2]

    template<std::size_t Dim, typename Type>
    uint_t where_last(const vec<Dim,Type>& v); // [3]


complement
----------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u complement(const vec<Dim,Type>& v, const vec1u& ids);


match
-----

.. code-block:: c++

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    void match(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2, vec1u& id1, vec1u& id2);


intersection_set, union_set
---------------------------

.. code-block:: c++

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> intersection_set(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2);

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> union_set(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2);


unique_ids, unique_values
-------------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u unique_ids(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type>
    vec1u unique_ids(const vec<Dim,Type>& v, const vec1u& sid); // [2]

    template<std::size_t Dim, typename Type>
    vec1u unique_ids_sorted(const vec<Dim,Type>& v); // [3]

    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values(const vec<Dim,Type>& v); // [4]

    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values(const vec<Dim,Type>& v, const vec1u& sid); // [5]

    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values_sorted(const vec<Dim,Type>& v); // [6]


is_any_of
---------

.. code-block:: c++

    template<typename Type1, std::size_t Dim2, typename Type2>
    bool is_any_of(const Type1& v1, const vec<Dim2,Type2>& v2); // [1]

    template<std::size_t Dim1, typename Type1, std::size_t Dim2 = Dim1, typename Type2>
    vec<Dim1,bool> is_any_of(const vec<Dim1,Type1>& v1, const vec<Dim2,Type2>& v2); // [2]


bounds, lower_bound, upper_bound
--------------------------------

.. code-block:: c++

    template<typename T, std::size_t Dim, typename Type>
    uint_t lower_bound(const vec<Dim,Type>& v, T x); // [1]

    template<typename T, std::size_t Dim, typename Type>
    uint_t upper_bound(const vec<Dim,Type>& v, T x); // [2]

    template<typename T, std::size_t Dim, typename Type>
    std::array<uint_t,2> bounds(const vec<Dim,Type>& v, T x); // [3]

    template<typename T, typename U, std::size_t Dim, typename Type>
    std::array<uint_t,2> bounds(const vec<Dim,Type>& v, T x1, U x2); // [4]


equal_range
-----------

.. code-block:: c++

    template<typename T, std::size_t Dim, typename Type>
    vec1u equal_range(const vec<Dim,Type>& v, T x);


astar_find
----------

.. code-block:: c++

    bool astar_find(const vec2b& map, uint_t& x, uint_t& y);
