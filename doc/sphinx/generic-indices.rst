Index manipulation
==================

Defined in header ``<vif/utility/generic.hpp>``.


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

This function converts a "flat" index ``i`` into an array of multidimensional indices, following the provided dimensions ``dims`` ([1] and [2]) or the provided vector ``v`` ([3] and [4]). The ``flat_id`` function does the inverse job.

[2] and [4] are the vectorized version of [1] and [3], respectively. The return value is a 2D vector of indices: the first dimension contains as many elements as ``dims`` ([2]) or the dimensions of ``v`` ([4]), and the second dimension contains as many elements as the provided index vector ``i``.

**Example:**

.. code-block:: c++

    vec2i v(2,3);
    mult_ids(v,0); // {0,0}
    mult_ids(v,1); // {0,1}
    mult_ids(v,2); // {0,3}
    mult_ids(v,3); // {1,0}
    v[3] == v(1,0); // true


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

This function converts a set of multidimensional indices into a "flat" index, following the provided dimensions ``dims`` ([1] and [2]) or the provided vector ``v`` ([3] and [4]). The ``mult_ids`` function does the inverse job.

In [1] and [3], the multidimensional indices are provided as separate arguments to the function. In [2] and [4] they are grouped inside an index vector.

**Example:**

.. code-block:: c++

    vec2i v(2,3);
    flat_id(v,0,0); // 0
    flat_id(v,0,1); // 1
    flat_id(v,0,2); // 2
    flat_id(v,1,0); // 3
    v(1,0) == v[3]; // true


increment_index_list
--------------------

.. code-block:: c++

    void increment_index_list(vec1u& ids, const uint_t& n); // [1]

    void increment_index_list(vec1u& ids, const vec1u& n); // [2]

These functions perform *one* increment on the set of multidimensional indices ``ids``, following the order in memory (i.e., last dimension is incremented first). In [1], each dimension has the same size ``n``, while in [2] the dimensions may be different and are provided as a vector ``n``. These functions allow the full traversal of the multidimensional space in a single loop, and are typically used to iterate on a multidimensional data set when the number of dimensions is either too large or unknown at compile time.

If called on the last allowed index, the function will set ``ids`` to zero, hence come back to the first index.


**Example:**

.. code-block:: c++

    // Say we got some multidimensional data from a file
    vec1d data;
    vec1u dims = /* read from a file */;

    uint_t nelem = 1;
    for (uint_t d : dims) nelem *= d;

    data.resize(nelem);

    // Initialize the index vector to zero (first index)
    vec1u ids(dims.size());

    // Iterate in one single loop
    for (uint_t i : range(nelem)) {
        // data[i] is the element at index (ids[0],ids[1],...)

        // Increment using [2]
        increment_index_list(ids, dims);
    }
