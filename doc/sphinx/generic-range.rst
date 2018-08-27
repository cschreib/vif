Range-based iteration
=====================

Defined in header ``<vif/core/range.hpp>``.

range
-----

.. code-block:: c++

    template<std::size_t D, typename T>
    /*...*/ range(const vec<D,T>& v); // [1]

    /*...*/ range(uint_t n); // [2]

    /*...*/ range(uint_t i0, uint_t n); // [3]

This function returns a C++ *range*, that is, an object that can be used inside the C++ range-based ``for`` loop. This range will generate integer values starting from ``0`` (in [1], [2]) or ``i0`` (in [3]) to ``v.size()`` (in [1]) or ``n`` (in [2], [3]), that last value being *excluded* from the range. This nice way of writing an integer ``for`` loop actually runs as fast as (if not faster than) the classical way, and is less error prone.

The return value is a proxy class that holds the starting and ending point of the range, and offers ``begin()`` and ``end()`` function for iteration. Its type is of little importance.

**Example:**

.. code-block:: c++

    vec1i v = {4,5,6,8};

    // First version
    for (uint_t i : range(v)) { // [1]
        // 'i' goes from 0 to 3
        v[i] = ...;
    }

    // Note that the loop above generates
    // *indices* inside the vector, while:
    for (int i : v) { /* ... */ }
    // ... generates *values* from the vector.

    // Second version
    for (uint_t i : range(3)) { // [2]
        // 'i' goes from 0 to 2
        v[i] = ...;
    }

    // Third version
    for (uint_t i : range(1,3)) { // [3]
        // 'i' goes from 1 to 3
        v[i] = ...;
    }
