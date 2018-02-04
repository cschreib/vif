Range-based iteration
=====================

In ``phypp/core/range.hpp``.

range
-----

.. code-block:: c++

    template<std::size_t D, typename T>
    /*...*/ range(const vec<D,T>& v); // 1

    /*...*/ range(uint_t n); // 2

    /*...*/ range(uint_t i0, uint_t n); // 3
