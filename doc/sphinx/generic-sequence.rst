Integer sequences
=================

Defined in header ``<vif/utility/generic.hpp>``.

indgen
------

.. code-block:: c++

    template<typename T = uint_t, typename ... Dims>
    vec</*...*/,T> indgen(Dims&& ... ds);

This functions will create a new vector with values starting at ``0`` and increment linearly by steps of ``1`` until the end of the vector. Internally, the values are generated with the standard function ``std::iota```. The number of dimensions of the resulting vector depends on the types ``Args`` of the arguments:

* Each argument of type ``uint_t`` increases the number of dimensions by one.
* Each argument of type ``std::array<uint_t,D>`` increases the number of dimensions by ``D``.

The type of the values in the resulting vector is determined by the template parameter ``T``, which defaults to ``uint_t`` if none is provided.

**Example:**

.. code-block:: c++

    vec1u v = indgen(5);   // {0,1,2,3,4}
    vec2u w = indgen(3,2); // {{0,1}, {2,3}, {4,5}}
