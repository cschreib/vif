Integer sequences
=================

Defined in header ``<phypp/utility/generic.hpp>``.

uindgen, indgen, findgen, dindgen
---------------------------------

.. code-block:: c++

    template<typename ... Dims>
    vec</*...*/,uint_t> uindgen(Dims&& ... ds); // [1]

    template<typename ... Dims>
    vec</*...*/,int_t> indgen(Dims&& ... ds); // [2]

    template<typename ... Dims>
    vec</*...*/,float> findgen(Dims&& ... ds); // [3]

    template<typename ... Dims>
    vec</*...*/,double> dindgen(Dims&& ... ds); // [4]

These functions will create a new vector with values starting at ``0`` and increment linearly by steps of ``1`` until the end of the vector. The number of dimensions of the resulting vector depends on the types ``Args`` of the arguments:

* The starting dimension is ``0``.
* Each argument of type ``uint_t`` increases the final dimension by one.
* Each argument of type ``std::array<uint_t,D>`` increases the final dimension by ``D``.

**Example:**

.. code-block:: c++

    vec1i v = indgen(5);    // {0,1,2,3,4}
    vec2u w = uindgen(3,2); // {{0,1}, {2,3}, {4,5}}
