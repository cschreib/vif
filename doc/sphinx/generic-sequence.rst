Integer sequences
=================

In ``phypp/utility/generic.hpp``.

uindgen, indgen, findgen, dindgen
---------------------------------

.. code-block:: c++

    template<typename ... Dims>
    vec<D,uint_t> uindgen(Dims&& ... ds); // 1

    template<typename ... Dims>
    vec<D,int_t> indgen(Dims&& ... ds); // 2

    template<typename ... Dims>
    vec<D,float> findgen(Dims&& ... ds); // 3

    template<typename ... Dims>
    vec<D,double> dindgen(Dims&& ... ds); // 4

