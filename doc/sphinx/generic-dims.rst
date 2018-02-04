Modifying dimensions
====================

Defined in header ``<phypp/utility/generic.hpp>``.

flatten
-------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(vec<Dim,Type>&& v); // [2]

    template<std::size_t Dim, typename Type>
    vec<1,Type*> flatten(const vec<Dim,Type*>& v); // [3]

    template<std::size_t Dim, typename Type>
    vec<1,Type*> flatten(vec<Dim,Type*>&& v); // [4]

    template<typename T>
    T flatten(T&& t); // [5]


reform
------

.. code-block:: c++

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type> reform(const vec<Dim,Type>& v, Args&& ... args); // [1]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type> reform(vec<Dim,Type>&& v, Args&& ... args); // [2]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type*> reform(const vec<Dim,Type*>& v, Args&& ... args); // [3]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type*> reform(vec<Dim,Type*>&& v, Args&& ... args); // [4]
