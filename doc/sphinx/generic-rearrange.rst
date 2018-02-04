Rearranging elements
====================

Defined in header ``<phypp/utility/generic.hpp>``.

reverse
-------

.. code-block:: c++

    template<typename Type>
    vec<1,Type> reverse(vec<1,Type> v);


shift, inplace_shift
--------------------

.. code-block:: c++

    template<typename Type>
    vec<1,Type> shift(vec<1,Type> v, int_t n); // [1]

    template<typename Type>
    void inplace_shift(vec<1,Type>& v, int_t n); // [2]


transpose
---------

.. code-block:: c++

    template<typename Type>
    vec<2,Type> transpose(const vec<2,Type>& v);


replicate
---------

.. code-block:: c++

    template<typename Type, typename ... Args>
    vec</*...*/, meta::vtype_t<Type>> replicate(const Type& t, Args&& ... args); // [1]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/, meta::rtype_t<Type>> replicate(const vec<Dim,Type>& t, Args&& ... args); // [2]


sort, inplace_sort
------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u sort(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type, typename F>
    vec1u sort(const vec<Dim,Type>& v, F&& comp); // [2]

    template<std::size_t Dim, typename Type>
    void inplace_sort(vec<Dim,Type>& v); // [3]

    template<std::size_t Dim, typename Type, typename F>
    void inplace_sort(vec<Dim,Type>& v, F&& comp); // [4]

is_sorted
---------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    bool is_sorted(const vec<Dim,Type>& v);


append, prepend
---------------

.. code-block:: c++

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2>
    void append(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2); // [1]

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2>
    void prepend(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2); // [2]


remove, inplace_remove
----------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec<Dim,Type> remove(vec<Dim,Type> v, const vec1u& ids); // [1]

    template<std::size_t Dim, typename Type>
    void inplace_remove(vec<Dim,Type>& v, vec1u ids); // [2]
