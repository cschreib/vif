Rearranging elements and dimensions
===================================

Defined in header ``<vif/utility/generic.hpp>``.

flatten
-------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type>
    vec<1,Type> flatten(vec<Dim,Type>&& v); // [2]

    template<typename T>
    T flatten(T&& v); // [3]

This function transforms a multidimensional vector into a 1D vector ([1] and [2]; [3] is a no-op overload for scalars). The content in memory remains exactly the same, so the operation is fast. In particular, if the argument of this function is a temporary ([2]), this function is extremely cheap as it produces no copy. The ``reform()`` function does the inverse job, and more.

The provided argument ``v`` is unchanged ([1] and [3]).

**Example:**

.. code-block:: c++

    vec2i v = {{1,2,3}, {4,5,6}};
    vec1i w = flatten(v); // {1,2,3,4,5,6}


reform
------

.. code-block:: c++

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type> reform(const vec<Dim,Type>& v, Args&& ... args); // [1]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/,Type> reform(vec<Dim,Type>&& v, Args&& ... args); // [2]

This function transforms a vector into another vector, simply changing its dimensions. The content in memory remains exactly the same, so the operation is fast. In particular, if the argument of this function is a temporary ([2]), this function is extremely cheap as it produces no copy. However, the new dimensions have to sum up to the same number of elements as that in the provided vector. The ``flatten()`` function is a special case of ``reform()`` where all dimension are reformed into one, resulting in a 1D vector.

The provided argument ``v`` is unchanged in [1], but not [2]. The number of dimensions of the resulting vector depends on the types ``Args`` of the arguments:

* The starting dimension is ``0``.
* Each argument of type ``uint_t`` increases the final dimension by one.
* Each argument of type ``std::array<uint_t,D>`` increases the final dimension by ``D``.

**Example:**

.. code-block:: c++

    vec1i v = {1,2,3,4,5,6};
    vec2i w = reform(v, 2, 3); // {{1,2,3}, {4,5,6}}


reverse
-------

.. code-block:: c++

    template<typename Type>
    vec<1,Type> reverse(vec<1,Type> v);

This function will return a copy of the provided vector, in which the order of all the elements is reversed. The original vector is unchanged. Only works with 1D vectors or views.

**Example:**

.. code-block:: c++

    vec1i v = {1,2,3,4,5,6};
    vec1i w = reverse(v); // {6,5,4,3,2,1}


shift, inplace_shift
--------------------

.. code-block:: c++

    template<typename Type>
    vec<1,Type> shift(vec<1,Type> v, int_t n); // [1]

    template<typename Type>
    void inplace_shift(vec<1,Type>& v, int_t n); // [2]

``shift()`` ([1]) returns a copy of the provided vector ``v`` where the elements are moved by circular shift of ``n`` elements. If ``n`` is positive, elements that would go beyond the bounds of the vector after the shift are moved to the beginning, with their order preserved. If ``n`` is negative, elements that would go beyond the beginning of the vector are placed at the end, with their order preserved. This function calls ``std::rotate()``. The original vector is unchanged. Only works with 1D vectors or views.

``inplace_shift()`` ([2]) performs the same operation as ``shift()`` but operates directly on the provided vector, which is therefore modified, but no copy is made so the operation is faster.

**Example:**

.. code-block:: c++

    vec1i v = {1,2,3,4,5};

    // [1]
    vec1i sr1 = shift(v, 2);  // {4,5,1,2,3}
    vec1i sr2 = shift(v, -2); // {3,4,5,1,2}

    // [2]
    inplace_shift(v, 2);
    // v = {4,5,1,2,3}


transpose
---------

.. code-block:: c++

    template<typename Type>
    vec<2,Type> transpose(const vec<2,Type>& v);

This function will transpose the provided 2D vector so that its dimensions are swapped. In other words, ``v(i,j)`` becomes ``v(j,i)``. This is a matrix transposition. The original vector is unchanged.

**Example:**

.. code-block:: c++

    vec2i v = {{1,2}, {3,4}, {5,6}};
    vec2i w = transpose(v); // {{1,3,5}, {2,4,6}}
    // now w(i,j) == v(j,i)


replicate
---------

.. code-block:: c++

    template<typename Type, typename ... Args>
    vec</*...*/, meta::vtype_t<Type>> replicate(const Type& t, Args&& ... args); // [1]

    template<std::size_t Dim, typename Type, typename ... Args>
    vec</*...*/, meta::rtype_t<Type>> replicate(const vec<Dim,Type>& t, Args&& ... args); // [2]

This function will take the provided scalar ([1]) or vector ([2]), and replicate it multiple times according to the provided additional parameters, to generate additional dimensions.

The number of dimensions of the resulting vector depends on the types ``Args`` of the arguments:

* The starting dimension is ``0`` ([1]) or ``Dim`` ([2]).
* Each argument of type ``uint_t`` increases the final dimension by one.
* Each argument of type ``std::array<uint_t,D>`` increases the final dimension by ``D``.

**Example:**

.. code-block:: c++

    // [1]
    vec1i v = replicate(2, 5);
    // v = {2,2,2,2,2}, or 5 times 2

    vec2i w = replicate(2, 3, 2);
    // w = {{2,2},{2,2},{2,2}}, or 3 x 2 times 2

    vec3u x = replicate(1u, w.dims, 5);
    // equivalent to:
    // x = replicate(1u, 3, 2, 5);

    // [2]
    vec2i z = replicate(vec1i{1,2}, 3);
    // z = {{1,2},{1,2},{1,2}}, or 3 times {1,2}

    // Note that it is not possible to just use a plain initializer list
    // since its type cannot be deduced with current C++ rules
    vec2i z = replicate({1,2}, 3); // error


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

``sort()`` returns a vector of indices for the provided vector ``v``, ordered such that the pointed values are sorted by increasing value ([1]) or following the provided comparison function ([2]). The number of returned indices is the same as the number of values in ``v``. The original vector is not modified. If two elements of ``v`` compare equal, their respective order in the vector will be unchanged (this function uses ``std::stable_sort()``).

``inplace_sort()`` directly modifies the order of the values inside the vector, and returns nothing. It is fastest, but less powerful.

In [2] and [4], the comparator function ``comp(x,y)`` must return ``true`` if ``x`` should be placed after ``y`` after the sort.

.. warning:: The comparison function ``comp`` must provide a *strict total ordering*, otherwise the behavior of the function is undefined. See `cppreference.com <http://en.cppreference.com/w/cpp/concept/Compare>`_ for more information. In brief, this means that any value can only be "equal", "lesser", or "greater" than any other value. With a comparison function returning simply ``x < y``, this requirement is not met for ``float`` and ``double`` because of the special value "not-a-number", ``NaN``, which is neither. [1] and [3] use the default comparator for vif vectors, in which this issue is solved by considering ``NaN`` as "greater than" positive infinity. ``NaN`` values will thus be placed at the end of a sorted vector. To take advantage of this implementation, use ``vec<Dim,Type>::comparator_less{}(x,y)`` and ``vec<Dim,Type>::comparator_greater{}(x,y)`` instead of ``x < y`` and ``x > y`` inside your custom comparison functions. This is unnecessary for integer types and strings.

**Example:**

.. code-block:: c++

    // [1]
    vec1i v = {1,5,6,3,7};
    vec1u id = sort(v); // {0,3,1,2,4}
    // v[id] = {1,3,5,6,7} is sorted

    // Now, 'id' can also be used to modify the order of
    // another vector of the same dimensions.

    // [3]
    inplace_sort(v);
    v; // {1,3,5,6,7} is sorted

    // [4]
    vec1f v1 = {1.0,2.0,3.0,4.0, 5.0,6.0};
    vec1f v2 = {3.0,0.5,1.0,fnan,0.0,0.0};

    // Sort 'v1+v2'
    vec1u id = uindgen(v1.size());
    inplace_sort(id, [&](uint_t i1, uint_t i2) {
        return vec1f::comparator_less{}(v1[i1]+v2[i1], v1[i2]+v2[i2]);
    });

    // (v1+v2)[id] = {2.5,4,4,5,6,nan}
    // v1[id]      = {2.0,1,3,5,6,4}
    // v2[id]      = {0.5,3,1,0,0,nan}

    // Sort first by 'v2', then 'v1'
    id = uindgen(v1.size());
    inplace_sort(id, [&](uint_t i1, uint_t i2) {
        if (vec1f::comparator_less{}(v2[i1], v2[i2])) {
            return true;
        } else if (vec1f::comparator_greater{}(v2[i1], v2[i2])) {
            return false;
        } else {
            return vec1f::comparator_less{}(v1[i1], v1[i2]);
        }
    });

    // v1[id] = {5,6,2.0,3,1,4}
    // v2[id] = {0,0,0.5,1,3,nan}


is_sorted
---------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    bool is_sorted(const vec<Dim,Type>& v);

This function just traverses the whole input vector and checks if its elements are sorted by increasing value.

**Example:**

.. code-block:: c++

    // First version
    vec1i v = {1,5,6,3,7};
    is_sorted(v); // false
    inplace_sort(v);
    // v = {1,3,5,6,7}
    is_sorted(v); // true


append, prepend
---------------

.. code-block:: c++

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2>
    void append(vec<Dim,Type1>& v, const vec<Dim,Type2>& t); // [1]

    template<std::size_t N, std::size_t Dim, typename Type1, typename Type2>
    void prepend(vec<Dim,Type1>& v, const vec<Dim,Type2>& t); // [2]

These functions behave similarly to ``vec::push_back()``, in that they will add new elements at the end ([1]), but also at the beginning ([2]) of the provided vector ``v``. However, while ``vec::push_back()`` can only add new elements from a vector that is one dimension *less* than the original vector (or a scalar, for 1D vectors), these functions will add new elements from a vector of the *same* dimension. These functions are also more powerful than ``vec::push_back``, because they allow you to choose along which dimension the new elements will be added using the template parameter ``N`` (note that this parameter is useless and therefore does not exist for 1D vectors). The other dimensions must be otherwise identical.

The first argument ``v`` cannot be a view.

**Example:**

.. code-block:: c++

    // For 1D vectors
    vec1i v = {1,2,3};
    vec1i w = {4,5,6};
    append(v, w);
    // v = {1,2,3,4,5,6}
    prepend(v, w);
    // v = {4,5,6,1,2,3,4,5,6}

    // For multidimensional vectors
    vec2i x = {{1,2}, {3,4}};          // x is (2x2)
    vec2i y = {{0}, {0}};              // y is (2x1)
    vec2i z = {{5,6,7}};               // z is (1x3)
    append<1>(x, y);
    // x = {{1,2,0}, {3,4,0}}          // x is (2x3)
    prepend<0>(x, z);
    // x = {{5,6,7}, {1,2,0}, {3,4,0}} // x is (3x3)


remove, inplace_remove
----------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec<Dim,Type> remove(vec<Dim,Type> v, const vec1u& ids); // [1]

    template<std::size_t Dim, typename Type>
    void inplace_remove(vec<Dim,Type>& v, vec1u ids); // [2]

``remove()`` ([1]) will return a copy of the provided vector ``v`` with the elements at the indices provided in ``id`` removed. ``inplace_remove()`` ([2]) removes values directly from the provided vector, and is therefore faster.

The first argument ``v`` cannot be a view. The values in ``ids`` are checked to ensure they represent valid indices in ``v``; if not, a run time error is raised.

**Example:**

.. code-block:: c++

    // [1]
    vec1i v = {4,5,2,8,1};
    vec1i w = remove(v, {1,3}); // {4,2,1}

    // [2]
    inplace_remove(v, {1,3});
    // v = {4,2,1}
