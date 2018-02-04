Finding elements
================

Defined in header ``<phypp/utility/generic.hpp>``.

where
-----

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u where(const vec<Dim,Type>& v);

This function will scan the ``bool`` vector (or view) provided in argument, will store the "flat" indices of all the elements which are ``true``, and will return all these indices in a vector. This is a very useful tool to filter and selectively modify vectors, and probably one of the most used function of the whole library.

**Example:**

.. code-block:: c++

    vec1i v = {4,8,6,7,5,2,3,9,0};

    // We want to select all the elements which are greater than 3.
    // We use where() to get their indices:
    vec1u id = where(v > 3); // {0,1,2,3,4,7}

    // Now we can check:
    v[id]; // {4,8,6,7,5,9}, good!

    // The argument of where() can be quite complex:
    id = where(v < 3 || (v > 3 && v % 6 < 2)); // now guess

    // It can also involve multiple vectors, as long as they have
    // the same dimensions.
    vec1i w = {9,8,6,1,-2,0,8,5,1};
    id = where(v > w || (v + w) % 5 == 0);

    // The returned indices are then valid for both v and w.
    v[id]; // {8,6,7,5,2,9}
    w[id]; // {8,6,1,-2,0,5}

Note that, when called on a view, ``where()`` will return indices inside the view itself, and *not* inside the viewed vector:

**Example:**

.. code-block:: c++

    vec1i v = {4,8,6,7,5,2,3,9,0};

    // Select all the elements greater than 3.
    vec1u id1 = where(v > 3);

    for (uint_t i : range(1)) {
        // We want to apply another selection on top of the first one.
        // Say we now want only those elements which are even (when i=0)
        // or odd (when i=1):
        vec1u id2 = where(v[id1] % 2 == i);

        // Here 'id2' points inside 'v[id1]', not 'v'!
        // To access the correspond values in 'v', one must write:
        v[id1[id2]] = /* ... */;

        // This is dangerous, because 'v[id2]' is perfectly valid,
        // yet makes absolutely no sense.
    }

In general, such situations can be avoided by only calling ``where()`` at the last possible moment, as shown in the example below.

**Example:**

.. code-block:: c++

    vec1b base = v > 3; // only a 'bool' vector for now

    for (uint_t i : range(1)) {
        vec1u id = where(base && v % 2 == i);

        // Now we can access 'v' directly:
        v[id] = /* ... */;
    }

This is often more readable and less error prone, however it may also be less efficient, particularly if the first ``where()`` reduced significantly the number of elements to work with (e.g., if ``v`` contained millions of elements, and only a few are greater than ``3``). A typical case where these situations arise is when binning multidimensional data, for example to build a 2D histogram. Then, a much more efficient approach is to use the ``histogram()`` function and its siblings. Despite the name, these powerful function can be used for purposes other than histograms.


where_first, where_last
-----------------------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    uint_t where_first(const vec<Dim,Type>& v); // [1]

    template<std::size_t Dim, typename Type>
    uint_t where_last(const vec<Dim,Type>& v); // [2]

These functions will scan the ``bool`` vector (or view) provided in argument, and return the "flat" index of the first ([1]) or last ([2]) element which is ``true``, or ``npos`` if all are ``false``.

When used with views, the same caution applies as for ``where()``: the returned index points inside the view itself, not inside the viewed vector.

**Example:**

.. code-block:: c++

    vec1i v = {4,8,6,7,5,2,3,9,0};
    // We want to select the first element which is greater than 3
    uint_t id;
    id = where_first(v > 3); // 0
    v[id];                   // 4
    id = where_last(v > 3);  // 7
    v[id];                   // 9


complement
----------

.. code-block:: c++

    template<std::size_t Dim, typename Type>
    vec1u complement(const vec<Dim,Type>& v, const vec1u& ids);

This function works in parallel with ``where()``. Given a vector ``v`` and a set of "flat" indices ``id``, it will return the complementary set of indices inside this vector, i.e., all the indices of ``v`` that are *not* present in ``id``. The values of ``v`` are unused, only its dimensions are read.

**Example:**

.. code-block:: c++

    vec1i v = {1,5,6,3,7};
    vec1u id = where(v > 4); // {1,2,4}
    vec1u cid = complement(v, id); // {0,3}


match
-----

.. code-block:: c++

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    void match(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2, vec1u& id1, vec1u& id2);

This function returns the indices of the elements with equal values in ``v1`` and ``v2``. In practice, it traverses ``v1`` and, for each value in ``v1``, looks for elements in ``v2`` that have the same value. If one is found, the index of the element of ``v1`` is added to ``id1``, and the index of the element of ``v2`` is added to ``id2``. If other matches are found in ``v2`` for this same value, they are ignored, therefore only the *first* match is returned. Then the function goes on to the next value in ``v1``. The two vectors ``v1`` and ``v2`` need not be the same size.

When used with views, the same caution applies as for ``where()``: the returned indices point inside the views themselves, not inside the viewed vectors.

**Example:**

.. code-block:: c++

    vec1i v = {7,6,2,1,6};
    vec1i w = {2,6,5,3};
    vec1u id1, id2;
    match(v, w, id1, id2);
    id1; // {1,2,4}
    id2; // {1,0,1}
    v[id1] == w[id2]; // always true


set_intersection, set_intersection_sorted, set_union, set_union_sorted
----------------------------------------------------------------------

.. code-block:: c++

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> set_intersection(vec<D1,Type1> v1, vec<D2,Type2> v2); // [1]

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> set_intersection_sorted(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2); // [2]

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> set_union(vec<D1,Type1> v1, vec<D2,Type2> v2); // [3]

    template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
    vec<1,/*...*/> set_union_sorted(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2); // [4]

These functions return a 1D vector containing all the values that exist in both ([1], [2]) or either of ([3], [4]) ``v1`` and ``v2``. These values are returned sorted. If a value is found ``n1`` times in ``v1`` and ``n2`` times in ``v2``, [1] and [2] will return this value ``min(n1,n2)`` times, while [3] and [4] will return this value ``max(n1,n2)`` times. The two vectors do not need to have the same size, and the functions are symmetric: returned values are the same if ``v1`` and ``v2`` are swapped.

The algorithms used internally (``std::set_intersection()`` and ``std::set_union()``) operate on sorted vectors. [1] and [3] will automatically sort the input vectors, so there are no pre-requirement on their ordering. [2] and [4] will assume that the input vectors are already sorted, all will thus be faster.

The type of the returned vector is the common type between ``v1`` and ``v2``, namely, whatever type results of an operation like ``v1[0]+v2[0]``.

**Example:**

.. code-block:: c++

    vec1i v = {1,2,3,3,3,4,5};
    vec1i w = {2,3,3,4,6};

    set_insersection(v, w); // {2,3,3,4}
    set_union(v, w);        // {1,2,3,3,3,4,5,6}


unique_ids, unique_ids_sorted, unique_values, unique_values_sorted
------------------------------------------------------------------

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

These functions will traverse the provided vector ``v`` and find all the unique values. Functions [1] to [3] store the *indices* of these values and return them inside an index vector. If a value is present more than once, the index of the first one will be returned. Functions [4] to [6] directly return the values themselves, rather than indices.

Internally, the algorithm needs to sort ``v``. To optimize execution, several versions of these functions are provided which handle the sorting differently. Functions [1] and [4] will automatically sort the vector, so there is no pre-requirement on ``v``. Functions [2] and [5] takes a second argument ``id`` that contains indices that will sort ``v``. In particular, ``id`` can be the return value of ``sort(v)``. Lastly, functions [3] and [6] assume that ``v`` is already sorted, and are thus the fastest of the three.

When used with views, the same caution applies for functions [1] to [3] as for  ``where()``: the returned indices point inside the views themselves, not inside the viewed vectors.

**Example:**

.. code-block:: c++

    // For an non-sorted vector [1]
    vec1i w = {5,6,7,8,6,5,4,1,2,5};
    vec1u u = unique_ids(w); // {7,8,6,0,1,2,3}
    w[u]; // {1,2,4,5,6,7,8} only unique values

    // Providing a sorting vector [2]
    vec1u s = sort(w);
    vec1u u = unique_ids(w, s); // {7,8,6,0,1,2,3}
    w[u]; // {1,2,4,5,6,7,8} only unique values

    // For a sorted vector [3]
    vec1i v = {1,1,2,5,5,6,9,9,10};
    vec1u u = unique_ids_sorted(v); // {0,2,3,5,6,8}
    v[u]; // {1,2,5,6,9,10} only unique values


is_any_of
---------

.. code-block:: c++

    template<typename Type1, std::size_t Dim2, typename Type2>
    bool is_any_of(const Type1& v1, const vec<Dim2,Type2>& v2); // [1]

    template<std::size_t Dim1, typename Type1, std::size_t Dim2, typename Type2>
    vec<Dim1,bool> is_any_of(const vec<Dim1,Type1>& v1, const vec<Dim2,Type2>& v2); // [2]

Function [1] looks inside ``v2`` if there is any value that is equal to ``v1``. If so, it returns ``true``, else it returns ``false``. Function [2] is the vectorized version of [1], and executes this search for each value of ``v1``, then returns a ``bool`` vector containing the results.

There are no pre-requirements on ``v1`` or ``v2``.

**Example:**

.. code-block:: c++

    vec1i v = {7,4,2,1,6};
    vec1i d = {5,6,7};
    vec1b b = is_any_of(v, d); // {true, false, false, false, true}


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

These functions use a binary search algorithm to locate the element in the input vector ``v`` that is equal to or closest to the provided value ``x``, which must be a scalar. The binary search assumes that the elements in the input vector are *sorted* by increasing value. This algorithm also assumes that, if the input vector contains floating point numbers, none of them is ``NaN``.

``lower_bound()`` ([1]) locates the last element in ``v`` that is *less or equal* to ``x``. If no such element is found, ``npos`` is returned.

``upper_bound()`` ([2]) locates the first element in ``v`` that is *greater than* ``x``. If no such element is found, ``npos`` is returned.

``bounds()`` ([3]) combines what ``lower_bound()`` and ``upper_bound()`` do, and returns both indices in an array. The second overload of ``bounds()`` ([4]) calls ``lower_bound()`` to look for ``x1``, and ``upper_bound()`` to look for ``x2``.

**Example:**

.. code-block:: c++

    vec1i v = {2,5,9,12,50};
    bounds(v, 0);   // {npos,0}
    bounds(v, 9);   // {2,3}
    bounds(v, 100); // {4,npos}


equal_range
-----------

.. code-block:: c++

    template<typename T, std::size_t Dim, typename Type>
    vec1u equal_range(const vec<Dim,Type>& v, T x);

This function uses a binary search algorithm to locate all the values in the input vector ``v`` that are equal to ``x``. The binary search assumes that the elements in the input vector are *sorted* by increasing value. This algorithm also assumes that, if the input vector contains floating point numbers, none of them is ``NaN``.

The function returns the indices of all the values equal to ``x``. If no such value is found, an empty vector is returned.

If ``v`` is not sorted, the only alternative is to call ``where(v == x)``; this will be slower than ``equal_range()``, but it should still be faster than sorting ``v``.

**Example:**

.. code-block:: c++

    vec1i v = {2,2,5,9,9,9,12,50};
    equal_range(v, 9); // {3,4,5}

    // The above is a faster version of:
    where(v == 9);


astar_find
----------

.. code-block:: c++

    bool astar_find(const vec2b& map, uint_t& x, uint_t& y);

This function uses the A* ("A star") algorithm to look inside a 2D boolean map ``m`` and, starting from the position ``x`` and ``y`` (i.e. ``m(x,y)``), find the closest point that has a value of ``true``. Once this position is found, its indices are stored inside ``x`` and ``y``, and the function returns ``true``. If no element inside ``m`` is ``true``, then the function returns ``false``.

**Example:**

.. code-block:: c++

    // Using 'o' for true and '.' for false, assume we have the following boolean map,
    // and that we start at the position indicated by 'S', the closest point whose coordinates
    // will be returned by astar_find() is indicated by an 'X'

    vec2b m;

    //   0123456789  13
    // 0 .................
    // 1 .................
    // 2 .................
    // 3 ...ooooo.........
    // 4 ...ooooo.........
    // 5 ...ooooo.........
    // 6 ...ooooo.........
    // 7 ...ooooX.........
    // 8 .............S...
    // 9 .................
    //   .................

    uint_t x = 13, y = 8;
    astar_find(m, x, y);
    x; // 7
    y; // 7
