.. _Views:

Views
=====

Overview
--------

As shown in :ref:`Operator overloading`, instead of accessing each element of a vector individually to perform some operation, we can use operator overloading to act on all the elements of the vector at once: ``v *= 2``. As shown in :ref:`Indexing`, we can also, like with ``std::vector``, modify each element individually using their indices: ``v[2] *= 2``.

The *view* is a generalization of this concept, allowing you to access and operate on an arbitrary number of elements from an existing vector. Each element of the view is actually a *reference* to an element in this vector, and performing operations on the elements of the view will modify the elements in the vector. The interface of views is almost indistinguishable from that of vectors, and both can be used interchangeably in almost all cases. Generic codes and functions that work with vectors will therefore also work with views.

Creating a view is simple: instead of indexing a vector with integer values, as in ``v[2]``, you would index the vector using *another vector* containing multiple indices:

.. code-block:: c++

    // Create a simple vector.
    vec1f w = {1,2,3,4,5,6};

    // We want to "view" only the second, the third and the fifth
    // elements. So we first create a vector containing the
    // corresponding indices.
    vec1u id = {1,2,4};

    // We create the view and multiply all the elements by two.
    w[id] *= 2;
    // 'w' now contains {1,4,6,4,10,6};

In the example above, the type of the epxression ``w[id]`` is ``vec<1,float*>``. Views have the same elements type as the vector they point to, but the template parameter of ``vec<...>`` is specified as a *pointer* (``T`` -> ``T*``) to distinguish them from vectors. The dimensions of the view are set by the vector that was used for *indexing* (regardless of the dimensions of the pointed vector). For example, if ``v`` is a 1D vector and we create a 2D array of indices ``id``, then ``v[id]`` will be a 2D view.

.. note:: Since a view keeps *references* to the elements of the original vector, the lifetime of the view must not exceed that of the original vector. Else, it will contain *dangling* references, pointers to unused memory, and this should be avoided at *all cost*. For this reason, views are not meant to be stored into named variables, but should only be used in temporary expressions as above. This point is discussed also in :ref:`Known issues and problems`. The only notable exception to this rule is when passing views to functions (see :ref:`Generic function guidelines`).

As for regular indexing, views can only be created using vectors of *integer* indices (signed, or unsigned). Bounds checking will be done on each element of the index vector, to make sure that no index goes past the end of the vector. If you know this cannot happen, and therefore that this bounds checking is superfluous, you may want to use "safe" indexing (see :ref:`Indexing`) to improve performances.

Lastly, views on different vectors can be involved in the same expression, as long as their dimensions are the same:

.. code-block:: c++

    // Create two vectors.
    vec1f x = {1,2,3,4,5,6};
    vec1f y = {6,5,4,3,2,1};

    // Indices
    vec1u idx = {1,2,4};
    vec1u idy = {4,0,5};

    // Do some operation
    vec1f z = x[idx] + y[idy]; // {2,3,5} + {2,6,1} = {4,9,6}

Mixing views (or vectors) of different sizes will trigger an error at runtime.


Range indexing
--------------

Sometimes, one will want to use views to access all the elements at once, for example to set all the elements of a vector to a specific value. This can be done with a loop, of course, but the whole point of phy++ is to avoid explicit loops whenever possible. An alternative is to use a view, with an index vector that contains all the indices of the target vector:

.. code-block:: c++

    vec1i v = {1,2,3,4};
    vec1u id = {0,1,2,3}; // all the indices of 'v'
    v[id] = 12;           // all the values are now equal to 12

    // Note that, by design, the following will not compile (too error prone):
    v = 12; // "error: no viable overloaded '='"

However, not only is this not very practical to write, it is error prone and not very clear. If someday we decide to add an element to ``v``, we also have to modify ``id``. Not only this, but it will most likely be slower than writing the loop directly, because the compiler may not realize that you are accessing all the elements contiguously, and will fail to optimize it properly.

The optimal way to do this in phy++ is to use the "placeholder" symbol, defined as a single underscore ``_``. When used as an index, it means "all the indices in the range". Coming back to our example:

.. code-block:: c++

    vec1i v = {1,2,3,4};
    v[_] = 12; // it cannot get much shorter!

This placeholder index can be used in all situations, with both flat and multidimensional indexing:

.. code-block:: c++

    vec2f img(128,128);
    img(0,_) = 12; // accessing the first row of the image

    // Any combination is allowed
    vec4f crazy(5,4,12,8);
    crazy(5,_,2,_) = 5.0; // this creates a 2D view of shape 4x8

    // The above is equivalent to:
    for (uint_t i : range(crazy.dims[1]))
    for (uint_t j : range(crazy.dims[3])) {
        crazy(5,i,2,j) = 5.0;
    }

This can be further refined to only encompass a fraction of the whole range, using a specific syntax:

.. code-block:: c++

    vec1i v = {1,2,3,4};
    v[_-2] = 12;   // only access the indices from 0 to 2 (included)
    v[2-_] = 12;   // only access the indices from 2 to 3 (the last, included)
    v[1-_-2] = 12; // only access the indices from 1 to 2 (included)

    // Watch out, this is *not* range indexing!
    v[1-2] = 12;   // only access index 1-2 = -1


Filtering and selecting elements
--------------------------------

In the previous section we have seen that a view can be created using a vector of indices. In most cases, such vector is not created manually, as in the examples above, but comes from a *filtering* function, ``where()``. This function is part of the support library, but it is important enough to be mentioned here.

``where()`` accepts a vector of ``bool`` (of any dimension) as single argument, and returns all the *flat* indices where the vector values are ``true``. This can be combined with views to perform complex operations on vectors. For example:

.. code-block:: c++

    // Set all negative values to zero
    vec1f v1 = {-1.01, 2.0, 5.0, -2.1, 6.5};
    v1[where(v1 < 0.0)] = 0.0;
    v1;     // { 0.0,  2.0, 5.0,  0.0, 6.5}

    // Add one to all values between 0 and 6
    vec2f v2 = {{-1.0, 2.0}, {8.0, 3.4}};
    v2[where(v2 > 0.0 && v2 < 6.0)] += 1.0;
    v2;     // {{-1.0, 3.0}, {8.0, 4.4}}


Differences between views and vectors
-------------------------------------

While views are mostly compatible with vectors in terms of interface, by design some features of vectors are not available for views:

* Initialization: views can only be created as described above.
* Assignment and resizing: assigning anything to the view will affect the target vector, not the view itself. Therefore once a view is created, you cannot change which elements it points to.


Constant views and views on constant data
-----------------------------------------

There are two ways that views can have "constant" semantics, where it is only possible to *read* the viewed data and not modify it. The first way is when constructing a view from a constant vector, in which case the view carries the ``const`` qualifier in its data type (``vec<1,const int*>``):

.. code-block:: c++

    const vec1i v = {1,2,3,4};
    v[_] = 12; // error: cannot modify values of vec<1,const int*>

The second way arises when views are function parameters (see :ref:`Generic function guidelines` for more detail):

.. code-block:: c++

    void set_values(const vec<1,int*>& v) {
        v[_] = 12; // error: cannot modify values of const vec<1,int*>
    }

There is no difference between these two cases: "a constant view on non-constant data" and "a view on constant data", ``const vec<1,float*>`` is semantically identical to ``vec<1,const float*>``. This is different from raw pointers, because a pointer can be modified to point to a different value, while views cannot (by design).


Aliasing
--------

The implementation of vectors and views in phy++ is such that aliasing *never* occurs in vectorized operations. More precisely, any assignment of the form ``x = y`` (or ``x += y``, etc.) occurs *as if* executed in the following order:

1. ``y`` (the right-hand-side) is evaluated,
2. the values of ``y`` are copied in a temporary vector,
3. ``x`` (the left-hand-side) is evaluated,
4. the values of the temporary vector are assigned to the elements of ``x``.

In practice, the creation of the temporary vector (step 2) may be dropped for optimization purposes, but only in cases where it would not change the outcome of the operation, that is, when aliasing is guaranteed not to occur. The following illustrates when aliasing *could* occur, and describes in practice how it is avoided in phy++.

Because views hold *references* to existing data, there is the possibility of the same data being read and modified in the same expression. This is, essentially, what is called "aliasing":

.. code-block:: c++

    vec1i v = {1,2,3,4};
    vec1u id = {1,2,3,0};
    v[id] = v; // what happens here?

This can create confusing situations, like the above, where it matters in which *order* the operations are performed. These situations are identified using a check, made prior to every assignment between a vector and view, a view and a vector, or two views. Each view carries a pointer to the original vector: if this pointer matches the vector involved in the assignment (or the pointer of the other view), then aliasing is detected. In such cases, the data on the *right* side of the equal sign is copied to a temporary vector, which is then assigned to the data on the *left* side of the equal sign. In all other cases, aliasing is ignored and no temporary is created to avoid the performance hit.

So, the example above first creates a copy of ``v``, then assigns it to itself following the order in the view. The vector then contains the values ``{4, 1, 2, 3}``, as one would expect if the data on the right side of the equal sign originated from another vector. If aliasing had not been detected, one possible outcome would have been ``{1, 1, 1, 1}``, as some of the vector's values would have been modified *before* being read.

A similar problem can arise without views:

.. code-block:: c++

    vec1i v = {1,2,3,4};
    v += v[0]; // what happens here?

Possible outcomes are ``{2,3,4,5}`` if ``v[0]`` is treated as the *value* ``1``, or ``{2,4,5,6}`` if ``v[0]`` is treated as the *reference* to the first element of ``v``, leading to aliasing. To avoid the latter, assigning operators such as ``+=`` always take scalar arguments by value. The outcome will therefore be ``{2,3,4,5}``.

This means that the above codes are *not* identical to their equivalent with explicit loops:

.. code-block:: c++

    vec1i v = {1,2,3,4};
    vec1u id = {1,2,3,0};

    for (uint_t i : range(v)) {
        v[id[i]] = v[i];
    }

    // v = {1,1,1,1}, aliasing *did* occur

    v = {1,2,3,4};

    for (uint_t i : range(v)) {
        v[i] += v[0];
    }

    // v = {2,4,5,6}, aliasing *did* occur

While this may cause confusion, experience has shown that aliasing is more often an unwanted nuisance than a feature. Furthermore, with the explicit loop it is immediately apparent that ``v[i]``, ``v[id[i]]``, or ``v[0]`` will be re-evaluated on each iteration, therefore that the corresponding value may change.
