.. _Indexing:

Indexing
========

Flat and multidimensional indices
---------------------------------

The ``std::vector``, which is used to implement vif vectors, is a purely linear container: one can access its elements using ``v[i]``, with ``i`` ranging from ``0`` up to ``std::vector::size()-1`` (included). In C++, all indexing is zero-based: the first value has an index of zero (contrary to R, or Fortran).

One can do exactly the same thing on a vif vector. In this context, ``i`` is called a "flat" index, since it conveys no particular structure about the data. The vif vectors go beyond this, and also allow N-dimensional indexing, that is, using the combination of multiple indices to identify one element. This is particularly useful to work on images, which can be seen as 2-dimensional objects where each pixel is identified by its coordinates ``x`` and ``y``. In this context, the pair ``x,y`` is called a "multidimensional" index.

Since regular vectors use the syntax ``v[i]`` to access 1-dimensional data, the natural syntax for 2-dimensional indexing would be ``img[x,y]``. This syntax is valid C++ code, but unfortunately will not do what you expect... This will call the dreaded *comma* operator, which evaluates both elements ``x`` and ``y`` and returns the last one, i.e., ``y``. So this code actually means ``img[y]``. If you ever make this mistake, most compilers should emit a warning, since in this context ``x`` is a useless statement. So watch out! Unfortunately, there is no sane way around it. The only valid syntax is therefore ``img(x,y)``.

Below is an example of manipulation of a 2D matrix:

.. code-block:: c++

    // Create a simple matrix.
    vec2f m = {{1,2,3}, {4,5,6}, {7,8,9}};

    // Index ordering is similar to C arrays: the last index is contiguous
    // in memory. Note that this is *opposite* to the IDL convention.
    m(0,0); // 1
    m(0,1); // 2
    m(1,0); // 4

    // It is also possible to access elements as they are laid out in memory,
    // using square brackets and "flat" indices.
    m[0]; // 1
    m[1]; // 2
    m[3]; // 4

    // ... but doing so with parenthesis will generate a compiler error:
    m(0); // error: wrong number of indices for this vector

The elements of a vector can therefore be accessed in two ways: either using a flat index (and square brackets), or using the appropriate multidimensional index (and parentheses). As illustrated above, when using the multidimensional index, you must provide as many indices as the number of dimensions in the vector, no more, no less.

Indexing a vector can only be done with integers. Indexing with *unsigned* integers is faster, because it removes the need to check if the index is negative, and it should therefore be preferred whenever possible. Negative indices are allowed, and they are interpreted as *reverse* indices, that is, ``-1`` refers to the last element of the vector, ``-2`` to the one before the last, etc.

.. _Safe indexing:

Bounds checking, and safe indexing
----------------------------------

All the indexing methods described above perform *bound checks* before accessing each element. In other words, the vector class makes sure that each index is smaller than either the total size of the vector (for flat indices) or the length of its corresponding dimension (for multidimensional indices). If this condition is not satisfied, an assertion is raised explaining the problem, with a backtrace, and the program is stopped immediately to prevent memory corruption.

.. code-block:: c++

    vec1f v(10):
    v[20] = 3.1415;

When executed, the code above produces:

.. code-block:: none

    error: operator[]: index out of bounds (20 vs. 10)

    backtrace:
     - vif_main(int, char**)
       at /home/cschreib/test.cpp:5


This bound checking has a small but sometimes noticeable impact on performances. In most cases, the added security is definitely worth it. Indeed, accessing a vector with an out-of-bounds index has very unpredictable impacts on the behavior of the program: sometimes it will crash, but most of the time it will not and memory will be silently corrupted. These problems are hard to notice, and can have terrible consequences... Identifying the root of the problem and fixing it may prove even more challenging. This is why these checks are enabled by default, even in "Release" mode.

However, there are cases where bound checking is superfluous, for example if we already know *by construction* that our indices will always be valid, and no check is required. Sometimes the compiler may realize that and optimize the checks away, but one should not rely on it. If these situations are computation-limited, i.e., a lot of time is spent doing some number crushing for each element, then the performance hit of bound checking will be negligible, and one should not worry about it. On the other hand, if very little work is done per element and most of the time is spent iterating from one index to the next and loading the value in the CPU cache, then bounds checking can take a significant amount of the total time.

For this reason, the vif vector also offers an alternative indexing interface, the so-called "safe" interface. It behaves exactly like the standard indexing interface described above, except that it does not perform bound checking. One can use this interface using ``v.safe[i]`` instead of ``v[i]`` for flat indexing, or ``v.safe(x,y)`` instead of ``v(x,y)`` for multidimensional indexing. The safe interface can also be used to create views. This interface is not meant to be used in daily coding, but rather for computationally intensive functions that you write once but use many times. Just be certain to use this interface when you know *for sure* what you are doing, and you can prove that the index will always be valid. A good strategy is to always use the standard interface and, only when the program runs successfully, switch to the safe interface for the most stable but time-consuming code sections.

.. note:: One may wonder why the word ``safe`` was used instead of ``unsafe``, since indexing without bounds checking is actually an "unsafe" operation. The reason why is that writing ``v.safe[i]`` can be understood as: "we are in a context where we know where the index ``i`` came from, we're *safe*, we can disable bounds checking". Perhaps another reason is that I would feel somewhat uncomfortable at writing ``unsafe`` everywhere in the core functions of the library, which is supposed to only contain safe code...


Loops, and traversing data
--------------------------

The fastest way to loop over all the elements of a vector is to use a range-based loop, since this avoids having to deal with indexing and bound checking altogether:

.. code-block:: c++

    for (float& v : m) {
        do_some_stuff_with(v);
    }

If you care not only about the values but also about their flat index in ``m``, then the fastest loop will be:

.. code-block:: c++

    for (uint_t i : range(m)) {
        do_some_more_stuff_with(i, m[i]);
    }

The ``range(m)`` function can only be used in ``for`` loops. It will generate integers from ``0`` to ``m.size()`` (excluded) if ``m`` is a vector, or from ``0`` to ``m`` (excluded) if ``m`` is an integer. It can also have a different starting value when called as ``range(i0,n)``, in which case the first value will be ``i0``.

Lastly, if you care about the multidimensional index, then you need to loop on each dimensions explicitly. When doing so, always loop on the *first* dimension in the outermost loop, and on the *last* dimension in the innermost loop, to make best use of memory locality and CPU caches:

.. code-block:: c++

    for (uint_t i : range(m.dims[0]))
    for (uint_t j : range(m.dims[1])) {
        do_even_more_stuff_with(i, j, m(i,j));
    }


Partial indexing
----------------

When dealing with multi-dimensional vectors, in some cases one may not want to access a single element, but instead all the elements along a certain dimension, or a handful of elements at once. For example, you may want to access an entire row of pixels in an image, or only the values which are greater than zero. This can be done with :ref:`Views`.
