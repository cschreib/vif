Vectors
=======

Overview
--------

A vector in vif is basically an enhanced ``std::vector`` (which, in fact, is used to implement the vif vectors), and it therefore shares most of its features and strengths.

In particular, a vector can contain zero, one, or as many elements as your computer can handle. Its size is defined at *runtime*, meaning that its content can vary depending on user input, and that it can change its total number of elements at any time. The elements of a vector are stored *contiguously* in memory, which provides optimal performances in most situations. Lastly, a vector is an *homogeneous* container, meaning that a given vector can only contain a single type of elements; for example ``int`` or ``float``, but not both.

On top of the ``std::vector`` interface, the vif vectors have some extra functionalities. The most important ones are :ref:`Operator overloading` (which allows writing ``v+w`` instead of writing a loop to sum all the elements one by one), structured :ref:`Indexing` for multi-dimensional data (i.e., images, data cubes, and objects of higher dimensions), and :ref:`Views` (which allow accessing and modifying subsets of existing vectors).

Like in most advanced C++ libraries, vif vectors are *template-based*. This means they are designed to work with *any* data type ``T`` for their elements, for example ``int``, ``float``, or ``std::string``. The type of a vector is therefore spelled ``vec<1,T>``, where ``T`` can be replaced by any type (it could be a vector itself). There is no explicit restriction regarding the data type ``T``, however some features may obviously not be available depending on the capabilities of your type. For example, if your type has no ``operator*`` (such as ``std::string``), you will not be able to multiply vectors of this type. Lastly, the vif vector shares the same restrictions as the ``std::vector`` regarding the *copyable* and *movable* capabilities of the stored type.

The number of dimensions of a vector is specified in its type; this is the ``1`` in ``vec<1,T>``. For example, a 2D image of ``float`` will be declared as ``vec<2,float>`` (or ``vec2f``, see :ref:`Type aliases` below), while a 1D tabulated data set of ``int`` will be ``vec<1,int>`` (or ``vec1i``). The fact that the number of dimensions is part of the type means that, while the number of *elements* in a vector is determined at runtime, the multi-dimensional nature of a vector is determined at *compile time*. In other words, a 1D vector cannot be turned into a 2D vector; you would have to create a new variable (using the ``reform()`` and ``flatten()`` functions). However it is possible to change, say, a 128x128 2D vector into a 256x256 2D vector. Only the *number* of dimensions is fixed, the length of each dimension is free to vary.

This restriction was imposed for two reasons: first, type safety, and second, performance. Since the compiler knows how many dimensions there are, it will output an error whenever you try to perform operations on two vectors of different dimensionality (for example, adding a 1D array to a 2D image; which would make no sense). In terms of performance, this also means that we also know at the time of compilation how many dimensions we need to deal with, so the compiler can more easily decide whether (or how) to unroll loops.

The hard limit on the number of dimensions depends on your compiler, as each dimension involves an additional level of template recursion. The C++ standard does not guarantee anything in this respect, but you should be able to go as high as 256 on all major compilers. Beyond this, you should probably see a therapist first.


.. _Type aliases:

Type aliases
------------

While templates are a fantastic tool for library writers, they can easily become a burden for the *user* of the library, because of the additional syntax complexity (the ``<...>`` in the name of the vector type). Since vif is a numerical analysis library, we know in advance what types will most often be stored inside the vectors, and we therefore introduce type aliases for the most common vector types:

* ``vec1f``: vector of ``float``,
* ``vec1d``: vector of ``double``,
* ``vec1cf``: vector of ``std::complex<float>``,
* ``vec1cd``: vector of ``std::complex<double>``,
* ``vec1i``: vector of ``int`` (precisely, ``int_t = std::ptrdiff_t``),
* ``vec1u``: vector of ``unsigned int`` (precisely, ``uint_t = std::size_t``),
* ``vec1b``: vector of ``bool``,
* ``vec1s``: vector of ``std::string``,
* ``vec1c``: vector of ``char``.

Such type aliases are provided for dimensions up to 6 (i.e., ``vec6f`` exists, but not ``vec7f``).


.. _Size and dimensions:

Size and dimensions
-------------------

For any vector ``v``, ``v.dims`` is an ``std::array`` that holds the length of each dimension of the vector: ``v.dims[0]`` is the number of elements along the first dimension, etc. This array is set automatically when the vector is created or resized, and should not be modified manually (except in very specific and rare circumstances).

The total number of elements in a vector can be obtained using ``v.size()``, and this is equal to the product of all the dimensions of the vector. This is a cheap function, as the size is stored internally and does not need to be computed for each call. One can check if a vector is empty (i.e., contains no element) using ``v.empty()``, which returns either ``true`` or ``false`` (this does *not* empty the vector).


.. _Initialization:

Initialization and assignment
-----------------------------

There are several ways to initialize and assign data to a vector:

* The "default" initialization, where the vector is empty.
* The "size" initialization, where the vector contains ``n`` default-constructed elements.
* The "list" initialization, where the vector is assigned a set of values.
* The "copy" initialization, where the vector contains a copy of the data from another vector.
* The "move" initialization, where the vector steals the data from another vector.

The last three can be used both for initialization and assignment (using ``operator=``).

The "default" initialization is very cheap, since it involves no (or very little) allocation:

.. code-block:: c++

    vec1f v; // empty

The "size" initialization pre-allocates memory for the data, which is very useful if you know in advance how many elements your vector needs to contain. The allocated data consists of elements which are default-constructed; this means a value of ``0`` for arithmetic types, ``false`` for ``bool``, and empty strings for ``std::string``.

.. code-block:: c++

    vec1f v(10); // 10 zeros
    vec2f w(10,20); // 200 zeros, arranged in a 2D shape of 10x20
    vec3f z(w.dims,4); // 800 zeros, arranged in a 3D shape of 10x20x4

As shown in the last example (initializing ``z``), the arguments of the constructor can be a mixture of integer values and ``std::array``; the arrays can come, for example, from the dimensions of other existing vectors.

The "list" initialization explicitly specifies a set of values to be stored in the vector. This uses initializer lists, which can be nested for multidimensional vectors.

.. code-block:: c++

    vec1f v = {1, 2, 3, 4); // values from 1 to 4
    vec2f w = {{1, 2}, {3, 4}, {5, 6}}; // values from 1 to 6 arranged in a 2D shape of 3x2

The "copy" initialization trivially copies (and optionally converts, see :ref:`Type conversion`) the values of a vector into another one.

.. code-block:: c++

    vec1f v = {1, 2, 3, 4); // values from 1 to 4
    vec1f w = v;            // also contains values from 1 to 4

The "move" initialization will "steal" the values of another vector. The vector from which the values are "stolen" then becomes empty, and can be reused for other purposes later. This will usually be much faster than the "copy" initialization above, if you do not mind the side effect.

.. code-block:: c++

    vec1f v = {1, 2, 3, 4); // values from 1 to 4
    vec1f w = std::move(v); // also contains values from 1 to 4, but now 'v' is empty


.. _Resizing and adding elements:

Resizing and adding elements
----------------------------

The dimensions and size of a vector can be modified in three main ways.

First, ``v.clear()`` will erase all the values from the vector and set all its dimensions to zero. This will set the vector in the state of a default-initialized vector (see above).

Second, ``v.resize(...)`` will change the dimensions of the vector. The parameters it accepts are the same as the "size" constructor (see above), i.e., either integral values for individual dimensions, or an ``std::array`` containing multiple dimensions, or any combination of these. While the total number of elements can be modified at will, the number of *dimensions* of the vector must remain unchanged.

.. code-block:: c++

    vec2f w;            // empty, zero elements
    w.resize(20,10);    // 20x10, 200 elements
    w.resize(200,10);   // 200x10, 2000 elements
    w.resize(200,10,5); // error: cannot change a 2D vector into a 3D vector

Once the vector has been resized, its previous content is left in an undefined state, i.e., you can generally assume the previous values (if any) have been lost and replaced by meaningless garbage. The only exception is for 1D vectors. If the resize operation *decreases* the total number of elements, then values are erased from the end of the vector and the rest remains untouched. On the other hand, if the resize operation *increased* the total number of elements, new elements are inserted at the end of the vector, default constructed (i.e., zeros for integral types, etc.). This is the same behavior as ``std::vector``.

Third, ``v.push_back(...)`` will add new values at the end of the vector, increasing its size. The behavior of this function is different for 1D and multidimensional vectors. For 1D vectors, this function appends a new element at the end of the vector, and therefore takes for argument a single scalar value. For multidimensional vectors, this function takes for argument another vector of ``D-1`` dimensions, and whose lengths match the *last* ``D-1`` dimensions of the first vector. The new vector is inserted after the existing elements in memory, and the *first* dimension of the first vector is increased by one.

.. code-block:: c++

    vec1i v = {1,2,3};
    v.push_back(4); // {1,2,3,4}

    vec2i w = {{1,2,3}, {4,5,6}}; // shape 2x3
    w.push_back({7,8,9});         // shape 3x3
    w.push_back({7,8});           // error: dimensions do not match


For optimization, the ``push_back(...)`` function will generally be used in conjunction with ``v.reserve()``. This function is identical to ``std::vector::reserve()``. To understand what this function actually does, one needs to know the internal behavior of ``std::vector``. At any instant, the ``std::vector`` only has enough memory to hold ``N`` elements, which is usually larger than the actual size of the vector. ``N`` is called the capacity of the vector. Once the allocated memory is full, and a new ``push_back()`` is called, ``std::vector`` allocates a larger amount of memory (typically ``2*N`` elements), copies the existing elements inside this new memory, and frees the old memory. This strategy allows virtually unlimited growth of a given vector; it is quite efficiently tuned, but it remains an expensive operation. Performances can be greatly improved if one knows *in advance* (and even approximatively) the total number of objects that need to be stored in the vector, so that the right amount of memory is allocated from the start, and no further reallocation is required. This function does just that: it tells ``std::vector`` how many elements it will (or might) contain at some point, so that the vector can already allocate enough memory to store them contiguously. Later, if you have reserved way too much memory, you can always ask the vector to free the surplus by calling ``std::vector::shrink_to_fit()``, which will result in an additional reallocation but will free some unused memory.


.. _Type conversion:

Type conversion, and casting
----------------------------

The rules for converting a vector of a type ``T`` into a vector of another type ``U`` follow the rules for converting ``T`` itself into ``U``. If ``T`` is implicitly/explicitly convertible to ``U``, then it is always possible to implicitly/explicitly convert a ``vec<1,T>`` into ``vec<1,U>``. For example here with a conversion from ``vec1f`` to ``vec1i``:

.. code-block:: c++

    vec1f v1 = {1.5, -2.2, 100.0};
    vec1i v2 = v1; // this works

There is one notable exception to this rule, which is for vectors of type ``bool``. In C++, ``bool`` can be implicitly converted to (and from) any other arithmetic type (such as ``int`` or ``float``). While implicit conversion is very convenient in most cases, in the case of ``bool`` the risk of unwanted narrowing conversion (where data is lost) is much greater, while the actual use cases for implicit conversion are rarer; ``bool`` indeed carries a very different semantic compared to the other arithmetic types. For this reason, in vif it was decided to disable implicit conversion to and from ``bool``. If needed, the conversion is still possible at no extra cost by using an explicit cast:

.. code-block:: c++

    vec1f v1 = {1.5, -2.2, 100.0};
    vec1b v2 = v1;        // does *not* work! compiler error
    vec1b v2 = vec1b{v1}; // this works


.. _Operator overloading:

Operator overloading
--------------------

When dealing with ``std::vector``, the only thing you can do to operate on all the elements of an ``std::vector`` is to iterate over these elements explicitly, either using a C++11 range-based loop, or using indices:

.. code-block:: c++

    // Goal: multiply all elements by two.
    std::vector<float> v = {1,2,3,4};

    // Either using a range-based loop,
    for (float& x : v) {
        x *= 2;
    }

    // ... or an index-based loop.
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] *= 2;
    }

While this is fairly readable (especially the first version), it is still not very concise and expressive. For vif vectors, we have *overloaded* the usual mathematical operators to make it possible to write the above code in a much simpler way:

.. code-block:: c++

    // Using vif vector.
    vec1f v = {1,2,3,4};
    v *= 2;

Not only this, but we can also perform operations on a pair of vectors:

.. code-block:: c++

    // Goal: sum the content of the two vectors.
    vec1f x = {1,2,3,4}, y = {4,3,2,1};
    vec1f z = x + y;
    // z: {5,5,5,5}

Almost all the mathematical and logical operators are overloaded. Therefore, as a rule of thumb, if you can do an operation with a type ``T``, you can do it with ``vec<1,T>`` as well. The one notable exception are bitwise operators: ``|``, ``&``, and ``^``. The reason is twofold: first, these are not so commonly used in data analysis, and second, the ``^`` operator can be mistakenly interpreted as the exponentiation operator, that some other languages possess (if you need to do exponentiation, use ``pow()``).

.. note:: In Python+NumPy and some C++ libraries (like xtensor_), it is possible to mix vectors with different numbers of dimensions in a given operation, if some conditions are satisfied. For example, this allows multiplying the same 1D vector to each row or column of a 2D vector without having to write an explicit loop. This mechanism is called *broadcasting*. It is *not* implemented in vif, and will likely never be. Indeed, it is perceived that the benefits are not worth the costs, both in terms of making the vif codebase more complex, but also in introducing more complex rules for the user of the library. Therefore, in vif, one can only do arithmetics on vectors which have the *exact* same dimensions. If you require an operation similar to what broadcasting provides, you can always write the loop explicitly.

.. note:: Contrary to some other C++ libraries with vectorized arithmetic (such as Eigen_, blazelib_, or xtensor_), vif does not use *expression templates*. Instead, each operation is executed immediately (no lazy evaluation) and operates if necessary on temporary intermediate vectors. While this may appear to be a sub-optimal implementation, vif was tuned to makes good use of return value optimization, move semantics, and for reusing the memory of temporaries in chained expressions. As a result, performance was found to be on par with expression templates in the most common situations, but memory consumption is generally higher in vif. This is of course dependent on the precise calculation to perform. The benefit of not using expression templates is a reduced compilation time, and a much simpler code base.

.. _Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
.. _blazelib: https://bitbucket.org/blaze-lib/blaze
.. _xtensor: https://xtensor.readthedocs.io/en/latest/

Constant vectors
----------------

It is not possible to create a vector of constant elements, ``vec<1, const int>`` will not compile. The correct way to create a constant vector is to make the vector itself constant: ``const vec<1,int>``.

