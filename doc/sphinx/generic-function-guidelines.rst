.. _Generic function guidelines:

Guidelines for writing generic functions
========================================

This section contains more advanced technical details about the implementation of vectors and views in phy++. It also includes tips and tricks for writing correct, efficient, and generic functions.


What is a generic function?
---------------------------

A "generic" function can operate on vectors regardless of the precise type of their elements. For example, a function to shuffle the values inside a vector does not care whether these values are integers, strings, or potatoes, it just needs to know how many values there are.

In C++, such generic functions are written using *template metaprogramming*:

.. code-block:: c++

    template<typename T>
    void multiply_second_row(T& v) {
        v(0,_) *= 2;
    }

    vec2f v1 = {{1.0}, {5.0}, {-1.0}}; // float, shape 3x1
    multiply_second_row(v1);

    vec2i v2 = {{1, 5}, {-1, 0}}; // int, shape 2x2
    multiply_second_row(v2);


Expressing constraints on function arguments
--------------------------------------------

In the example above, the type of the function's argument ``v`` is ``T``, and is totally unconstrained. It could be anything. However this specific function has an *implicit* constraint on the type ``T``: it must be possible to write ``v(0,_) *= 2``. This means ``v`` *must* be a 2D vector, and the elements must be of arithmetic type. If you try to use this function on a value for which does not satisfy this implicit constraint, the compiler will generate an error:

.. code-block:: c++

    vec1f v4 = {1.0, 5.0, 6.0};
    multiply_second_row(v4); // error!

The program will not compile, which is good. However the error message will be *nasty* (180 lines of errors with Clang), and the error will point to code *inside* the function. This is far from ideal, because the user of the function should not need to understand what the function does internally to fix the error. This can be fixed by adding constraints on the type ``T``:

.. code-block:: c++

    template<typename T>
    void multiply_second_row(vec<2,T>& v) {
        v(1,_) *= 2;
    }

Here, we state that ``v`` must be a 2D vector, and we leave the type of the elements unconstrained. Using this new definition of ``multiply_second_row()``, the error message in the case above becomes much smaller (4 lines of errors), and explicitly says that there is no matching function call for ``multiply_second_row(v4)``. This is much better. Therefore you should always make sure to specify as many constraints as possible in the *signature* of the function (i.e., the type of its arguments).

We are not done though. Indeed, we still left the type of the elements unconstrained, while we need elements of arithmetic types to be able to write ``v(1,_) *= 2``. For example, using a vector of strings would be an error:

.. code-block:: c++

    vec2s v5 = {{"I", "am"}, {"a", "string"}};
    multiply_second_row(v5); // error!

Again, this emits a length error from inside the function (20 lines of errors). We can fix this by adding extra constraints on the type ``T`` of the elements. One possibility is to force it to be some "common" type, like ``double``:

.. code-block:: c++

    // Note: parameter is fully constrained, it not generic anymore
    void multiply_second_row(vec<2,double>& v) {
        v(1,_) *= 2;
    }

This makes the error much easier to understand (7 lines of errors), but it has the important downside that the function is no longer generic: it *needs* a vector of ``double``. If you try to call it on a vector of ``float``, it will not compile. It will also fail to work on *views* (see below). So unless you know the function should only be used with ``double`` values, this is not the right solution. Instead, we can leave the type of elements to be generic, and use ``std::enable_if<>`` to express a constraint on this type, in this case ``std::is_arithmetic<T>``:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void multiply_second_row(vec<2,T>& v) {
        v(1,_) *= 2;
    }

With this version of the function, the error when called on vectors of strings becomes much clearer (4 lines of errors) and says that you cannot call the function on strings. Again, much better!

So, that's it? Not quite. There is one last implicit requirement when we write ``v(1,_)``: the first dimension of ``v`` must have at least two elements. There is no way to check this at the time of compilation, so the program will compile:

.. code-block:: c++

    vec2i v6;
    multiply_second_row(v6); // compiles, but runtime error!

It will fail at runtime though. The backtrace will show that the error happened in ``multiply_second_row()``, but with a rather cryptic error message:

.. code-block:: c++

    error: operator(): index out of bounds (1 vs. 0)

The solution here is to perform an explicit check inside the function, and emit a clearer error message using the ``phypp_check()`` function:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void multiply_second_row(vec<2,T>& v) {
        phypp_check(v.dims[0] >= 2, "vector must have at least two elements along first dimension ",
            "(got ", v.dims[0], ")");
        v(1,_) *= 2;
    }

The error shown to the user then becomes clear:

.. code-block:: c++

    error: vector must have at least two elements along first dimension (got 0)

**Note:** Since we now do an explicit check that the index ``1`` is valid before accessing the vector, we no longer need the vector to perform bounds checking. Therefore we can use the "safe" indexing interface:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void multiply_second_row(vec<2,T>& v) {
        phypp_check(v.dims[0] >= 2, "vector must have at least two elements along first dimension ",
            "(got ", v.dims[0], ")");
        v.safe(1,_) *= 2;
    }

This is the optimal way to write this function, and it is clearly not as pretty as the very first version. This shows that, while writing generic functions is easy, writing them *well* is much harder. For this reason, always check in the support library if a function already exists before writing your own.

It should be said, however, that the very first version we wrote actually does the work we expect it to do. It is not "incorrect"; its only defect is that it will not be very helpful when things go wrong.


Supporting both vectors and views
---------------------------------

For a vector of type ``vec<D1,T>``, a view will have a type ``vec<D2,T*>`` or ``vec<D2,const T*>``. The number of dimensions can be different, and the data type is a pointer to the type of elements in the vector. The ``const`` qualifier is used to propagate const-correctness if the original vector was declared ``const``.

This makes it relatively easy to write function that work on both vectors and view, but this distinction means that there are a number a details to keep in mind. Consider this generic function that computes the sum of all the elements in a vector:

.. code-block:: c++

    template<std::size_t D, typename T>
    T total(const vec<D,T>& v) {
        T ret = 0;
        for (const T& val : v) {
            ret += val;
        }

        return ret;
    }

This implementation works for all vectors, but it will fail for views. Indeed, if called on a view of type ``vec<1,int*>``, then ``T = int*``, and the return value is not an integer but an (invalid!) pointer to an integer. Fortunately, it will not even compile because the loop will try to assign the values of ``v`` to a ``const int*&``, which will fail. Therefore, the type ``T`` should never be used directly like this.

Instead, you should apply the transform ``meta::rtype_t<T>``, which essentially transforms ``T*`` into ``T`` and removes const qualifiers, and use ``auto`` whenever possible to let the type system make the right decisions for you:

.. code-block:: c++

    template<std::size_t D, typename T>
    meta::rtype_t<T> total(const vec<D,T>& v) {
        meta::rtype_t<T> ret = 0;
        for (const auto& val : v) {
            ret += val;
        }

        return ret;
    }


Vectorizing scalar functions
----------------------------

Most function created in C++ thus far, including those in the C++ standard library, are *scalar* functions which operate on one single value. The best example of this are all the mathematical functions, ``sqrt()``, ``pow()``, ``ceil()``, etc. These functions can be *vectorized* to operate directly on vector data without having to write a loop. The phy++ support library contains a large number of such vectorized functions:

.. code-block:: c++

    double v1 = 2.0;
    sqrt(v1); // 1.41...
    vec1d v2 = {2.0, 4.0, 6.0};
    sqrt(v2); // {1.41..., 2.0, 2.45...}

However the phy++ support library cannot contain *all* functions that ever existed, and you may create your own scalar functions that you wish to vectorize. This can be achieved using the preprocessor macro ``PHYPP_VECTORIZE()```:

.. code-block:: c++

    float myfunc(float v) {
        return sqrt(3*v + 5.0); // whatever you wish to do
    }

    PHYPP_VECTORIZE(myfunc)

The macro must be called in the global scope, inside a namespace, or a the root scope of a class. It *cannot* be called inside a function. This macro emits two additional functions with the same name. The first function is the most generic vectorized version of the scalar version, which will get used most of the time.

The second version offers an interesting optimization opportunity when the return type is the same as the argument type (as is the case for ``myfunc``), and when the function is called on a temporary vector (not views). This optimized version reuses the memory of the temporary vector instead of returning a brand new vector. This offers important optimizations in case of chained calls:

.. code-block:: c++

    vec1d v1 = {1.0, 1.2, 1.5};
    vec1d v2 = myfunc(sqrt(v1));

In this example, ``sqrt(v1)`` creates a temporary vector, and ``myfunc()`` applies ``myfunc()`` in-place on the values of the temporary vector. It is equivalent to this:

.. code-block:: c++

    vec1d v1 = {1.0, 1.2, 1.5};

    vec1d tmp(v1.dims);
    for (uint_t i : range(v1)) {
        tmp[i] = sqrt(v1[i]);
    }
    for (double& v : tmp) {
        v = myfunc(v);
    }

    vec1d v2 = std::move(tmp);

Without this optimization, the chained call would have created two temporaries:

.. code-block:: c++

    vec1d v1 = {1.0, 1.2, 1.5};

    vec1d tmp1(v1.dims);
    for (uint_t i : range(v1)) {
        tmp1[i] = sqrt(v1[i]);
    }
    vec1d tmp2(tmp1.dims);
    for (uint_t i : range(tmp1)) {
        tmp2[i] = myfunc(tmp1[i]);
    }

    vec1d v2 = std::move(tmp2);

The optimal version would avoid the extra loop:

.. code-block:: c++

    vec1d v1 = {1.0, 1.2, 1.5};

    vec1d tmp(v1.dims);
    for (uint_t i : range(v1)) {
        tmp[i] = myfunc(sqrt(v1[i]));
    }

    vec1d v2 = std::move(tmp);

This is only possible using expression templates, which phy++ does not currently support for the sake of simplicity. Therefore, if performances are critical you may want to write the loop explicitly (following the guidelines in :ref:`Indexing` for optimal performance). An cleaner alternative is to use ``vectorize_lambda_first()``, which transforms a lambda function into a functor with overloaded call operator that works on both vector and scalar values. It also supports the optimization for chained calls. Contrary to the ``PHYPP_VECTORIZE()`` macro, ``vectorize_lambda_first()`` can be called in any scope, including inside other functions:

.. code-block:: c++

    auto chained = vectorize_lambda_first([](float f) { return myfunc(sqrt(f)); });

    vec1d v1 = {1.0, 1.2, 1.5};
    vec1d v2 = chained(v1);

Both ``PHYPP_VECTORIZE()`` and ``vectorize_lambda_first()`` will vectorize the function/lambda on the *first* argument only. Other arguments will simply be forwarded to all the calls, so ``foo(v,w)`` will call ``foo(v[i],w)`` for each index ``i`` in ``v``.

If instead you need to call ``foo(v[i],w[i])``, you should use ``vectorize_lambda()``. This is an alternative implementation that will support vector or scalars for *all* its arguments, and will assume that the vectors all have the same size and should be jointly iterated. The downside of this implementation is that the chaining optimization is not available.
