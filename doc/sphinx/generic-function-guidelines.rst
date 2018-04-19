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
    void foo(const T& v) {
        print(v(1,_)*2.5);
    }

    vec2f v1 = {{1.0}, {5.0}, {-1.0}}; // float, shape 3x1
    foo(v1); // prints {12.5}

    vec2i v2 = {{1, 5}, {-1, 1}}; // int, shape 2x2
    foo(v2); // prints {-2.5, 2.5}

.. note:: This example uses the ``print()`` function from the phy++ support library, which simply displays its arguments on the terminal.


Summary of guidelines
---------------------

The functions we use as examples here can be somewhat silly, but they will serve to illustrate a number of important "rules" which one should follow when writing generic functions. These rules are explained in detail below, and can be summarized as follows:

* Use the most specific type possible for the function arguments (e.g., ``vec<1,T>`` instead of just ``T``).
* Even when the function is supposed to work on only one specific data type, leave the data type of vectors unconstrained in order to support both vectors and view (e.g., ``vec<D,T>`` instead of ``vec<D,int>``).
* Use ``std::enable_if<>`` to express any remaining constraints on the type.
* Do not use the ``T`` in ``vec<D,T>`` to form new variables, use ``meta::rtype_t<T>`` instead.
* Provide default values for template arguments whenever it makes sense, to enable support for initializer lists.
* Avoid output or input/output parameters whenever possible, else use universal references ``T&&`` and ``std::enable_if<>`` as described in the guidelines below.
* Use the available helper tools to vectorize existing functions.
* Use ``phypp_check()`` to express any constraints on the data that can only be checked at run time (number of elements, value ranges, etc.).


Expressing constraints on function arguments
--------------------------------------------

In the example above, the type of the function's argument ``v`` is ``T``, and is totally unconstrained. It could be anything. However this specific function has an *implicit* constraint on the type ``T``: it must be possible to write ``v(1,_)*2.5``. This means ``v`` *must* be a 2D vector, and the elements must be of arithmetic type. If you try to use this function on a value which does not satisfy this implicit constraint, the compiler will generate an error:

.. code-block:: c++

    vec1f v4 = {1.0, 5.0, 6.0};
    foo(v4); // error! 1D vector

The program will not compile, which is good. However the error message will be *nasty* (180 lines of errors with Clang), and the error will point to code *inside* the function. This is far from ideal, because the user of the function should not need to understand what the function does internally to fix the error. This can be fixed by adding constraints on the type ``T``:

.. code-block:: c++

    template<typename T>
    void foo(const vec<2,T>& v) {
        print(v(1,_)*2.5);
    }

Here, we state that ``v`` must be a 2D vector, and we leave the type of the elements unconstrained. Using this new definition of ``foo()``, the error message in the case above becomes much smaller (4 lines of errors), and explicitly says that there is no matching function call for ``foo(v4)``. This is much better. Therefore you should always make sure to specify as many constraints as possible in the *signature* of the function (i.e., the type of its arguments).

We are not done though. Indeed, we still left the type of the elements unconstrained, while we need elements of arithmetic types to be able to write ``v(1,_)*2.5``. For example, using a vector of strings would be an error:

.. code-block:: c++

    vec2s v5 = {{"I", "am"}, {"a", "string"}};
    foo(v5); // error! vector of strings

Again, this emits a lengthy error from inside the function (20 lines of errors). We can fix this by adding extra constraints on the type ``T`` of the elements. One possibility is to force it to be some "common" type, like ``double``:

.. code-block:: c++

    // Note: parameter is fully constrained, it is not generic anymore
    void foo(const vec<2,double>& v) {
        print(v(1,_)*2.5);
    }

This makes the error much easier to understand (7 lines of errors), but it has the important downside that the function is no longer generic: it *needs* a vector of ``double``. If you try to call it on a vector of ``float``, it will first have to make a copy of that vector and convert all values to ``double`` before calling the function, which is not optimal. It will also fail to work on *views* (see below). So unless you know the function should only be used with ``double`` values, this is not the right solution. Instead, we can leave the type of elements to be generic, and use ``std::enable_if<>`` to express a constraint on this type, in this case ``std::is_arithmetic<T>``:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void foo(const vec<2,T>& v) {
        print(v(1,_)*2.5);
    }

With this version of the function, the error when called on vectors of strings becomes much clearer (4 lines of errors) and says that you cannot call the function on strings. Again, much better!

So, that's it? Not quite. There is one last implicit requirement when we write ``v(1,_)``: the first dimension of ``v`` must have at least two elements. There is no way to check this at the time of compilation, so the faulty program below will compile:

.. code-block:: c++

    vec2i v6;
    foo(v6); // compiles, but runtime error! empty vector

It will fail at runtime though. The backtrace will show that the error happened in ``foo()``, but with a rather cryptic error message:

.. code-block:: c++

    error: operator(): index out of bounds (1 vs. 0)

The solution here is to perform an explicit check inside the function, and emit a clearer error message using the ``phypp_check()`` function:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void foo(const vec<2,T>& v) {
        phypp_check(v.dims[0] >= 2, "vector must have at least two elements along first dimension ",
            "(got ", v.dims[0], ")");
        print(v(1,_)*2.5);
    }

The error shown to the user then becomes clear:

.. code-block:: c++

    error: vector must have at least two elements along first dimension (got 0)

Now that we do an explicit check that the index ``1`` is valid before accessing the vector, we no longer need the vector to perform automatic bounds checking. Therefore we can use the "safe" indexing interface:

.. code-block:: c++

    template<typename T,
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void foo(const vec<2,T>& v) {
        phypp_check(v.dims[0] >= 2, "vector must have at least two elements along first dimension ",
            "(got ", v.dims[0], ")");
        print(v.safe(1,_)*2.5);
    }

This is the optimal way to write this function, and it is clearly not as pretty as the very first version. This shows that, while writing generic functions is easy, writing them *well* is much harder. For this reason, always check in the support library if a function already exists before writing your own.

It should be said, however, that the very first version we wrote actually does the work we expect it to do (save for the fact that it does not support initializer lists, see below). It is not "incorrect"; its only defect is that it will not be very helpful when things go wrong.


Supporting initializer lists
----------------------------

There is one last modification we can do to make the ``foo()`` function "as good as it gets". Indeed, even with the last version, we cannot use initializer lists directly as function arguments:

.. code-block:: c++

    foo({{1, 5}, {-1, 1}});

This generates an error because the compiler is not smart enough to infer the type ``T`` of the vector from this initializer list. Unfortunately, in general we cannot do this *perfectly* and support any type in the initializer list.

But we can still make it work. The trick is to specify a default value for the template parameter ``T``, for example ``double``. This way, the initializer list will automatically be used to form a vector of ``double``, and the code will compile and run. This is not a perfect solution because the *true* type of the values in the initializer list is lost, but in most cases it is possible to identify a "safe" type (such as ``double``) that will be able to do the job anyway.

In this particular case, ``double`` is actually a perfect choice because we multiply the values of the vector by ``2.5``, which requires a conversion to ``double`` anyway. So converting the values of the initializer list to ``double`` will not change the final result. The definite, final version of our function is thus:

.. code-block:: c++

    template<typename T = double, // use a default value here
        typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    void foo(const vec<2,T>& v) {
        phypp_check(v.dims[0] >= 2, "vector must have at least two elements along first dimension ",
            "(got ", v.dims[0], ")");
        print(v.safe(1,_)*2.5);
    }


Supporting both vectors and views
---------------------------------

For a vector of type ``vec<D1,T>``, a view will have a type ``vec<D2,T*>`` or ``vec<D2,const T*>``. The number of dimensions can be different, and the data type is a pointer to the type of elements in the vector. The ``const`` qualifier is used to propagate const-correctness if the original vector was declared ``const``.

This makes it relatively easy to write function that work on both vectors and view, but this distinction means that there are a number a details to keep in mind. Consider this generic function that computes the sum of all the elements in a vector:

.. code-block:: c++

    template<std::size_t D = 1, typename T = double>
    T sum_it_all(const vec<D,T>& v) {
        T ret = 0;
        for (const T& val : v) {
            ret += val;
        }

        return ret;
    }

.. note:: Such a function already exists in the phy++ support library, and is called ``total()`` (for integers and floating point values) or ``count()`` (for boolean values). Their return type is determined in a smarter way than we discuss here, to prevent overflow and underflow.

This implementation works for all vectors, but it will fail for views. Indeed, if called on a view of type ``vec<1,int*>``, then ``T = int*``, and the return value is not an integer but an (invalid!) pointer to an integer. Fortunately, it will not even compile because the loop will try to assign the values of ``v`` to a ``const int*&``, which will fail. Therefore, the type ``T`` should never be used directly like this.

Instead, you should apply the transform ``meta::rtype_t<T>``, which essentially transforms ``T*`` into ``T`` and removes const qualifiers, and use ``auto`` whenever possible to let the type system make the right decisions for you:

.. code-block:: c++

    template<std::size_t D = 1, typename T = double>
    meta::rtype_t<T> sum_it_all(const vec<D,T>& v) {
        meta::rtype_t<T> ret = 0;
        for (const auto& val : v) {
            ret += val;
        }

        return ret;
    }

There are a few, rarer corner cases to keep in mind when both view and vectors need to be supported. The case of output parameters, in particular, is described further below.


Vectorizing scalar functions
----------------------------

Most function created in C++ thus far, including those in the C++ standard library, are *scalar* functions which operate on one single value. The best example of this are all the mathematical functions, ``sqrt()``, ``pow()``, ``ceil()``, etc. These functions can be *vectorized* to operate directly on vector data without having to write a loop. The phy++ support library contains a large number of such vectorized functions:

.. code-block:: c++

    double v1 = 2.0;
    sqrt(v1); // 1.41...
    vec1d v2 = {2.0, 4.0, 6.0};
    sqrt(v2); // {1.41..., 2.0, 2.45...}

However the phy++ support library cannot contain *all* functions that ever existed, and you may create your own scalar functions that you wish to vectorize. This can be achieved using the preprocessor macro ``PHYPP_VECTORIZE()``:

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

This is only possible using expression templates, which currently phy++ does not support for the sake of simplicity. Therefore, if performances are critical you may want to write the loop explicitly (following the guidelines in :ref:`Indexing` for optimal performance). An cleaner alternative is to use ``vectorize_lambda_first()``, which transforms a lambda function into a functor with overloaded call operator that works on both vector and scalar values. It also supports the optimization for chained calls. Contrary to the ``PHYPP_VECTORIZE()`` macro, ``vectorize_lambda_first()`` can be called in any scope, including inside other functions:

.. code-block:: c++

    auto chained = vectorize_lambda_first([](float f) { return myfunc(sqrt(f)); });

    vec1d v1 = {1.0, 1.2, 1.5};
    vec1d v2 = chained(v1);

Both ``PHYPP_VECTORIZE()`` and ``vectorize_lambda_first()`` will vectorize the function/lambda on the *first* argument only. Other arguments will simply be forwarded to all the calls, so ``foo(v,w)`` will call ``foo(v[i],w)`` for each index ``i`` in ``v``.

If instead you need to call ``foo(v[i],w[i])``, you should use ``vectorize_lambda()``. This is an alternative implementation that will support vector or scalars for *all* its arguments, and will assume that the vectors all have the same size and should be jointly iterated. The downside of this implementation is that the chaining optimization is not available.


Output arguments and views
--------------------------

In general, the only output of a function must be its return value. Output arguments should only be used when: a) the function must return multiple values, and b) it would be inefficient or impractical to return them by value. Otherwise, one may wish to use input/output arguments for functions that have no return value but only modify the content of an existing vector. As you will see below, writing functions with output or input/output vector arguments is possible but a bit nasty, so make sure you really need them before diving in.

The typical example where output arguments are needed is the following function which converts a string to a value of another type (e.g., an integer):

.. code-block:: c++

    template<typename T>
    bool from_string(const std::string& s, T& v) {
        std::istringstream ss(s);
        return ss >> v;
    }

.. note:: In C++ there is no difference between purely output parameters (only used to store a result) and input/output parameters (used to read data and write results back). As a result, even though the discussion here is centered on output parameters, the same principles apply to input/output parameters as well.

This function returns a flag to let the user know whether the conversion was successful, and the output value is stored in the argument ``v``, which is a reference (``T&``). The function is then used as follows:

.. code-block:: c++

    int v;
    if (from_string("42", v)) {
        // do whatever with 'v'
    } else {
        error("could not convert the string");
    }

.. note:: In C++17, one may wish to return an ``std::optional<T>`` instead, which would be the optimal solution for the scalar case. However this solution does not vectorize well. Currently, ``vec<D,std::optional<T>>`` is not supported; it may work, but use it at your own risk.

The vectorization of such functions cannot be done with the automatic vectorization tools described above, so we will have to do it manually. It is rather simple, right? We only need to use a reference to an output vector ``vec<D,T>&``:

.. code-block:: c++

    template<std::size_t D = 1, typename U = std::string,
        typename T, typename enable = typename std::enable_if<
        std::is_same<meta::rtype_t<U>, std::string>::value
    >::type>
    vec<D,bool> from_string(const vec<D,U>& s, vec<D,T>& v) {
        vec<Dim,bool> res(s.dims);
        v.resize(s.dims);
        for (uint_t i : range(s)) {
            res.safe[i] = from_string(s.safe[i], v.safe[i]);
        }

        return res;
    }

In this particular case, we use ``std::enable_if<>`` to make sure the input type is either a vector of strings or a view on such vector. We then return a vector of ``bool`` so the user can check the success of the conversion for each individual value separately. The function is then used as follows:

.. code-block:: c++

    vec1s s = {"5", "-6", "9", "42"};

    vec1i v;
    vec1b r = from_string(s, v);

    for (uint_t i : range(s)) {
        if (r[i]) {
            // do whatever with 'v[i]'
        } else {
            error("could not convert the string");
        }
    }

The catch here is to support *views* as output arguments. Indeed, one may want to convert only part of a string vector with ``from_string()`` and store the result in a view, in which case our current implementation fails:

.. code-block:: c++

    vec1s s = {"5", "-6", "9", "42"};

    // Say we only want to convert values with 2 characters
    vec1u id = where(length(s) == 2);

    // This does *not* work:
    vec1i v(s.dims); // resize output vector beforehand
    vec1b r = from_string(s[id], v[id]);

    // error: 'v[id]' is an r-value, cannot bind it to a reference 'vec<D,T>&'

    // But this works:
    vec1i v(s.dims); // resize output vector beforehand
    vec1i tmp;       // create a temporary
    vec1b r = from_string(s[id], tmp);
    v[id] = tmp;     // assign temporary values to 'v'

.. note:: This issue also affects IDL, in a nastier way since IDL will not throw any error. It will store the output values in a automatically generated temporary vector, which is then discarded, so the values of the view are not modified... Oops! In IDL, this can only be solved by explicitly introducing a temporary vector and assigning it back to the view, as done in the example above. But C++ is smarter, and we can make this work! Read on.

Such type of problem arises whenever you write a function that takes a non-constant reference to a vector in order to modify its values. To support this type of usage with views, we need an argument type that can be either an "l-value" (a reference to a vector) or an "r-value" (a temporary view). This is exactly what the "universal reference" is for: ``T&&``. Unfortunately, this universal reference requires an unconstrained type ``T``. This means we loose all the implicit constraints on the type: it is no longer ``vec<D,T>``, but simply ``T``. Therefore we will have to specify these constraints explicitly using ``std::enable_if<>``. And there are a lot of constraints! We want to make sure:

* that ``T`` is a vector or a view,
* that the number of dimensions of ``T`` is the same as the input vector of strings,
* that if ``T`` is an r-value, it must be a non-constant view,
* that if ``T`` is an l-value, it must be a non-constant reference (to a vector or a view).

Since these basic requirements will be the same for every vectorized function with output parameters, a specific trait is provided in phy++ to express all these constraints: ``meta::is_compatible_output_type<In,Out>``. It is used in the following way:

.. code-block:: c++

    template<std::size_t D = 1, typename U = std::string, typename T,
        typename enable = typename std::enable_if<
        std::is_same<meta::rtype_t<U>, std::string>::value && // this was there before
        meta::is_compatible_output_type<vec<D,U>,T>::value    // this is the new trait
    >::type>
    vec<D,bool> from_string(const vec<D,U>& s, T&& v) {
        // ...
    }

In addition, here we need to differentiate the behavior of the function for the two cases: we want the "vector" version to automatically resize the output vector to the dimensions of the input vector, and the "view" version to simply check that the view has the same dimensions as the input vector. This is expected to be a common behavior for functions with output parameters, therefore a helper function is provided in phy++ to do just that: ``meta::resize_or_check(v, d)``. This function resizes the vector ``v`` to the dimensions ``d``, or, if ``v`` is a view, checks that its dimensions match ``d``. The final, fully generic, vectorized function is therefore:

.. code-block:: c++

    template<std::size_t D = 1, typename U = std::string, typename T,
        typename enable = typename std::enable_if<
        std::is_same<meta::rtype_t<U>, std::string>::value &&
        meta::is_compatible_output_type<vec<D,U>,T>::value
    >::type>
    vec<D,bool> from_string(const vec<D,U>& s, T&& v) {
        vec<Dim,bool> res(s.dims);
        meta::resize_or_check(v, s.dims);
        for (uint_t i : range(s)) {
            res.safe[i] = from_string(s.safe[i], v.safe[i]);
        }

        return res;
    }

If you are in a case where there is no "input" vector to consider, and you simply want to write a function that modifies an existing vector's values (i.e., an input/output parameters), use the simpler ``meta::is_output_type<T>`` trait:

.. code-block:: c++

    template<typename T, typename enable = typename std::enable_if<
        std::is_vec<T>::value && meta::is_output_type<T>::value
    >::type>
    void twice(T&& v) {
        v[_] *= 2;
    }

This traits only checks the last two conditions of ``meta::is_compatible_output_type``, namely:

* that if ``T`` is an r-value, it must be a non-constant view,
* that if ``T`` is an l-value, it must be a non-constant reference (to a vector or a view).


Creating views
--------------

TODO


Metaprogramming helpers
-----------------------

TODO
