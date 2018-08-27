.. _Known issues and problems:

Known issues and problems
=========================

This section describes issues and limitations in vif. Some of these are design issues which should be solved by me, the author of the library, but which I actually *cannot* solve, or haven't found the time solve yet. Some are not issues but *features*, conscious choices I made that may surprise some of you.

In any case, all these may require you, the user, to pay special attention to some corner cases, and you should therefore make sure you are familiar with them.


Dangling views
--------------

**The problem 1.**

A "dangling" reference is a invalid reference that points to an object that no longer exists. Creating such dangling references is a common programming mistake, which compilers are fortunately well equipped to detect:

.. code-block:: c++

    int& foo() {
        int i = 0;
        return i; // warning: reference to local variable ‘i’ returned
    }

vif views suffer from the same problem: a dangling view can be created that points to a vector that no longer exists. Unfortunately, compilers are not aware of it:

.. code-block:: c++

    vec<1,int*> foo() {
        vec<1,int> v = {0, 1, 2, 3};
        return v[_];
    }

This code will unfortunately compile without warning, and there is no (efficient) programmatic way to identify it at run time. Calling ``foo()`` will therefore not throw any immediate error, but create an "undefined behavior"; it may do anything, and that's very bad.

The case above can be spotted by looking carefully at the function itself: 1) the function returns a view, 2) this view points to a vector that is created in the function, 3) therefore it's a dangling view. However there are more subtle cases where the issue is not as obvious...

**The problem 2.**

One such case where dangling views are hard to spot is when returning a view from a lambda, rather than from a function with a specified return type. In this case, the return type is *inferred* from the returned expression, and this causes some surprising results. To avoid this, C++ has some special rules for references:

.. code-block:: c++

    auto foo = []() {
        int i = 0;
        return i;
    }

This lambda actually returns *a copy* of ``i``, to avoid silently creating a dangling reference. Again, unfortunately this special behavior does not apply to vif views:

.. code-block:: c++

    auto foo = []() {
        vec1i v = {0, 1, 2, 3};
        return v[_];
    }

This lambda will return a dangling view, and because it would not do that for normal types and references like ``int`` and ``int&``, it is *not* obvious to spot.


**The solution.**

If you are a seasoned C++ programmer, you know how to avoid dangling references within the current C++ language rules. Use the same caution with views, and you will avoid most "Problem 1" (as described above) without relying on the compiler to throw warnings at you.

To solve "Problem 2", take extra precautions whenever you write a lambda function (or a function with ``auto`` return type in C++14) to ensure you are not accidentally returning a view. If you do not trust yourself to do this, then make sure you *always* specify the expected return type of your lambda functions:

.. code-block:: c++

    auto foo = []() -> vec1i {
        vec1u v = {0, 1, 2, 3};
        return v[_];
    }

Note that smart people are currently thinking of adding new C++ rules that will allow me (and other library authors who experience similar problems) to modify the view class such that it will benefit from all the good magic that C++ currently applies to references. This will fix "Problem 2", and some cases of "Problem 1". In the mean time, just be careful!


Invalid views
-------------

**The problem.**

With ``std::vector<T>``, any operation that modifies the size of the vector *invalidates* all the iterators that point to this vector:

.. code-block:: c++

    std::vector<int> v;
    auto b = v.begin();

    v.resize(10);
    // b is now invalid!

The same is true for views: if a view points to a vector and this vector is later resized or re-assigned, the view becomes invalid and *must not* be used any more.

.. code-block:: c++

    vec1u vec = {1,2,3,4};
    vec<1,int*> view = vec[_];

    v = {1,2,3,4,5,6};
    // the view is now invalid!

The reason why is that the view stores *pointers* to the values in ``vec``, not indices. These pointers may become invalid themselves if the values of ``vec`` are moved to another spot in the computer's memory.


**The solution.**

There is a reason why shortcut types are provided for vectors (``vec1i`` instead of ``vec<1,int>``) and not for views: *views are only meant to be temporaries*, they should not be saved into named variables like in the above. If you feel it is necessary to do this for performance reasons, simply avoid using views altogether and manipulate indices explicitly, this will be faster.
