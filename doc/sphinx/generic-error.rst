Error checking
==============

Defined in header ``<phypp/core/error.hpp>``.

phypp_check
-----------

.. code-block:: c++

    template<typename ... Args>
    void phypp_check(bool b, Args&& ... args);

This function makes error checking easier. When called, it checks the value of ``b``. If ``b`` is ``true``, then nothing happens and the values of ``args`` are not evaluated. However if ``b`` is ``false``, then the current file and line are printed to the standard output, followed by the other arguments of this function, which are supposed to compose an error message explaining what went wrong, and the program is immediately stopped. This function is used everywhere in the phy++ library to ensure that certain pre-conditions are met before doing a calculation, and it is *essential* to make the program stop in case something is unexpected rather than letting it run hoping for the best (and often getting the worst).

Since this function is actually implemented by a preprocessor macro, one should not worry about its performance impact beside the cost of evaluating ``b``. Performances will only be affected when something goes wrong and the program is about to stop anyway.

In addition, if the phy++ library was configured accordingly, this function can print the "backtrace" of the program that lead to the error. This backtrace lists which functions or lines of code the program was executing when the error occured, tracing back the offending line all the way up from the ``main()`` function. This can be very useful to identify the source of the problem, but is only available if debugging informations are stored inside the compiled program. Note that this only affects the size of the executable on disk: debugging informations do not alter performances. Lastly, the backtrace may not be sufficient to understand what went wrong, and one may need to use the debugger.

.. code-block:: c++

    // Suppose 'v' is read from the command line arguments.
    vec1i v = /* read from somewhere unsafe */;

    // The rest of the code needs at least 3 elements in the
    // vector 'v', so we need to check that first.
    phypp_check(v.size() >= 3, "this algorithm needs at least 3 values in the input vector, "
        "but only ", v.size(), " were found");

    // If we get past this point, we can proceed safely to use
    // the first three elements
    print(v[0]+v[1]+v[3]);

    // ... and we can safely disable index checking using the
    // 'safe' interface of the vector
    print(v.safe[0]+v.safe[1]+v.safe[3]);
