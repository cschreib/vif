Generic functions
=================

The vector and view classes are useful on their own. However, there are a number of tasks that one may needs to do quite often, like generating a sequence of indices, or sorting a vector, which would be tedious to rewrite each time they are needed. For this reason, the phy++ support library comes with a large set of utility functions to traverse, sort, rearrange, and select data inside vectors. In this section we list these various functions and describe the corresponding algorithms.

The support library introduces a global constant called ``npos``. This is an unsigned integer whose value is the largest possible integer that the ``uint_t`` type can hold. It is meant to be used as an error value, or the value to return if no valid integer value would make sense. It is very similar in concept to the ``std::string::npos`` provided by the C++ standard library. In particular, it is worth noting that converting ``npos`` to a *signed* integer produces the value ``-1``.

We now describe the functions provided by this module, sorted by categories.

.. toctree::
    :maxdepth: 2
    :caption: Categories

    generic-range
    generic-indices
    generic-sequence
    generic-rearrange
    generic-find
    generic-error
