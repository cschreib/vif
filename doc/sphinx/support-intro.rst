Introduction
============

Here we describe the set of helper functions that are part of the phy++ support library. These functions are not *essential* to the use of the library, in that the vector and view classes can be used on their own, but they are definitely *useful*. The functions in the support library are arranged inside modular components, which one may choose to use or not. This is the first repository in which one should look for existing, extensively tested functions, which are written to be as generic and efficient as possible.

In this documentation, all these functions are sorted into categories to help you discover new functions and algorithm. Alternatively, if you know the name of a function and would like to read its documentation, you may use the search feature.

In all the following sections, each function is presented and described separately, with its multiple overloads (if any) or "sibling" functions which share a similar role. Usage examples are also provided.

The signature of each function is given in a simplified form, both for conciseness and readability. In particular, we do not display the various template meta-programming tricks used to build the function's interface (``std::enable_if<>``), only the basic type of the arguments. The documention should therefore make it clear what the function's interface is, and specify what type of arguments are allowed and rejected (i.e., are views allowed? are non-arithmetic vector types rejected?). In that, we roughly follow the conventions of `cppreference.com <http://en.cppreference.com/w/>`_.
