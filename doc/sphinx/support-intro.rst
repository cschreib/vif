Introduction
============

Here we describe the set of helper functions that are part of the phy++ support library. These functions are not *essential* to the use of the library, in that the vector and view classes can be used on their own, but they are definitely *useful*. The functions in the support library are arranged inside modular components, which one may choose to use or not. This is the first repository in which one should look for existing, extensively tested functions, which are written to be as generic and efficient as possible.

In this documentation, all these functions are sorted into categories to help you discover new functions and algorithm. Alternatively, if you know the name of a function and would like to read its documentation, you may use the search feature.

In all the following sections, each function is presented and described separately, with its multiple overloads (if any) or "sibling" functions which share a similar role. Usage examples are also provided.

The signature of each function is given in a simplified form, both for conciseness and readability. In particular, we do not display the various template meta-programming tricks used to build the function's interface (``std::enable_if<>``), only the basic type of the arguments. The documention should therefore make it clear what the function's interface is, and specify what type of arguments are allowed and rejected (i.e., are views allowed? are non-arithmetic vector types rejected?). In that, we roughly follow the conventions of `cppreference.com <http://en.cppreference.com/w/>`_.


Interface conventions
---------------------

To preserve consistency and help users ``guess`` the correct name or behavior of a function without having to look at the documentation all the time, there are a few interface conventions that all functions and types in phy++ should adhere to. Not all of these conventions are strictly enforced, and exceptions are allowed when there is a good justification.

* **Case.** All names must use lower case letters, with words separated by underscores ``_``. For example ``to_string()`` is good, but not ``ToString()`` or ``TO_STRING()``. Capital letters are reserved for macros only. This rule will be strictly enforced.
* **Language.** All names must be written in US English. This rule will be strictly enforced.
* **Spelling.** Whenever reasonable, words should not be abbreviated or truncated. Acronyms should be used sparingly. For example ``to_string()`` is good, but not ``to_str()``. Counter example: ``sexagesimal_to_degrees()`` is long to type and one has to remember it is ``degrees`` and not ``degree``; this could be reasonably abbreviated to ``sex_to_deg()``. But ``sex2deg()`` should be avoided.
* **Actions vs. tests.** Functions performing an "action" (to modify or create things) should include the verb of that action as the first word. For example ``make_mask()``, ``mask_circular_area()``, or ``begin_calculation()``. Functions performing "tests" (by returning a ``bool``) must start with the testing verb conjugated to the 3rd person. For example ``is_finite()`` or ``begins_with_prefix()``, and not ``finite()`` or ``begin_with_prefix()``, which could both indicate an "action" function instead. The counter example here is ``vec::empty()``, which is a test that should have been spelled ``vec::is_empty()``. Unfortunately the spelling ``empty()`` is used in the C++ standard library for *all* the containers, and therefore we choose to follow this spelling since most C++ programmers will be expecting it.
* **Classes vs. functions.** The naming convention is the same for classes, structs, or functions. Class names are not required to start with ``C`` (as in ``CObject``), in fact this style is discouraged.
* **Class vs. struct.** Whenever possible and reasonable, use a ``struct`` over a ``class``. The underlying semantic is: keep member values public for easy inspection by the user, make it clear what their purpose is and if/how/when they should be modified, and avoid excessive abstraction (no getters/setters). This may not apply to complex classes, or classes which must handle an insecure resource (like a raw pointer from a C library), in which case a ``class`` semantic (with private interface, etc.) is more adequate.
* **Member values.** Member values of a class or struct should follow the same naming conventions as functions. If part of the private interface of a class, they may end with an underscore to highlight that they are not part of the public interface (for example: ``cached_``).
* **Free functions vs. member functions.** There is no definite convention, but favor free functions for implementing elaborate algorithms and data processing that use this class/struct, and prefer member functions for simple modifications or queries of the state of that class/struct.
* **Variables and function arguments.** There is no convention for the names of local variables, although lower case should usually be preferred.
* **Namespaces.** All functions and classes must reside in the ``phypp`` namespace. Nested namespaces are allowed and encouraged.


