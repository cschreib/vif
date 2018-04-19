Hash
====

Defined in header ``<phypp/utility/string.hpp>``.


hash
----

.. code-block:: c++

    template<typename ... Args>
    std::string hash(const Args& ... args); // [1]

This function scans all the arguments that are provided, and returns the hexadecimal representation of the SHA-1 "hash" of this argument list. The hash is a string such that: 1) all further calls of ``hash(...)`` with arguments that have the exact same value (perhaps when the program is executed a second time later) will always return the same string, and 2) the probability is very small that the function returns the same string for another set of arguments, or arguments with different values. Although this algorithm was created in 1995, the first "collision" (two different data sets producing the same hash) was found in 2017.

This is useful for example to cache the result of some heavy computation: once the computation is done, the *input* parameters of the computation can be fed to ``hash()`` to give a "sort-of-unique" identifier to the "input+result" pair. The result of the computation can then be saved somewhere with the hash as an identifier. Later on, if the computation is requested with a new set of parameters, these parameters are fed to ``hash()`` and the resulting string is compared to all the identifiers of the cached results: if a match is found, then the associated pre-computed result can be re-used, else the computation must be executed anew.

**Example:**

.. code-block:: c++

    std::string s;

    // With a single argument
    s = hash("hello world!");
    s; // "da52a1357f3c973e1ffc1b694d5308d0abcd9845"
    s = hash("hello world?")
    s; // "793e673d04e555f8f0b38033d5223c525a040719"
    // Notice how changing a single character gives a completely
    // different hash string

    // With multiple arguments
    s = hash(1, 2, 3);
    s; // "570331ab965721aae8a8b3c628cae57a21a37560"
    s = hash("123");
    s; // "0e898437b29ec20c39ca48243e676bcb177d4632"
    s = hash(1.0, 2.0, 3.0);
    s; // "9c45014f7c7943cb7860f3db4b885fb44b510ec8"
    // Notice how the hash is different even though we would
    // consider these different sets of values to be equivalent.
