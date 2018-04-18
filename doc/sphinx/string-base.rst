Basic string operations
=======================

Defined in header ``<phypp/utility/string.hpp>``.

empty
-----

.. code-block:: c++

    bool empty(const std::string& s); // [1]

    template<std::size_t D>
    vec<D,bool> empty(const vec<D,std::string>& s); // [2]

The function [1] will return ``true`` if the provided string does not contain *any* character (including spaces), and ``false`` otherwise. This is a synonym for ``s.empty()``. The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    vec1s str = {"", "abc", "   "};
    vec1b b = empty(str);
    b; // {true, false, false}
    // Not to be confused with the vec::empty() function
    str.empty(); // false
    str = {""};
    str.empty(); // false
    str = {};
    str.empty(); // true


length
------

.. code-block:: c++

    uint_t length(const std::string& s); // [1]

    template<std::size_t D>
    vec<D,uint_t> length(const vec<D,std::string>& s); // [2]

The function [1] will return the length of the provided string, i.e., the number of character it contains (including spaces). If the string is empty, the function will return zero. The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    vec1s str = {"", "abc", " a b"};
    vec1u n = length(str);
    n; // {0, 3, 4}


keep_first, keep_last
---------------------

.. code-block:: c++

    std::string keep_first(std::string s, uint_t n = 1); \\ [1]

    std::string keep_last(std::string s, uint_t n = 1); \\ [2]

    template<std::size_t D>
    vec<D,std::string> keep_first(vec<D,std::string> s, uint_t n = 1); \\ [3]

    template<std::size_t D>
    vec<D,std::string> keep_last(vec<D,std::string> s, uint_t n = 1); \\ [4]

These functions will return the first ([1]) or last ([2]) ``n`` characters of the string ``s`` and discard the rest. If ``n`` is larger than the size of ``s``, the whole string is returned untouched. Functions [3] and [4] are the vectorized versions of [1] and [2], respectively.

**Example:**

.. code-block:: c++

    vec1s v = {"p1_m2.txt", "p3_c4.fits", "p1_t8.fits"};
    vec1s s = keep_first(v, 2);
    s; // {"p1", "p3", "p1"}
    s = keep_last(v, 4);
    s; // {".txt", "fits", "fits"}


distance
--------

.. code-block:: c++

    uint_t distance(const std::string& s1, const std::string& s2); // [1]

    template<std::size_t D>
    vec<D,uint_t> distance(const vec<D,std::string>& s1, const std::string& s2); // [2]


The function [1] computes the *lexicographic distance* between two strings. The definition of this distance is the following. If the two strings are exactly identical, the distance is zero. Else, each character of the shortest string are compared to the corresponding character at the same position in the other string: if they are different, the distance is increase by one. Finally, the distance is increased by the difference of size between the two strings.

The goal of this function is to identify *near* matches in case a string could not be found in a pre-defined list. This is useful to suggest corrections to the user, who may have misspelled it.

**Example:**

.. code-block:: c++

    vec1s s = {"wircam_K", "hawki_Ks", "subaru_B"};
    vec1u d = distance(s, "wirkam_Ks");
    d; // {2, 8, 8}

    // Nearest match
    std::string m = s[min_id(d)];
    m; // "wircam_K"

