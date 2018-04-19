.. _Formatting:

Formatting
==========

Defined in header ``<phypp/utility/string.hpp>``.


trim
----

.. code-block:: c++

    std::string trim(std::string s, const std::string& c = " \t"); // [1]

    template<std::size_t D>
    vec<D,std::string> trim(vec<D,std::string> s, const std::string& c = " \t"); // [2]

The function [1] will look at the *beginning* and *end* of the string ``s`` for any of the characters that is present in ``c`` (order is irrelevant), and remove them. This procedure is repeated until no such character is found. The net effect of this function is that the provided string ``s`` is *trimmed* from any of the characters listed in ``c``. This is useful for example to remove leading and trailing spaces of a string (which is what the default value of ``c`` does), or to removes quotes, leading zeroes, etc. The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    vec1s str = {"", "abc", " a b", " a b c  "};
    vec1s t = trim(str, " "); // trim spaces
    t; // {"", "abc", "a b", "a b c"}

    str = {"", "(a,b)", "((a,b),c)"};
    t = trim(str, "()"); // trim parentheses
    t; // {"", "a,b", "a,b),c"}


to_upper, to_lower
------------------

.. code-block:: c++

    std::string to_upper(std::string s); // [1]

    std::string to_lower(std::string s); // [2]

    template<typename T>
    vec<D,std::string> to_upper(vec<D,std::string> s); // [3]

    template<typename T>
    vec<D,std::string> to_lower(vec<D,std::string> s); // [4]

These functions will transform all characters of the string to be upper case ([1]) or lower case ([2]). It has no effect on non-alphabetic characters such as numbers, punctuation, of special characters. Functions [3] and [4] are the vectorized versions of [1] and [2], respectively.

**Example:**

.. code-block:: c++

    vec1s str = {"", "abc", "AbCdE", "No, thanks!"};
    vec1s t = to_upper(str);
    t; // {"", "ABC", "ABCDE", "NO, THANKS!"}
    t = to_lower(str);
    t; // {"", "abc", "abcde", "no, thanks!"}


align_left, align_right, align_center
-------------------------------------

.. code-block:: c++

    std::string align_left(std::string s, uint_t w, char f = ' '); // [1]

    std::string align_right(std::string s, uint_t w, char f = ' '); // [2]

    std::string align_center(std::string s, uint_t w, char f = ' '); // [3]

These functions will pad the provided string with the character ``f`` (default to a space) so that the total width the returned string is equal to ``w``. If the provided string is larger than ``w``, it is returned untouched. Padding characters will be appended at the end of the string ([1]), at the beginning of the string ([2]), or equally to both ([3]), so the string will be aligned left, right, and centered, respectively.

**Example:**

.. code-block:: c++

    std::string s = "5.0";
    std::string n = align_left(s, 6);
    n; // "5.0   "
    n = align_right(s, 6);
    n; // "   5.0"
    n = align_center(s, 6);
    n; // " 5.0  "

    // Another padding character can be used
    n = align_left(s, 6, '0');
    n; // "5.0000"
