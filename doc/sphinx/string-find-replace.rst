Find/replace
============

Defined in header ``<phypp/utility/string.hpp>``.

find
----

.. code-block:: c++

    uint_t find(const std::string& s, const std::string& p); // [1]

    template<std::size_t D>
    vec<D,uint_t> find(const vec<D,std::string>& s, const std::string& p); // [2]

The function [1] returns the position in the string ``s`` of the first occurrence of the sub-string ``p``. If no such occurrence is found, the function returns ``npos``. The function [2] is the vectorized version.

**Example:**

.. code-block:: c++

    vec1s v = {"Apple", "please", "complementary", "huh?"};
    vec1u p = find(v, "ple");
    p; // {2, 0, 3, npos}


replace
-------

.. code-block:: c++

    std::string replace(std::string s, const std::string& p, const std::string& r); // [1]

    template<std::size_t D>
    vec<D,std::string> replace(vec<D,std::string> s, const std::string& p, const std::string& r); // [2]

The function [1] will look inside the string ``s`` for occurrences of the pattern ``p`` and replace each of them with the replacement ``r``. The string is unchanged if no occurrence is found. In particular, this function can also be used to remove all the occurrences of ``p`` simply by setting ``r`` equal to an empty string. The function [2] is the vectorized version.

**Example:**

.. code-block:: c++

    vec1s str = {"I eat apples", "The apple is red"};
    vec1s r = replace(str, "apple", "pear");
    r; // {"I eat pears", "The pear is red"}

    str = {"a:b:c", "g::p"};
    r = replace(str, ":", ",");
    r; // {"a,b,c", "g,,p"};


begins_with, ends_with
----------------------

.. code-block:: c++

    bool begins_with(const std::string& s, const std::string& p); // [1]

    bool ends_with(const std::string& s, const std::string& p); // [2]

    vec<D,bool> begins_with(const vec<D,std::string>& s, const std::string& p); // [3]

    vec<D,bool> ends_with(const vec<D,std::string>& s, const std::string& p); // [4]

These functions will return ``true`` if the beginning ([1]) or the end ([2]) of the string ``s`` exactly matches the string ``p``. The functions [3] and [4] are the vectorized versions.

**Example:**

.. code-block:: c++

    vec1s v = {"p1_m2.txt", "p3_c4.fits", "p1_t8.fits"};
    vec1b b = begins_with(v, "p1");
    b; // {true, false, true}

    b = ends_with(v, ".fits");
    b; // {false, true, true}


erase_begin, erase_end
----------------------

.. code-block:: c++

    std::string erase_begin(std::string s, uint_t n); // [1]

    std::string erase_begin(std::string s, const std::string& p); // [2]

    std::string erase_end(std::string s, uint_t n); // [3]

    std::string erase_end(std::string s, const std::string& p); // [4]

    vec<D,std::string> erase_begin(vec<D,std::string> s, uint_t n); // [5]

    vec<D,std::string> erase_begin(vec<D,std::string> s, const std::string& p); // [6]

    vec<D,std::string> erase_end(vec<D,std::string> s, uint_t n); // [7]

    vec<D,std::string> erase_end(vec<D,std::string> s, const std::string& p); // [8]

These functions will erase characters from the beginning ([1] and [2]) or the end ([3] and [4]) of the string ``s``.

Functions [1] and [3] will remove ``n`` characters. If ``n`` is larger than the size of ``s``, the returned string will be empty. Functions [2] and [4] first check that the string begins or ends with the other string ``p`` provided as second argument: if it does, it removes this substring from ``s``; if it does not, an error is reported and the program stops.

Functions [5] to [8] are the vectorized versions of functions [1] to [4], respectively.

**Example:**

.. code-block:: c++

    vec1s v = {"p1_m2.txt", "p3_c4.fits", "p1_t8.fits"};
    std::string s = erase_begin(v[0], "p1_");
    s; // "m2.txt"
    s = erase_begin(v[1], "p1_");
    // will trigger an error
    s = erase_begin(v[2], "p1_");
    s; // "t8.fits"

    s = erase_end(v[0], ".fits");
    // will trigger an error
    s = erase_end(v[1], ".fits");
    s; // "p3_c4"
    s = erase_end(v[2], ".fits");
    s; // "p1_t8"

    vec1s t = erase_begin(v, 3);
    t; // {"m2.txt", "c4.fits", "t8.fits"}
    t = erase_end(v, 5);
    t; // {"p1_m", "p3_c4", "p1_t8"}


replace_block, replace_blocks
-----------------------------

.. code-block:: c++

    template<typename T>
    std::string replace_block(std::string v, const std::string& b, const std::string& e, T f); // [1]

    template<typename T>
    std::string replace_blocks(std::string v, const std::string& b, const std::string& s, const std::string& e, T f); // [2]

Function [1] looks in the string ``v`` and identifies all "blocks" that start with ``b`` and end with ``e``. The content of each block is fed to the user-supplied function ``f`` which does any kind of conversion or operation on that content, and must returns a new string as a replacement. This new string is then inserted in ``v`` and replaces the entire block.

Function [2] does the same thing, except that each block can have multiple "components" that are separed with the separator ``s``. In this case, the function extracts all these "components" and stores them inside a string vector, and feeds this vector to the conversion function ``f``.

See ``regex_replace`` for a more powerful (but also more complex and possibly slower) alternative.

**Example:**

.. code-block:: c++

    // We want to modify the content inside <b>...</b> to be upper case
    std::string s = "This is a <b>whole</b> lot of <b>money</b>";
    std::string ns = replace_block(s, "<b>", "</b>", [](std::string t) {
        return to_upper(t);
    });

    ns; // "This is a WHOLE lot of MONEY"

    // We want to convert this LaTeX link into HTML
    s = "Look at \url{http://www.google.com}{this} link.";
    ns = replace_blocks(s, "\url{", "}{", "}", [](vec1s t) {
        return "<a href=\""+t[0]+"\">"+t[1]+"</a>";
    });

    ns; // "Look at <a href="http://www.google.com">this</a> link."

