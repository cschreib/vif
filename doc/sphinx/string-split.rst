Split and combine
=================

Defined in header ``<vif/utility/string.hpp>``.


split, split_any_of
-------------------

.. code-block:: c++

    vec1s split(const std::string& s, const std::string& p); // [1]

    vec1s split_if_any_of(const std::string string& s, const std::string& c); //[2]

Function [1] will split the string ``s`` into a vector of sub-strings each time the pattern ``p`` occurs. If no such pattern is found in ``s``, the function returns a vector containing a single element which is the whole string ``s``. It should be used to parse lists of values separated by a fixed pattern, like a coma (``','``).

Function [2] will split the string ``s`` into a vector of sub-strings each time any of the characters listed in ``c`` occurs. If no such character is found in ``s``, the function returns a vector containing a single element which is the whole string ``s``. It should be used to isolate values that are separated by a variable amount of characters, such as spaces.

**Example:**

.. code-block:: c++

    vec1s str = split("foo,bar,bob", ",");
    str; // {"foo", "bar", "bob"};

    // Difference between split and split_any_of:
    std::string s;

    s = "  this is the   end";
    str = split(s, " ");
    str; // {"", "", this", "is", "the", "", "", "end"};
    str = split_any_of(s, " ");
    str; // {"this", "is", "the", "end"};

    s = "foo, bar ,  bob";
    str = split(s, ",");
    str; // {"foo", " bar ", "  bob"};
    str = split_any_of(s, ", ");
    str; // {"foo", "bar", "bob"};

    // Use case: split a line of text into words
    str = split_any_of(/* ... */, " \t\n\r");

cut
---

.. code-block:: c++

    vec1s cut(const std::string& s, uint_t n); // [1]

This function will split the string ``s`` into a vector of sub-strings (or "lines") that are exactly ``n`` characters long (except possibly the last line, which may be shorter). Contrary to the function ``wrap()``, this function does not care about spaces and word boundaries.

**Example:**

.. code-block:: c++

    vec1s str = cut("this is the end", 5);
    str; // {"this ", "is th", "e end"};


wrap
----

.. code-block:: c++

    vec1s wrap(const std::string& s, uint_t w, const std::string& i = "", bool e = false); // [1]

This function will split the string ``s`` into a vector of sub-strings (or "lines") that are at most ``w`` characters long. Contrary to the function ``cut()``, this function takes care of not splitting the line in the middle of a word. When this would occur, the cut is shifted back to just before the beginning of the word, and that word is flushed to the next line. If a word is larger than ``w``, then it will be left alone on its own line. Alternatively, if ``e`` is set to ``true``, the word is truncated and the last characters are lost and replaced by an ellipsis (``"..."``) to notify that the word has been truncated. Finally, the parameter ``i`` can be used to add indentation: these characters are added at the beginning of each line and are taken into account when calculating the line length. In this case, the first line is *not* automatically indented, to allow using a different header. This function is useful to display multi-line messages on the terminal.

**Example:**

.. code-block:: c++

    std::string str = "This is an example text with many words. Just "
        " for the sake of the example, we are going to write a "
        "veryyyyyyyyyyyyyyyyyyyyyyyyyy long word.";

    vec1s s = wrap(str, 23);
    s[0]; // "This is an example text"
    s[1]; // "with many words. Just "
    s[2]; // "for the sake of the"
    s[3]; // "example, we are going"
    s[4]; // "to write a"
    s[5]; // "veryyyyyyyyyyyyyyyyyyyyyyyyyy"
    s[6]; // "long word."

    vec1s s = wrap(str, 23, "", true);
    s[0]; // "This is an example text"
    s[1]; // "with many words. Just "
    s[2]; // "for the sake of the"
    s[3]; // "example, we are going"
    s[4]; // "to write a"
    s[5]; // "veryyyyyyyyyyyyyyyyy..."
    s[6]; // "long word."

    vec1s s = wrap(str, 23, "  ", true);
    s[0]; // "This is an example text"
    s[1]; // "  with many words. Just"
    s[2]; // "  for the sake of the"
    s[3]; // "  example, we are going"
    s[4]; // "  to write a"
    s[5]; // "  veryyyyyyyyyyyyyyy..."
    s[6]; // "  long word."


collapse
--------

.. code-block:: c++

    template<std::size_t D>
    std::string collapse(vec<D,std::string> v, const std::string& s = ""); // [1]

This function will concatenate together all the strings present in the vector ``v`` to form a single string. A separator can be provided using the argument ``s``, in which case the string ``s`` will be inserted between each pair of strings of ``v`` before concatenation.

**Example:**

.. code-block:: c++

    vec1s v = {"a", "b", "c"};
    std::string s = collapse(v);
    s; // "abc"

    s = collapse(v, ", ");
    s; // "a, b, c"
