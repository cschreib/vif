Regular expressions (regex)
===========================

Defined in header ``<phypp/utility/string.hpp>``.


regex_match
-----------

.. code-block:: c++

    bool regex_match(const std::string& s, const std::string& r); // [1]

    template<std::size_t D>
    bool regex_match(const vec<D,std::string>& s, const std::string& r); // [2]

Function [1] will return ``true`` if the string ``s`` matches the regular expression (or "regex" for short) ``r``. This regular expression can be used to identify complex patterns in a string, far more advanced than just matching the existence of a sub-string. This particular implementation uses POSIX regular expressions. The syntax is complex and not particularly intuitive, but it has become a well known standard and is therefore understood by most programmers. A gentle tutorial can be found here_. If the regular expression is invalid, an error will be reported to diagnose the problem, and the program will stop. Function [2] is the vectorized version.

.. _here: http://www.zytrax.com/tech/web/regex.htm

.. note:: For regular expression we advise to use C++ *raw string literals*. Indeed, in these expressions one often needs to use the backslash character (``\``), but this character is also used in C++ to form special characters, such as ``'\n'`` (new line). For this reason, to feed the backslash to the regex compiler one actually needs to escape it: ``"\\"``. But the POSIX rules also use  ``\`` as an escape character, and so this can quickly cause head aches... (to have a regex matching the backslash character itself, one would have to write ``"\\\\"``!) To avoid this inconvenience, one can enclose the regular expression in ``R"(...)"``, and not need to worry about escaping characters from the C++ compiler.

**Example:**

.. code-block:: c++

    vec1s v = {"abc,def", "abc:def", "956,fgt", "9g5,hij", "ghji,abc"};

    // We want to find which strings have the "XXX,YYY" format
    // where "XXX" can be any combination of three letters or numbers
    // and "YYY" can be a combination of three letters
    vec1b b = regex_match(v, R"(^[a-z0-9]{3},[a-z]{3}$)");

    // This regular expression can be read like:
    //  '^'   : at the beginning of the string, match...
    //  '['   : any character among...
    //  'a-z' : letters from 'a' to 'z'
    //  '0-9' : numbers from 0 to 9
    //  ']'
    //  '{3}' : three times
    //  ','   : followed by a coma, then...
    //  '['   : any character among...
    //  'a-z' : letters from 'a' to 'z'
    //  ']'
    //  '{3}' : three times
    //  '$'   : and then the end of the string

    b[0]; // true
    b[1]; // false, there is no ","
    b[2]; // true
    b[3]; // true
    b[4]; // false, too many characters before the ","


regex_extract
-------------

.. code-block:: c++

    vec2s regex_extract(const std::string& s, const std::string& r); // [1]

This function will analyze the string ``s``, perform regular expression matching (see ``regex_match`` above) using the regular expression ``r``, and will return a vector containing all the extracted substrings. To extract one or more substrings in the regular expression, just enclose the associated patterns in parentheses. The returned vector is two dimensional: the first dimension corresponds to the number of times the regular exception was matched in the provided string, the second dimension corresponds to each extracted substring.

**Example:**

.. code-block:: c++

    std::string s = "array{5,6.25,7,28}end45.6ddfk 3.1415";

    // We want to find all floating point numbers in this mess, and
    // extract their integer and fractional parts separately.
    vec2s sub = regex_extract(s, R"(([0-9]+)\.([0-9]+))");

    // The regular expression can be read like:
    // '('   : open a new sub-expression containing...
    // '['   : any character among...
    // '0-9' : the numbers 0 to 9
    // ']'
    // '+'   : with at least one such character
    // ')'   : end of the sub-expression
    // '\.'  : followed by a dot (has to be escaped with '\')
    // '('   : open a new sub-expression containing...
    //         ... exactly the same pattern as the first one
    // ')'   : end of the sub-expression

    // So we are looking for two sub-expressions, the first is the
    // integral part of the floating point number, and the second is the
    // fractional part.

    // It turns out that there are three locations in the input string
    // that match this pattern:
    sub(0,_); // {"6",  "25"}
    sub(1,_); // {"45", "6"}
    sub(2,_); // {"3",  "1415"}


regex_replace
-------------

.. code-block:: c++

    template<typename T>
    std::string regex_replace(std::string s, const std::string& reg, T fun); // [1]

This function will search in the string ``s`` using the regular expression ``reg`` (see ``regex_match`` above) to locate some expressions. Parts of these expressions can be *captured* by enclosing them in parenthesis. These captured sub-expressions are extracted from ``s``, stored inside a string vector, and fed to the user-supplied "replacement function" ``fun``. In turn, this function analyzes and/or modifies the captured sub-expressions to produce a new replacement string that will be inserted in place of the matched expression. The function is call for each match of the regular expression in ``s``.

If no sub-expression is captured, then the string vector that is fed to ``fun`` will be empty. Expressions are found and replaced in the order in which they appear in the input string ``s``.

This function is very similar to the ``sed`` program.

**Example:**

.. code-block:: c++

    std::string s = "a, b, c=5, d=9, e, f=2, g";

    // First a simple example.
    // We want to find all the "X=Y" expressions in this string and
    // add parentheses around "Y".
    std::string n = regex_replace(
        s,                    // the string to analyze
        R"(([a-z])=([0-9]))", // the regular expression
        // the replacement function
        [](vec1s m) {
            // "X" is in m[0], "Y" is in m[1].
            return m[0]+"=("+m[1]+")";
        }
    );

    // The regular expression can be read like:
    // '('   : open a new sub-expression containing...
    // '['   : any character among...
    // 'a-z' : the letters 'a' to 'z'
    // ']'
    // ')'   : end of the sub-expression
    // '='   : followed by the equal sign
    // '('   : open a new sub-expression containing...
    // '['   : any character among...
    // '0-9' : the numbers 0 to 9
    // ']'
    // ')'   : end of the sub-expression

    // The result is:
    n; // "a, b, c=(5), d=(9), e, f=(2), g"

    // A second, more complex example.
    // We take the same example as above, but this time we want to
    // change "X" to upper case and increment "Y" by one.
    std::string n = regex_replace(
        s,                    // the string to analyze
        R"(([a-z])=([0-9]))", // the regular expression
        // the replacement function
        [](vec1s m) {
            // Again, "X" is in m[0], "Y" is in m[1].

            // Read the integer "Y" and increment it
            uint_t k;
            from_string(m[1], k);
            ++k;

            // Build the replacement string
            return to_upper(m[0])+"="+to_string(k);
        }
    );

    // The result is:
    n; // "a, b, C=6, D=10, e, F=3, g"
