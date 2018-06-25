.. _String conversions:

String conversions
==================

Defined in header ``<phypp/core/string_conversion.hpp>``.

to_string, to_string_vector
---------------------------

.. code-block:: c++

    template<typename Type>
    std::string to_string(const Type& v); // [1]

    template<std::size_t Dim, typename Type>
    vec<Dim,std::string> to_string_vector(const vec<Dim,Type>& v); // [2]

The function [1] will convert the value ``v`` into a string. This value can be of any type, as long as it is convertible to string. In particular, ``v`` can be a vector, in which case the output string will contain all the values of the vector, separated by commas, and enclosed inside curly braces: ``"{a,b,c,d,....}"``.

**Example:**

.. code-block:: c++

    to_string(2);            // "2"
    to_string(true);         // "1"
    to_string("foo");        // "foo"
    to_string(vec1i{2,5,9}); // "{2, 5, 9}"

The function [2] allows you to convert *each value* of a vector into a separate string, and store them in a string vector.

.. code-block:: c++

    to_string(vec1i{2,5,9});        // "{2, 5, 9}"
    to_string_vector(vec1i{2,5,9}); // {"2", "5", "9"}

Internally, the argument is converted to a string using the ``std::ostream`` ``operator<<``. This means that most types from the standard C++ or external C++ libraries will be convertible to a string out of the box. If you encounter some errors for a particular type, this probably means that the ``operator<<`` is missing and you have to write it yourself. Here is how you would do that:

.. code-block:: c++

    // We want to make this structure printable
    struct test {
        std::string name;
        int i, j;
    };

    // We just need to write this function
    std::ostream& operator<< (std::ostream& o, const test& t) {
        o << t.name << ":{i=" << i << " j=" << j << "}";
        return o; // do not forget to always return the std::ostream!
    }

    // The idea is always to rely on the existence of an operator<<
    // function for the types that are contained by your structure
    // or class. In our case, 'std::string' and 'int' are already
    // conertible to string, This is the standard C++ way of doing
    // string conversions, but it can be annoying to use regularly
    // because the "<<" are taking a lot of screen space.
    // 'to_string()' makes it much more convenient.

    // Now we can convert!
    test t = {"toto", 5, 12};
    std::string s = to_string(t);
    s; // "toto:{i=5 j=12}"


format::precision, format::scientific, format::fixed
----------------------------------------------------

.. code-block:: c++

    template<typename Type>
    /* ... */ format::scientific(const Type& v); // [1]

    template<typename Type>
    /* ... */ format::precision(const Type& v, uint_t ndigit); // [2]

    template<typename Type>
    /* ... */ format::fixed(const Type& v); // [3]

The ``to_string()`` and ``to_string_vector()`` functions adopt a default format for converting numbers into strings. While integers have a unique and natural string representation, floating point numbers often require a choice regarding the number of significant digits, and whether scientific notation should be used. By default, these functions follow the behavior of ``std::ostream``, which is to only use scientific notation when the number would be "too big" (or "too small").

Function [1], ``format::scientific()``, will specify that it's argument ``v`` *must* be formated using the scientific notation.

**Example:**

.. code-block:: c++

    double v = 0.15;
    to_string(v);                     // "0.15"
    to_string(format::scientific(v)); // "1.500000e-01"

Function [2], ``format::precision()``, will specify that it's first argument ``v`` *must* be formated using ``ndigit`` digits. Normally, "digits" include numbers on either side of the decimal separator, so ``"3.15"``, ``"31.5"``, and ``"315"`` are all three digits. When not in scientific format, trailing zeros after the decimal separator will still be removed, so the total number of digits may still be less than ``ndigit``.

**Example:**

.. code-block:: c++

    double v = 0.15;
    to_string(v);                       // "0.15"
    to_string(format::precision(v, 8)); // "0.15"

    v = 0.123456789123456789;
    to_string(v);                       // "0.123457"
    to_string(format::precision(v, 8)); // "0.12345679"

Function [3], ``format::fixed()``, will format the value with a fixed number of digits after the decimal separator. Trailing zeroes will not be removed. This is best used in combination with ``format::precision``, which then specifies how many digits to keep after the decimal separator (digits before the separator do not count).

**Example:**

.. code-block:: c++

    double v = 0.15;
    to_string(v);                                      // "0.15"
    to_string(format::fixed(format::precision(v, 8))); // "0.15000000"

    v = 3150.15;
    to_string(v);                                      // "3150.15"
    to_string(format::fixed(format::precision(v, 8))); // "3150.15000000"

    v = 0.123456789123456789;
    to_string(v);                                      // "0.123457"
    to_string(format::fixed(format::precision(v, 8))); // "0.12345679"

Note that all these functions can be used in other contexts than just ``to_string()`` and ``to_string_vector()``, essentially whenever a conversion to string is performed. See for example ``ascii::write_table()``.


from_string
-----------

.. code-block:: c++

    template<typename Type>
    bool from_string(const std::string& s, const Type& v); // [1]

    template<std::size_t D, typename Type>
    vec<D,bool> from_string(const vec<D,std::string>& v, vec<D,Type>& v); // [2]

The function [1] tries to convert the string ``s`` into a C++ value ``v`` and returns ``true`` in case of success. If the string cannot be converted into this value, for example if the string contains letters and the value has an arithmetic type, or if the number inside the string is too big to fit inside the C++ value, the function will return ``false``. In this case, the value of ``v`` is undefined.

The version [2] will try to convert each value inside the string vector ``s``, and will store the converted values inside the vector ``v``. It will automatically take care or resizing the vector ``v``, so you can pass an empty vector in input. The return value is an array of boolean values, corresponding to the success or failure of conversion for each individual value inside ``s``. If an element of ``s`` failed to convert, the corresponding value in ``v`` will be undefined.

**Example:**

.. code-block:: c++

    float f;
    bool b = from_string("3.1415", f);
    b; // true
    f; // 3.1415

    b = from_string("abcdef", f);
    b; // false;
    f; // ??? could be 3.1415, or NaN, or anything else

    vec1f fs;
    vec1b bs = from_string({"1", "1.00e5", "abc", "1e128", "2.5"}, fs);
    bs; // {true, true, false, false, true}
    fs; // {1,    1e5,  ???,   ???,   2.5}
