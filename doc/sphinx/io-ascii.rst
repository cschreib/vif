.. _ASCII tables:

ASCII tables
============

Defined in header ``<phypp/io/ascii.hpp>``.

ascii::read_table
-----------------

.. code-block:: c++

    template<typename ... Args>
    void ascii::read_table(std::string f, uint_t s, Args&& ... args); // [1]

    template<typename ... Args>
    void ascii::read_table(std::string f, auto_find_skip_tag t, Args&& ... args); // [2]

These functions read the content of the ASCII table whose path is given in ``f``, and stores the data inside the vectors listed in ``args``. For [1], the second argument ``s`` gives the number of lines to skip before starting reading data, and which can be used to ignore the "header" of the table. If the length of this header is unknown at the time of compilation, one can use function [2] instead and provide a ``ascii::auto_skip(c)`` as second argument, where ``c`` is the leading character of header lines (defaults to ``'#'`` if not provided).

The vectors that will receive the table data must be listed in ``args``, one after the other. Each column of the file will be stored in a separate vector, in the order in which they are provided to the function. If there is more columns than vectors, the extra columns are ignored. If there is more vectors than columns, the program will stop and report an error.

**Example:**

.. code-block:: none

    # my_table.dat
    # id    x     y
       0   10    20
       5   -1   3.5
       6    0    20
       8    5     1
      22  6.5    -5

We can read this in C++ with the following code:

.. code-block:: c++

    // Declare the vectors that will receive the data
    vec1u id;
    vec1f x, y;

    // Read the data.
    // The two lines below are equivalent in this case, but [2] is more flexible
    // if we ever change the format of the file by adding more lines to the header.
    ascii::read_table("my_table.dat", 2,                  id, x, y); // [1]
    ascii::read_table("my_table.dat", ascii::auto_skip(), id, x, y); // [2]

    // Use the vectors
    print(id); // {0, 5, 6, 8, 22}


.. warning:: Beware that, with these functions, the *names* of the C++ vectors are not used to identify the data in the file, and the information contained in the header is plainly ignored. Only the *order* of the columns and vectors matters.


**Spacing, alignment.** Columns in the file can be separated by any number of spaces and tabulations. The data does not need to be perfectly aligned to be read by this function, even though it is recommended for better human readability. The table may also contain empty lines, they will simply be ignored. However, "holes" in the table are not supported (i.e., rows that only have data for some columns, but not all). In such cases, the "hole" would be considered as white space and ignored, and the data from this column would actually be read from the next one. This would eventually trigger an error when trying to read the last columns, because there won't be enough data on this line.


**Data type.** Values in ASCII tables are not explicitly typed, so a column containing integers can be read as a vector of integers, floats, or even strings. As long as the data in the table can be converted to a value of the corresponding C++ vector using ``from_string()`` (see :ref:`String conversions`), this function will be able to read it. There are specific rules though:

* Strings may not contain spaces (since they would be understood as column separators). Adding quotes ``"..."`` will *not* change that. If you need to read strings containing spaces, you will have to read the table data manually with an ``std::ifstream``.
* For all numeric columns, if the value to be read is too large to fit in the corresponding C++ variable, the program will stop report the error. This will happen for example when trying to read a number like ``1e128`` inside a ``float`` vector. In such cases, use a larger data type to fix this (e.g., ``double`` in this particular case).


**Skipping columns.** If you want to ignore a specific column, you can use the "placeholder" symbol ``_`` instead of providing an actual vector. The corresponding data in the table will not be read. If you want to ignore ``n`` columns, you can use ``ascii::columns(n,_)``. With the example table above:

.. code-block:: c++

    // Declare the vectors that will receive the data
    vec1u id;
    vec1f x, y;

    // Read the data, ignoring the 'x' column
    ascii::read_table("my_table.dat", 2, id, _, y);

    // Read the data, ignoring the 'id' and 'x' columns
    ascii::read_table("my_table.dat", 2, _, _, y);
    ascii::read_table("my_table.dat", 2, ascii::columns(2,_), y);


**Reading 2D columns.** With the interface describe above, if you need to read ``N`` columns, you need to list ``N`` vectors when calling the function. This can be cumbersome for tables with a large number of columns. In the cases where it makes sense, you can choose to combine ``n`` adjacent columns of the ASCII table into a single 2D vector. The first dimension (``v.dims[0]``) will be the number of rows, and the second dimension (``v.dims[1]``) will be the number of columns (``n``). This can be done by specifying ``ascii::columns(n,v)``. With the example table above:

.. code-block:: c++

    // Declare the vectors that will receive the data
    vec1u id;
    vec2f v;

    // Read the data, merging the 'x' and 'y' columns into the 2D vector 'v'
    ascii::read_table("my_table.dat", 1, id, ascii::columns(2,v));

    // The data in 'v' is stored such that
    vec1f x = v(_,0);
    vec1f y = v(_,1);


**Multiple 2D columns.** The trick of reading 2D columns can be extended to read several columns into multiple 2D vectors by following a pattern. A typical case is when you have, say, three quantities 'A', 'B', and 'C' listed in the table, each with their values and uncertainties:

.. code-block:: none

    # abc_data.dat
    # id    A  Aerr     B  Berr     C  Cerr
       0   10   1.0     1   0.1    -1     1
       5   -1   3.5     2   0.2     1     1
       6    0     6     3   0.2     1     2
       ...

This table can be read easily into two 2D vectors ``value`` and ``uncertainty`` by using ``ascii::columns(3,value,uncertainty)``. This is interpreted as "read 3 sets of columns, each containing value and uncertainty":

.. code-block:: c++

    // Declare the vectors that will receive the data
    vec1u id;
    vec2f value, uncertainty;

    // Read the data
    ascii::read_table("abc_data.dat", 2, id, ascii::columns(3,value,uncertainty));

    // The data is stored in 'value' and 'uncertainty' such that
    vec1f a    = value(_,0);
    vec1f aerr = uncertainty(_,0);
    vec1f b    = value(_,1);
    vec1f berr = uncertainty(_,1);
    vec1f c    = value(_,2);
    vec1f cerr = uncertainty(_,2);

This can also be mixed with the placeholder symbol `_` to skip column (see above):

.. code-block:: c++

    // If we only wanted to read the values, and not the uncertainties, we could write:
    ascii::read_table("abc_data.dat", 2, id, ascii::columns(3,value,_));


ascii::write_table, ascii::write_table_csv, ascii::write_table_hdr
------------------------------------------------------------------

.. code-block:: c++

    void ascii::write_table(std::string f, uint_t w, const Args& ... args); // [1]

    void ascii::write_table_csv(std::string f, const Args& ... args); // [2]

    void ascii::write_table_hdr(std::string f, uint_t w, const vec1s& hdr, const Args& ... args); // [3]

    void ascii::write_table_hdr(std::string f, uint_t w, ftable(...)); // [4]

These functions will write the data of the vectors listed in ``args`` into the file whose path is given in ``f``. The data will be formated in a human-readable form, commonly called "ASCII" format. With function [1], the data is written in separate columns of fixed width of ``w`` characters, and white spaces are used to fill the empty space between columns. With function [2], the columns will be separated by commas (``','``), as "CSV" stands for "comma-separated values". In all cases, all columns must have the same number of rows.

**Example:**

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec1i x = {125,568,9852,12,-51};
    vec1i y = {-56,157,2,99,1024};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", 8, id, x, y);
    // ... or a CSV file
    ascii::write_table_csv("my_table.csv", id, x, y);

The content of ``my_table.dat`` will be:

.. code-block:: none

           1     125     -56
           2     568     157
           3    9852       2
           4      12      99
           5     -51    1024

The content of ``my_table.csv`` will be:

.. code-block:: none

    1,125,-56
    2,568,157
    3,9852,2
    4,12,99
    5,-51,1024


.. note:: Human-readable formats are simple, and quite convenient for small files. But if the volume of data is huge, consider instead using :ref:`FITS files` instead. This will be faster to read and write, and will also occupy less space on your disk.

**Output format.** When providing vectors of floats or doubles, these functions will convert the values to strings using the default C++ format. See discussion in :ref:`String conversions`. When this is not appropriate, you can use the ``format::...`` functions to modify the output format, as you would with ``to_string()``.

**Writing 2D vectors.** These functions support writing 2D vectors as well. They are interpreted as containing multiple columns: the number of columns is the first dimension of the vector (``v.dims[0]``), and the number of rows is the second dimension (``v.dims[1]``). For them to be recognized by the function, you must wrap them in ``ascii::columns(n,v)``, where ``n`` is the number of columns. For example:

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec2f value(5,2);
    value(_,0) = {0.0, 1.2, 5.6, 9.5, 1.5};
    value(_,1) = {9.6, 0.0, 4.5, 0.0, 0.0};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", 8, id, ascii::columns(2,value));

The content of ``my_table.dat`` will be:

.. code-block:: none

           1       0     9.6
           2     1.2       0
           3     5.6     4.5
           4     9.5       0
           5     1.5       0


**Multiple 2D vectors.** As for ``ascii::read_table()``, you can use the above mechanism to write multiple 2D columns following a pattern by listing them in the ``ascii::columns()``. For example, ``ascii::columns(n, value, uncertainty)`` will write ``n`` pairs of columns, with ``value`` and ``uncertainty`` in each pair.

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec2f value(5,2);
    value(_,0) = {0.0, 1.2, 5.6, 9.5, 1.5};
    value(_,1) = {9.6, 0.0, 4.5, 0.0, 0.0};
    vec2f uncertainty(5,2);
    uncertainty(_,0) = {0.01, 0.03, 0.05, 0.09, 0.01};
    uncertainty(_,1) = {0.02, 0.05, 0.01, 0.21, 0.04};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", 8, id, ascii::columns(2,value,uncertainty));

The content of ``my_table.dat`` will be:

.. code-block:: none

           1       0    0.01     9.6    0.02
           2     1.2    0.03       0    0.05
           3     5.6    0.05     4.5    0.01
           4     9.5    0.09       0    0.21
           5     1.5    0.01       0    0.04


**Headers.** Functions [3] and [4] will adopt the same output format as function [1]. The only difference is that these functions all also write a "header" at the beginning of the file. This header is composed of two lines each starting with a hash ``'#'``. The first line lists the names of the columns, automatically aligned with the data, while the second line is empty. In C++, the header must be specified as a 1D vector of strings (``hdr``), each containing the name of a column. Therefore there must be as many elements in this vector as there are columns to write.

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec1i x = {125,568,9852,12,-51};
    vec1i y = {-56,157,2,99,1024};

    // Write these in a simple ASCII file
    ascii::write_table_hdr("my_table.dat", 8, {"id", "x", "y"}, id, x, y); // [3]

The content of ``my_table.dat`` will be:

.. code-block:: none

    #     id       x       y
    #
           1     125     -56
           2     568     157
           3    9852       2
           4      12      99
           5     -51    1024

If, as in the example above, the vectors in C++ have meaningful names that could be used directly as names for the columns in the file, you can use function [4] with the ``ftable()`` macro to automatically build the header. With this function, the example above becomes:

.. code-block:: c++

    // Write these in a simple ASCII file
    ascii::write_table_hdr("my_table.dat", 8, ftable(id, x, y)); // [4]

This also works for 2D vectors. In such cases, ``_i`` is appended to the name of the vector for each column ``i``. If you need better looking headers, you can always write them manually using function [3].
