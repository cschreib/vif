.. _ASCII tables:

ASCII tables
============

Defined in header ``<vif/io/ascii.hpp>``.

ascii::read_table
-----------------

.. code-block:: c++

    template<typename ... Args>
    void ascii::read_table(std::string f, Args&& ... args); // [1]

    template<typename ... Args>
    void ascii::read_table(std::string f, const ascii::input_format& o, Args&& ... args); // [2]

These functions read the content of the ASCII table whose path is given in ``f``, and stores the data inside the vectors listed in ``args``. Each column of the file will be stored in a separate vector, in the order in which they are provided to the function. If there is more columns than vectors, the extra columns are ignored. If there is more vectors than columns, the program will stop and report an error.

Function [1] assumes a number of default options regarding the layout of the table. In particular, it assumes the columns are separated by white spaces, and that there may be a header at the beginning of the table (lines starting with the character ``'#'``) that must be skipped before reading the data. See below for more detail on how the table data is read. Below is an example of such a table.

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

    // Read the data, asking for three columns.
    ascii::read_table("my_table.dat", id, x, y);

    // Use the vectors
    print(id); // {0, 5, 6, 8, 22}


.. warning:: Beware that, with these functions, the *names* of the C++ vectors are not used to identify the data in the file, and the information contained in the table header is plainly ignored. Only the *order* of the columns and vectors matters.

Function [2] allows you to fine tune how the table data is read using the option structure ``o``. This structure has the following members:

.. code-block:: c++

    struct input_format {
        bool        auto_skip    = true;
        std::string skip_pattern = "#";
        uint_t      skip_first   = 0;
        std::string delim        = " \t";
        bool        delim_single = false;
    };

* ``auto_skip`` and ``skip_pattern``. When ``auto_skip`` is set to ``true``, the function will automatically ignore all the lines starting with ``skip_pattern`` (typically, the header).
* ``skip_first``. This is an alternative way to skip a header, when the header has always the same number of lines (one or two, typically), but when the lines do not start with a specific character. By setting this option to a positive number, the function will skip the first ``skip_first`` lines before reading the data.
* ``delim`` and ``delim_single``. The string ``delim`` determines what characters are used to separate the columns in the file. When ``delim_single`` is ``false``, ``delim`` is interpreted as a list of characters that can be expected in between columns, in any number and order. For example, ``delim = " \t"; delim_single = false;`` states that columns can be separated by any number of white spaces and tabulations. On the other hand, when ``delim_single`` is ``true``, ``delim`` is interpreted as a fixed string that must be found between each column, and any other character is considered part of the column data itself. For example, ``delim = ","; delim_single = true;`` would specify a comma-separated table.

Some pre-defined sets of options are made available for simplicity:

* ``ascii::input_format::standard()``. This is the default behavior of ``read_table()``, when no options are given (see default values above). With this setup, columns in the file can be separated by any number of spaces and tabulations. The data does not need to be perfectly aligned to be read correctly, even though it is recommended for better human readability. The table may also contain empty lines, they will simply be ignored.

  However, "holes" in the table are not supported (i.e., rows that only have data for some columns, but not all). For example:

  .. code-block:: none

      # my_table.dat
      # id    x     y
         0   10    20
         5   -1   3.5
         6         20  # <-- not OK!
         8    5     1
        22  6.5    -5

  In such cases (see how ``x`` is missing a value on the third row), the "hole" would be considered as whitespace and ignored, and the data from this column would actually be read from the next one (so ``x`` would have a value of ``20`` for this row). This would eventually trigger an error when trying to read the last columns, because there won't be enough data on this line (there is no value for ``y``). Therefore, *every* row must have a value for *every* column. If data is missing, use special values such as ``-1`` or ``NaN`` to indicate it.

  This also means that string columns cannot contain spaces in them, since they would otherwise be understood as column separators. Adding quotes ``"..."`` will *not* change that. If you need to read strings containing spaces, you should use another table format (such as CSV, see below).

  Using this format is done as follows:

  .. code-block:: c++

      // Declare the vectors that will receive the data
      vec1u id;
      vec1f x, y;

      // Read the data, asking for three columns.
      ascii::read_table("my_table.dat", ascii::input_format::standard(), id, x, y);
      // This is equivalent to
      ascii::read_table("my_table.dat", id, x, y);

  You can also use it as a starting point to create customized options:

  .. code-block:: c++

      ascii::input_format opts = ascii::input_format::standard();
      opts.skip_pattern = "%"; // skip lines starting with '%'
      ascii::read_table("my_table.dat", opts, id, x, y);

* ``ascii::input_format::csv()``. This preset enables loading comma-separated values (CSV) tables. In these tables, columns are separated by a single comma (``','``). Contrary to the ``standard`` format, spaces are considered to be a significant part of the data, and will not be trimmed.


The information below applies to any type of table.

**Data type.** Values in ASCII tables are not explicitly typed, so a column containing integers can be read as a vector of integers, floats, or even strings. As long as the data in the table can be converted to a value of the corresponding C++ vector using ``from_string()`` (see :ref:`String conversions`), this function will be able to read it. Note that, for all numeric columns, if the value to be read is too large to fit in the corresponding C++ variable, the program will stop and report an error. This will happen for example when trying to read a number like ``1e128`` inside a ``float`` vector. In such cases, use a larger data type to fix this (e.g., ``double`` in this particular case).


**Skipping columns.** If you want to ignore a specific column, you can use the "placeholder" symbol ``_`` instead of providing an actual vector. The corresponding data in the table will not be read. If you want to ignore ``n`` columns, you can use ``ascii::columns(n,_)``. With the example table above:

.. code-block:: c++

    // Declare the vectors that will receive the data
    vec1u id;
    vec1f x, y;

    // Read the data, ignoring the 'x' column
    ascii::read_table("my_table.dat", 2, id, _, y);

    // Read the data, ignoring the 'id' and 'x' columns.
    // This can be done with two '_':
    ascii::read_table("my_table.dat", 2, _, _, y);
    // ... or with the function ascii::columns():
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


ascii::write_table
------------------

.. code-block:: c++

    void ascii::write_table(std::string f, const Args& ... args); // [1]

    void ascii::write_table(std::string f, const ascii::output_format& o, const Args& ... args); // [2]

    void ascii::write_table(std::string f, ftable(...)); // [3]

    void ascii::write_table(std::string f, const ascii::output_format& o, ftable(...)); // [4]

These functions will write the data of the vectors listed in ``args`` into the file whose path is given in ``f``. The data will be formated in a human-readable form, colloquially called "ASCII" format. In all cases, all columns must have the same number of rows, otherwise the function will report an error.

Function [1] uses a "standard" format, where the data is written in separate columns, separated and automatically aligned by white spaces. See below for more detail on how the table data is written. Here is a simple example.

**Example:**

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec1i x = {125,568,9852,12,-51};
    vec1i y = {-56,157,2,99,1024};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", id, x, y);

The content of ``my_table.dat`` will be:

.. code-block:: none

    1  125  -56
    2  568  157
    3 9852    2
    4   12   99
    5  -51 1024


.. note:: Human-readable formats are simple, and quite convenient for small files. But if the volume of data is huge, consider instead using :ref:`FITS files` instead. This will be faster to read and write, and will also occupy less space on your disk.

Function [2] allows you to change the output format by specifying a number of options in the option structure ``o``. This structure has the following members:

.. code-block:: c++

    struct output_format {
        bool        auto_width   = true;
        uint_t      min_width    = 0;
        std::string delim        = " ";
        std::string header_chars = "# ";
        vec1s       header;
    };

* ``auto_width``. When set to ``true`` (the default), the function will compute the maximum width (in characters) of each column before writing the data to the disk. It will then use this maximum width to nicely align the data in each column. Note that it also takes into account the width of the header string (see below). This two-step process reduces performances a bit, and for large data sets you may want to disable it by setting this option to ``false``. In this case, either the data is written without alignment (still readable by a machine, but not really by a human), or with a fixed common width if ``min_width`` is set to a positive value.
* ``min_width``. This defines the minimum width allowed for a column, in characters. The default is zero, which means columns can be as narrow as one single character if that is all the space they require.
* ``delim``. This string defines which character(s) should be used to separate columns in the file. The default is to use a single white space (plus any alignment coming from adjusting the column widths).
* ``header`` and ``header_chars``. These variables can be used to print a header at the beginning of the file, before the data. This header can be used by a human (or, possibly, a machine) to understand what kind of data is contained in the table. The header will be written on a single line, starting with ``header_chars`` (the header starting string). Then, each column written in the file must have its name listed in the ``header`` array, in the same order as given in ``args``.

Some pre-defined sets of options are made available for simplicity:

* ``ascii::output_format::standard()``. This is the default behavior of ``write_table()``, when no options are given (see default values above). With this setup, columns in the file are separated by at least one white space character (and possibly more, for alignment). The data in a given column is automatically aligned.

  Using this format is done as follows:

  .. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec1i x = {125,568,9852,12,-51};
    vec1i y = {-56,157,2,99,1024};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", output_options::standard(), id, x, y);
    // This is equivalent to
    ascii::write_table("my_table.dat", id, x, y);

  You can also use it as a starting point to create customized options:

  .. code-block:: c++

      ascii::output_format opts = ascii::output_format::standard();
      opts.header_chars = "% "; // begin header line with '% ' instead of '# '
      ascii::write_table("my_table.dat", opts, id, x, y);

* ``ascii::output_format::csv()``. This preset enables writing comma-separated values (CSV) tables. In these tables, columns are separated by a single comma (``','``), and the data is not aligned at all.

The information below applies to any type of table.

**Output format.** When providing vectors of floats or doubles, these functions will convert the values to strings using the default C++ format. See discussion in :ref:`String conversions`. When this is not appropriate, you can use the ``format::...`` functions to modify the output format, as you would with ``to_string()``.

**Writing 2D vectors.** These functions support writing 2D vectors as well. They are interpreted as containing multiple columns: the number of columns is the first dimension of the vector (``v.dims[0]``), and the number of rows is the second dimension (``v.dims[1]``). For them to be recognized by the function, you must wrap them in ``ascii::columns(n,v)``, where ``n`` is the number of columns. For example:

.. code-block:: c++

    // We have some data
    vec1u id = {1,2,3,4,5};
    vec2f value(5,2);
    value(_,0) = {0.0, 1.2, 5.6, 9.5, 1.5};
    value(_,1) = {9.6, 0.0, 4.5, 0.0, 0.0};

    // Write these in a simple ASCII file
    ascii::write_table("my_table.dat", id, ascii::columns(2,value));

The content of ``my_table.dat`` will be:

.. code-block:: none

    1   0 9.6
    2 1.2   0
    3 5.6 4.5
    4 9.5   0
    5 1.5   0


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
    ascii::write_table("my_table.dat", id, ascii::columns(2,value,uncertainty));

The content of ``my_table.dat`` will be:

.. code-block:: none

    1   0 0.01 9.6 0.02
    2 1.2 0.03   0 0.05
    3 5.6 0.05 4.5 0.01
    4 9.5 0.09   0 0.21
    5 1.5 0.01   0 0.04


**Easy headers.** Functions [3] and [4] will adopt the same output format as functions [1] and [2]. The only difference is that they will automatically create the header based on the names of the C++ variables that are listed in argument. To do so, you must use the ``ftable()`` macro to list the data to be written. For example:

.. code-block:: c++

    // Write these in a simple ASCII file
    ascii::write_table_hdr("my_table.dat", ftable(id, x, y)); // [4]

This also works for 2D vectors. In such cases, ``_i`` is appended to the name of the vector for each column ``i``. If you need better looking headers, you can always write them manually using function [2].
