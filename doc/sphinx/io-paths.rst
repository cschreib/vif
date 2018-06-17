File paths manipulation
=======================

Defined in header ``<phypp/io/filesystem.hpp>``.

file::directorize
-----------------

.. code-block:: c++

    std::string file::directorize(const std::string& p); // [1]

    template<std::size_t D>
    vec<D,std::string> file::directorize(const vec<D,std::string>& p); // [2]

The function [1] modifies the path given in argument to make sure that a file name can be appended to it and form a valid file path. In UNIX systems, for example, the function ensures that the path ends with a forward slash ``/``.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    std::string p;
    p = file::directorize("/some/path");     // "/some/path/"
    p = file::directorize("/another/path/"); // "/another/path/"


file::is_absolute_path
----------------------

.. code-block:: c++

    std::string file::is_absolute_path(const std::string& p); // [1]

    template<std::size_t D>
    vec<,D,std::string> file::is_absolute_path(const vec<,D,std::string>& p); // [2]

The function [1] returns ``true`` if its argument describes an absolute file path, and ``false`` otherwise. In UNIX systems, for example, this is equivalent to checking that the path starts with a forward slash ``/`` (referencing the root directory).

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    bool b = file::is_absolute_path("/some/path");           // true
    b = file::is_absolute_path("../sub/directory/file.txt"); // false


file::get_basename
------------------

.. code-block:: c++

    std::string file::get_basename(const std::string& p); // [1]

    template<std::size_t D>
    vec<,D,std::string> file::get_basename(const vec<,D,std::string>& p); // [2]

The function [1] extracts the name of a file from its full path given in argument. If this path is that of a directory, the function returns the name of this directory. This behavior is similar to the bahs function ``basename``.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    std::string n = file::get_basename("/some/path");      // "path"
    n = file::get_basename("/another/path/to/a/file.txt"); // "file.txt"


file::get_extension
-------------------

.. code-block:: c++

    std::string file::get_extension(const std::string& f); // [1]

    template<std::size_t D>
    vec<,D,std::string> file::get_extension(const vec<,D,std::string>& f); // [2]

The function [1] scans the provided string to look for a file extension. The "extension" is whatever is found at the end the string after the *last* dot (and including this dot), for example ``".cpp"``. If an extension is found, this function returns it (including the leading dot), else it returns an empty string.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    vec1s v = {"p1_m2.txt", "p3_c4.fits", "p1_t8.dat.fits", "readme"};
    vec1s s = file::get_extension(v); // {".txt", ".fits", ".fits", ""}


file::remove_extension
----------------------

.. code-block:: c++

    std::string file::remove_extension(const std::string& f); // [1]

    template<std::size_t D>
    vec<,D,std::string> file::remove_extension(const vec<,D,std::string>& f); // [2]

The function [1] scans the provided string to look for a file extension. The "extension" is whatever is found at the end the string after the *last* dot (and including this dot), for example ``".cpp"``. If an extension is found, this function returns the input string with this extension removed. If no extension is found, the input string returned unchanged.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    vec1s v = {"p1_m2.txt", "p3_c4.fits", "p1_t8.dat.fits", "readme"};
    vec1s s = file::remove_extension(v); // {"p1_m2", "p3_c4", "p1_t8.dat", "readme"}


file::split_extension
---------------------

.. code-block:: c++

    std::pair<std::string> file::split_extension(const std::string& f); // [1]

    template<std::size_t D>
    vec<D,std::pair<std::string>> file::split_extension(const vec<,D,std::string>& f); // [2]

The function [1] scans the provided string to look for a file extension. The "extension" is whatever is found at the end the string after the *last* dot (and including this dot), for example ``".cpp"``. If an extension is found, this function splits the input string into two substrings, the first being the string with the extension removed (see ``file::remove_extension()``), and the second being the extension itself (see ``file::get_extension()``).

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    auto p = file::split_extension("p1_m2.txt");
    p.first; // "p1_m2"
    p.second; // ".txt"


file::get_directory
-------------------

.. code-block:: c++

    string file::get_directory(const std::string& p); // [1]

    template<std::size_t D>
    vec<,D,std::string> file::get_directory(const vec<,D,std::string>& p); // [2]

The function [1] scans the path given in argument and returns the path to the parent directory. This behavior is similar to the bash function ``dirname``, except that here the returned path always ends with a forward slash ``/``.

**Example:**

.. code-block:: c++

    std::string n;
    n = file::get_directory("/some/path");                  // "/some/"
    n = file::get_directory("/another/path/to/a/file.txt"); // "/another/path/to/a/"
