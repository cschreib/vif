File system
===========

Defined in header ``<vif/io/filesystem.hpp>``.

file::exists
------------

.. code-block:: c++

    bool file::exists(const std::string& f); // [1]

    template<std::size_t D>
    vec<D,bool> file::exists(const vec<D,std::string>& f); // [2]

The function [1] returns ``true`` if a file (or directory) exists at the location given in ``f``, and ``false`` otherwise.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    bool b;
    b = file::exists("~/.vifrc");      // hopefully true
    b = file::exists("/i/do/not/exist"); // probably false


file::is_older
--------------

.. code-block:: c++

    bool file::is_older(const std::string& f1, const std::string& f2); // [1]

    template<std::size_t D>
    vec<D,bool> file::is_older(const vec<D,std::string>& f1, const std::string& f2); // [2]

The function [1] returns ``true`` if the file (or directory) ``f1`` is *older* than the file (or directory) ``f2``. The "age" of a file corresponds to the time spent since that file was last modified. If one of the two files does not exists, the function returns ``false``.

The function [2] is the vectorized version of [1], where all files in ``f1`` are compared against the same file ``f2``.

**Example:**

.. code-block:: c++

    bool b = file::is_older("~/.vifrc", "/usr/bin/cp"); // maybe false?


file::list_directories
----------------------

.. code-block:: c++

    vec1s file::list_directories(const std::string& d, const std::string& p = "");

This function returns the list of all the subdirectories inside the directory ``d`` that match the search pattern ``p``. Only the names of the directories are returned, not their full paths. When the search pattern is empty (default), all the directories are returned. Otherwise, the pattern can contain any number of "wildcard" character ``*`` to filter the output, as when listing files in bash (see examples below).

Hidden directories are ignored. An empty list is returned if there is no subdirectory matching the pattern, or if the operation could not be performed (that is, if ``d`` doest not exist, or is not a directory, or if you do not have read access to it). The function does not look inside subdirectories recursively. The *order* of the directories in the output list is undefined and should not be relied upon: if you need a sorted list, you have to sort it yourself.

**Example:**

.. code-block:: c++

    vec1s d = file::list_directories("./");
    d; // subdirectories of the working directory

    d = file::list_directories("/path/to/vif/");
    d; // {"cmake", "test", "doc", "bin", "include", "tools"}

    d = file::list_directories("/path/to/vif/", "t*");
    d; // {"test", "tools"}


file::list_files
----------------

.. code-block:: c++

    vec1s file::list_files(const std::string& d, const std::string& p = "");

This function returns the list of all the files inside the directory ``d`` that match the search pattern ``p``. Only the names of the files are returned, not their full paths. When the search pattern is empty (default), all the files are returned. Otherwise, the pattern can contain any number of "wildcard" character ``*`` to filter the output, as when listing files in bash (see examples below).

Hidden files are ignored. An empty list is returned if there is no file matching the pattern, or if the operation could not be performed (that is, if ``d`` doest not exist, or is not a directory, or if you do not have read access to it). The function does not look inside subdirectories recursively. The *order* of the files in the output list is undefined and should not be relied upon: if you need a sorted list, you have to sort it yourself.

**Example:**

.. code-block:: c++

    vec1s d = file::list_files("./");
    d; // files in the working directory

    d = file::list_files("/path/to/vif/doc");
    d; // {"vif.pdf", "compile.sh", "vif.tex"}

    d = file::list_files("/path/to/vif/doc", "*.tex");
    d; // {"vif.tex"}


file::explorer
--------------

.. code-block:: c++

    class file::explorer {
    public:
        struct file_data {
            std::string full_path;
            std::string name;
            uint_t size;
            bool is_hidden = false;
            bool is_dir = false;
        };

        // Constructors
        explorer(); // [1]
        explorer(const std::string& d, const std::string& p = ""); // [2]

        void open(const std::string& d, const std::string& p = ""); // [3]
        bool find_next(file_data& f); // [4]
        void close(); // [5]
    };

This class allows you to browse through the content of a directory, to list the files and other directories it contains. Its interface is similar to ``std::ifstream``: it can be default-constructed ([1]) then initialized with ``open()`` ([3]), or this can be achieved in a single step using the constructor [2], which takes the same arguments as ``open()``.

To use this class, you must first open a directory, either with [2] or [3]: the class will attempt to open the directory ``d`` and initialize a new search, optionally with a search pattern ``p``. If the directory does not exist or is not readable, ``open()`` will return ``false``, and the search will be aborted (subsequent calls to ``find_next()`` will return ``false``). The search pattern must constain at least one wildcard character ``*`` to indicate which part of the files (or directories) name is allowed to vary, like when listing files in bash.

Once the directory is open, you can iterate over its content using ``find_next()`` ([4]). This function take a pre-constructed ``file_data`` in argument, in which it will fill the details of the next file it found. If no more file is found (i.e., if the previous call to ``find_next()`` returned the last file), this function returns ``false`` and the ``file_data`` is not modified (should not be used).

The ``file_data`` object is a simple structure holding basic informations about the file (or directory): ``name`` is the name the file, ``full_path`` is the name appended to the search directory ``d``, ``size`` is the size of the file (in bytes), ``is_hidden`` is ``true`` for hidden files or directories, and ``is_dir`` is ``true`` for directories and ``false`` for files.

Once you are done with a search, you can let the ``explorer`` instance be destroyed at the end of its scope. This will call ``close()`` ([5]) automatically. If you need to close the access to the directory immediately, or if you wish to start another search, you can also call ``close()`` explicitly.

**Example:**

.. code-block:: c++

    // Create explorer
    file::explorer e;

    // Try to open directory
    if (e.open("some/dir", "*.cpp")) {
        // Success, now list the files/directories
        file::explorer::file_data f;
        while (e.find_next(f)) {
            // We found a file/directory, do something with it:
            print("found ",
                (f.is_hidden ? "hidden " : ""),
                (f.is_dir ?    "directory " : "file "),
                f.name, " (size: ", f.size, ")");
        }
    } else {
        // Failed
        error("could not open directory some/dir");
    }


file::mkdir
-----------

.. code-block:: c++

    bool file::mkdir(const std::string& d); // [1]

    template<std::size_t D>
    vec<D,bool> file::mkdir(const vec<D,std::string>& d); // [2]

The function [1] creates a new directory at the path given in argument (including all the parent directories, if necessary), and returns ``true``. If the directory could not be created (e.g., because of permission issues), the function returns ``false``. If the directory already exists, the function does nothing and returns ``true``. This function is equivalent to the bash function ``mkdir -p``.

The function [2] is the vectorized version of [1].

**Example:**

.. code-block:: c++

    bool b = file::mkdir("/path/to/vif/a/new/directory");
    // Will most likely create the directories:
    //  - /path/to/vif/a
    //  - /path/to/vif/a/new
    //  - /path/to/vif/a/new/directory
    b; // maybe true or false, depending on your permissions


file::copy
----------

.. code-block:: c++

    bool file::copy(const std::string& from, const std::string& to);

This function creates a copy of the file ``from`` at the location given in ``to`` and returns ``true``. If the new file could not be created (e.g., because of permission issues or because its parent directory does not exist), or if the file ``from`` could not be found or read, the function returns ``false``. If the file ``to`` already exists, it will be overwritten without warning. Copying directories is not presently supported. This function is equivalent to the bash function ``cp -f``.

**Example:**

.. code-block:: c++

    bool b = file::copy("/home/joe/.vifrc", "/home/bob/.vifrc");
    b; // maybe true or false, depending on your permissions


file::remove
------------

.. code-block:: c++

    bool file::remove(const std::string& f);

This function will delete the file (or directory) given in argument and return ``true`` on success, or if the file did not exist. It will return ``false`` only if the file exists but could not be removed (i.e., because you are lacking the right permissions). This function is equivalent to the bash function ``rm -rf``.

**Example:**

.. code-block:: c++

    // That's a bad idea, but for the sake of the example...
    bool b = file::remove("/home/joe/.vifrc");
    b; // probably true


file::to_string
---------------

.. code-block:: c++

    std::string file::to_string(const std::string& f);

This function reads the content of the file whose path is given in argument, stores all the characters (including line break characters and spaces) inside a string and returns it. If the file does not exist, the function returns an empty string.

.. warning::

    This is a very sub-optimal way of reading the content of a file, and it should only be attempted on short files. If you need to read a file line by line, use ``std::getline()`` instead. If you need to read a data table, use the dedicated functions in :ref:`ASCII tables`.

**Example:**

.. code-block:: c++

    std::string r = file::to_string("/etc/lsb-release");
    // 'r' now contains all the lines of the file, each
    // separated by a newline character '\n'.

