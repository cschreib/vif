Interacting with the operating system
=====================================

Defined in header ``<phypp/utility/os.hpp>``.

system_var
----------


.. code-block:: c++

    std::string system_var(const std::string& v, const std::string& d); // [1]

    template<typename T>
    T system_var(string v, T d); // [2]

These functions looks inside the operating system environment for a variable named ``v`` (using the C function ``getenv()``). If this variable exists, its value is returned as a string ([1]), or converted into a value of type ``T`` ([2]). If the conversion fails, or if the environment variable does not exist, the default value ``d`` is returned instead.

Environment variables are complementary to :ref:`Command line arguments`. They are mostly used to store constant, system-specific configurations that will typically change only between one machine or one user to another. Because they seldom change, it would be tedious to have to specify these configurations as command line arguments and provide them for *each* call of a given program. Instead, environment variables are set "once and for all" at the beginning of the user's session (on Linux this is usually done in the ``~/.bashrc`` file, or equivalent), and are read on demand by each program that needs them.

By convention and for portability issues, it is recommended to specify environment variable names in full upper case (i.e., ``"PATH"`` and not ``"Path"`` or ``"path"``).

**Example:**

.. code-block:: c++

    // One typical use case is to get the path of some external component
    std::string sed_library_dir = system_var("SUPERFIT_LIB_PATH");
    if (!sed_library_dir.empty()) {
        // The directory has been provided, look what is inside...
    } else {
        // This component is missing, try to do without or print an error
    }

    // It can also be used to modify generic behaviors, for example
    // configure how many threads we want to use by default in all the
    // programs of a given tool suite.
    uint_t nthread = system_var<uint_t>("MYTOOLS_NTHREADS", 1);
