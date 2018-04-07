Command line arguments
======================

read_args
---------

The biggest problem of C++ programs is that, while they are fast to execute, they can take a long time to *compile*. For this reason, C++ is generally not a productive language in situations where the code has to be written by *trial and error*, a process that involves frequently changing the behavior or starting point of a program.

There are ways around this issue. One in particular is called *data driven* programming: the behavior of a program depends on the data that are fed to it. The simplest way to use this paradigm is to control the program through "command line arguments". The C++ language provides the basic bricks to use command line arguments, but the interface is inherited from C and lacks severely in usability.

For this reason we introduce in phy++ a single function, ``read_args()``, that uses these bricks to provide a simple and concise interface to implement command line arguments in a program. The first two arguments of the function (``argc, argv``) must be the arguments of the ``main()`` function, in the same order. The following argument must be ``arg_list(...)``, inside of which one must list all the variables that can be modified through the command line interface. These variables can be of any type, as long as it is possible to convert a string into a variable of this type.


Usage examples
--------------

Let us illustrate this with an example. Assume that we want to build a simple program that will print to the terminal the first ``n`` powers of two, with ``n`` being specified by the user of the program. Here is how this would be done with ``read_args()``:

.. code-block:: c++

    # include <phypp.hpp>

    // This is the standard entry point of every C++ program.
    // The signature of the main function is imposed by the C++ standard
    int main(int argc, char* argv[]) {
        // Declare the main parameters of the program, in this case
        // the number of powers of two to display, 'n'.
        uint_t n = 1;

        // Then read command line arguments...
        read_args(argc, argv, arg_list(n));

        // Now we can go on with the program, using 'n' normally
        print(pow(2.0, findgen(n)+1));

        // And quit gracefully
        return 0;
    }


By just adding this line

.. code-block:: c++

    read_args(argc, argv, arg_list(n));

we exposed the variable ``n`` to the public: everyone that runs this program can modify the value of ``n`` to suit their need. Simple! Assuming the name of the compiled program is ``show_pow2``, then the program is ran the following way:

.. code-block:: bash

    # First try with no parameter.
    # 'n' is not modified, and keeps its default value of 1.
    ./show_pow2
    # output:
    # 2

    # Then we change 'n' to 5.
    ./show_pow2 n=5
    # output:
    # 2 4 8 16 32

The advantages of this approach are immediate. Instead of recompiling the whole program just to change ``n``, we exposed it in the program arguments. We then compiled the program *once*, and changed its behavior without ever recompiling. This can save a lot of time, for example when trying to figure out what is the best value of a parameter in a given problem (i.e., tweaking parameters of an algorithm). And of course, this is most useful when writing *tools* with configurable options.

Within the ``arg_list()``, one can put as many C++ variables as needed. The function will recognize their name, and when the program is executed it will try to find a command line argument that matches. If it finds one, it tries to convert its value into the type of the variable, and if successful, store this new value inside the variable. In all other cases, the variable is not modified. It is therefore important to give a meaningful default value to each variable!

In the example above, we chose to expose a simple integer. But in fact, the interface can be used to expose any type, provided that there is a way to convert a string into a value of this type. In particular, this is the case for vectors. The values of the vector must be separated by commas, (*without any space*, unless you put the whole argument inside quotes), and surrounded by brackets ``[...]``. Again, let us illustrate this with an example. We will modify the previous program to allow it to show not only powers of ``2``, but the powers of multiple, arbitrary numbers. Note: in the following, we will not repeat the whole ``main()`` function, just the important bits.

.. code-block:: c++

    // The number of powers of two to display
    uint_t n = 1;
    // The powers to display
    vec1f p = {2};

    // Read command line arguments
    read_args(argc, argv, arg_list(n, p));

    // Go on with the program
    for (float v : p) {
        print(pow(v, findgen(n)+1));
    }


The program can now change the powers it displays, for example:

.. code-block:: bash

    # We keep 'n' equal to 5, and we show the powers of 2, 3 and 5.
    ./show_pow2 n=5 p=[2,3,5]
    # output:
    # 2 4 8 16 32
    # 3 9 27 81 243
    # 5 25 125 625 3125

    # It is possible to use spaces inside the [...], but then you must add quotes:
    ./show_pow2 n=5 p="[2, 3, 5]"


Now, you may think that ``p`` is not a very explicit name for this last parameter. It would be clearer if we could call it ``pow``. Unfortunately, ``pow`` is already the name of a function in C++, so we cannot give this name to the variable. However, the ``read_args()`` interface allows you to manually give a name to any parameter using the ``name()`` function. Let us do that and modify the previous example.

.. code-block:: c++

    // The number of powers of two to display
    uint_t n = 1;
    // The powers to display, we still call it 'p' in the program
    vec1f p = {2};

    // Read command line arguments
    read_args(argc, argv, arg_list(n, name(p, "pow"));

    // Go on with the program
    for (float v : p) {
        print(pow(v, findgen(n)+1));
    }


Now we will write instead:

.. code-block:: bash

    ./show_pow2 n=5 pow=[2,3,5]
    # output:
    # 2 4 8 16 32
    # 3 9 27 81 243
    # 5 25 125 625 3125


Flags
-----

Often, command line options are "flags". These are boolean variables that are ``false`` by default, but can be changed to ``true`` to enable some specific functionality. For example, setting ``verbose=1`` can be used to tell the program to display information in the terminal about its progress. To simplify usage of these flags, ``read_args()`` allows an alternative syntax where specifying ``verbose`` without any equal sign in the arguments is equivalent to ``verbose=1``:

.. code-block:: bash

    ./my_program verbose
    # ... is equivalent to:
    ./my_program verbose=1

There is no shortcut for ``var=0``.


Alternative syntax
------------------

In the examples above, command lines arguments are specified as ``variable=value``. This is the tersest available syntax. However, most linux programs tend to use dashes (``-``) to identify command line arguments; for example ``-variable=value`` or ``--variable=value``. To avoid confusing users, ``read_args()`` supports both ways of writing command line arguments; dashes can be used but are not mandatory.
