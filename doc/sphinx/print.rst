Printing to the terminal
========================

Defined in header ``<phypp/core/print.hpp>``.

print
-----

.. code-block:: c++

    template<typename ... Args>
    void print(Args&& ...); // [1]

This function will "print" (or display) the content of all its arguments into the standard output (i.e., the terminal). This has *nothing* to do with printing things on paper (unless you have set *that* as your standard output). Arguments are printed one after the other on the same line, without any spacing or separator, and the "end of line" character is printer after the last argument so that the next call of ``print()`` will display on a new line. This function should only be used for debugging purposes, or to inform the user of the progress of the program (but see also below and ``print_progress()``).

**Example:**

.. code-block:: c++

    uint_t a = 5;
    vec1u x = {1,2,3,4,9,10};

    print("the program reached this location, a=", a, ", and x=", x);

The above code will display on the terminal:

.. code-block:: none

    the program reached this location, a=5, and x={1,2,3,4,9,10}


**A few words about formatting.**

First, since the text is meant to be sent to a simple terminal window, there are very few options to affect the look and feel of the output. In fact there is currently no option to change fonts, add boldface or italic, or change the color (some of this may be implemented in the future, but it will always remain fairly limited). The only thing you can actually control is when to start a new line. This is done with the special character ``'\n'``, and such a line break is always made at the end of each printed message. If you need a finer control of line jumps, you will have to come back to the standard C++ functions for output (i.e., ``std::cout`` and the likes).

Second, if you are used to C functions like ``prinft()``, note from the above example how printing in phy++ behaves quite differently. The ``print()`` function does not work with "format strings", hence the character ``'%'`` is not special and can be used safely without escaping. One can print a value that is not a string simply by adding it to the list of the function arguments, which will do the conversion to string automatically based on that argument's type. See also :ref:`String conversions` for more information on that process. The format of this conversion is defined by the C++ standard library, and is usually fine for most cases. If you need a finer control on the display format (for example for floating point numbers), look at :ref:`Formatting` and ``format::`` functions.


error, warning, note
--------------------

.. code-block:: c++

    template<typename ... Args>
    void error(...); // [1]

    template<typename ... Args>
    void warning(...); // [2]

    template<typename ... Args>
    void note(...); // [3]

These functions behave exactly like ``print()``. The only difference is that they automatically append a "prefix" to the message, namely: ``"error: "`` ([1]), ``"warning: "`` ([2]), or ``"note: "`` ([3]). In some operating systems, these prefixes can be colored to make them stand out better from the regular ``print()`` output. These three functions are indeed much more useful than ``print()``, since they are meant to talk to the *user* of a program rather than the *developer*. They can be used, for example, to tell the user that something went wrong ([1]), or could go wrong ([2]), or simply to inform them of what's going on ([3]).

More specifically, function [1] (``error()``) should be used to print *unrecoverable* and *serious* errors. In other words, this is a situation the program was not designed to handle, and it has to stop and tell the user why. Function [2] (``warning()``) is for *recoverable* errors (unusual situations for which a workaround is implemented), or situations where the user is *likely* to have made a mistake (but maybe not). Lastly, function [3] (``note()``) is for everything else which is not critical and does not require immediate attention. Typically this will be logging, displaying the progress of a calculation, etc. Since these types of messages are not critical, the user should be able to quickly glance over them and ignore them. Is it therefore good practice to offer an option to completely disable this non-critical output, so only errors and warnings are displayed.

**Example:**

.. code-block:: c++

    bool verbose = false; // let the user enable extra output only if needed

    std::string datafile1 = "toto1.txt";
    if (!file::exist(datafile1)) {
        // We cannot work without 'datafile1', print an error
        error("cannot open '", datafile1, "'");
        error("make sure the program is run inside the data directory");
        return 1; // typically, exit the program
    }

    std::string datafile2 = "toto2.txt";
    if (!file::exist(datafile2)) {
        // It's better is we have 'datafile2', but we can work without it.
        // So we print a warning and let the user know what are the
        // consequences of this non-critical issue.
        warning("cannot open '", datafile2, "', so calculation will be less accurate");
        warning("will use '", datafile1, "' as a fallback");
        datafile2 = datafile1;
    }

    if (verbose) {
        // Things are going fine, inform the user of what we are about to do
        note("analysing the data, please wait...");
    }

    // Do the stuff...

The above code, if the first file does not exist, will display:

.. code-block:: none

    error: cannot open 'toto1.txt'
    error: make sure the program is run inside the data directory


prompt
------

.. code-block:: c++

    template<typename T>
    bool prompt(const std::string& msg, T& v, const std::string& err = ""); // [1]

This function interacts with the user of the program through the standard output and input (i.e., the terminal). It first prints ``msg``, waits for the user to enter a value and press the Enter key, then try to load this value inside ``v``. If the value entered by the user is invalid and cannot be converted into the type ``T``, the program asks again and optionally writes an error message ``err`` to clarify the situation.

Currently, the function can only return after successfully reading a value, and always returns ``true``. In the future, it may fail and return ``false``, for example after the user has failed a given number of times. If possible, try to keep the possibility of failure into account.

**Example:**

.. code-block:: c++

    uint_t age;
    if (prompt("please enter your age: ", age,
        "it better just be an integral number...")) {
        print("your age is: ", age);
    } else {
        print("you will do better next time");
    }

Here is a possible interaction scenario with a (naive) user:

.. code-block:: none

    please enter your age: 15.5
    error: it better just be an integral number...
    please enter your age: what?
    error: it better just be an integral number...
    please enter your age: oh I see, it is 15
    error: it better just be an integral number...
    please enter your age: ok...
    error: it better just be an integral number...
    please enter your age: 15
    your age is: 15
