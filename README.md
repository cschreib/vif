# What is phy++ ?

phy++ is a set of library and tools built to provide user friendly data manipulation as offered in interpreted languages like [IDL] or its open source clone [GDL], but with the added benefit of C++, i.e. increased robustness and speed. The library offers the following features:

 - multidimensional arrays with overloaded mathematical operators and
   mathematical functions, providing easy manipulation of tabulated data
 - a whole library of general purpose functions to modify and analyze these data
 - a FITS module allowing read/write operations on FITS images and tables
 - an ASCII module allowing read/write operations on ASCII tables
 - an astrophysics module providing tools such as cosmological calculations,
   PSF fitting, unit conversions, catalog cross-matching, image manipulation, ...
 - a full featured plot module to output graphs in EPS format (WIP!)

It relies only on standard C++11, and a couple of well-known libraries (see INSTALL). A few features require compile-time reflection, which is not [yet] part of the C++ standard, so a small tool ('refgen') is used to fill the gap. These features are non-essential though, so reflection can be completely disabled if not needed.

Currently, there is one main drawback of using phy++ instead of IDL/GDL. Indeed, while I tried to be exhaustive in which functions are provided inside the library, it is certainly missing some of the functions or features available today in IDL. The only way to solve this is to write these missing functions with phy++! You can do it yourself and contribute them back to phy++ by submitting a pull request, or submit a feature request so that I take a look at it.

Another important thing to consider is that, to date, there is no stable and full-featured C++ interpreter (although [cling] is a promising candidate). This means that it is impossible (or hard) to do "real-time" interaction with your code, namely, stop the program and inspect the value of some arrays, plot them, modify them, then continue execution. Any C++ program must be recompiled whenever a change is made, and compilation itself usually takes a few precious seconds. This issue can be partly solved by making use of command line arguments to quickly modify the behavior of a program without having to recompile. But this is obviously more limited. 

For this reason, I recommend the use of the phy++ library for writing "programs" that deal with only one specific task (like reading a specific ASCII catalog into a FITS table, or compute a stellar mass function from a catalog with incompleteness corrections), rather than "scripts" that do many things at once and need to be edited often whenever a new task must be done. This also encourages writing reusable and modular code, that is written once and is generic enough that it can be used in many different situations, which a good practice!

[IDL]: http://www.exelisvis.com/ProductsServices/IDL.aspx
[GDL]: http://gnudatalanguage.sourceforge.net/
[cling]: http://root.cern.ch/drupal/content/cling
[yet]: https://groups.google.com/a/isocpp.org/forum/#!forum/reflection


# Using the 'phypp' library

 - create a new C++ file (foobar.cpp for example)
 - write some code:

```cpp
    // Include the phy++ headers
    #include <phypp.hpp>

    // Declare the main function, entry point of the program
    int phypp_main(int argc, char* argv[]) {
        // Read some arguments from the command line
        uint_t n = 10;
        read_args(argc, argv, arg_list(n));
        // Create a vector containing 1, 2, ..., n-1
        vec1d v = dindgen(n);
        // Print the square root of each element
        print(sqrt(v));
        return 0;
    }
```

 - invoke the phy++ compiler from the directory where the source code is
   located:

```bash
   phy++ foobar n=20
```

 - the program will compile and run automatically (if no compiler error
   occurred)
 - alternatively, you may also just compile the program using the cphy++
   compiler and run the program yourself afterwards:

```bash
   cphy++ foobar.cpp
   ./foobar n=20
```

# The different phy++ "compilers"

Internally, the "phy++ compilers" use your favorite C++ compiler to compile C++ files into binary executables. The various compilers listed below are merely scripts built around the real C++ compiler to simplify the compilation process. In particular:
 - they take care of setting all the necessary compiler flags, include
   directories, linked libraries, enabling useful warnings, etc.,
 - if reflection is enabled, they automatically generate the reflection data using
   the 'refgen' tool,
 - like 'make', they only recompile the code when a change has been made (either
   to the code itself or to any file of the phy++ library)

It is of course possible to compile your code with any other C++ compiler, provided that you set it up properly and that it supports all the language features that this library depends on. Just look at the code for the phy++ scripts to know what options, libraries and includes are needed.

List of compilers:
 - phy++:
    Standard compiler. Should compile fast but will not optimize the executable.
 - ophy++:
    Optimized compiler. Can take twice more time to compile, but will generate faster code.
 - gphy++:
    Debug compiler. As with 'phy++', but compiles with debug information and runs the program through 'gdb'.
 - gophy++:
    Optimized debug compiler. As with 'ophy++', but compiles with debug information and runs the program through 'gdb'. Note that, since optimizations can remove intermediary steps/variables, debugging an optimized program can be troublesome (some variables or functions do disappear). Only use this option when non-optimized binaries are too slow.
 - prophy++: (only if you have installed Google perftools)
    Profiling compiler. As with 'ophy++', but generates profiling data during execution. When the program ends, it parses these data to generate a graphical visualization of the time spent in the various functions, to help you identify the origin of performance issues.
 - cphy++:
    General purpose compiler. Contrary to all the compilers listed above, this one will not execute the program after compilation. Nor will it check for changes in the code or its dependencies: the program is always recompiled. You can set some options, like enabling optimizations, debug information, and specifying the output file name.

As a rule of thumb, prefer using 'ophy++'. If the compilation time becomes large and you need to recompile often for testing purposes, try with 'phy++' until your tests are done, then use 'ophy++' again. In any case, you want to use 'ophy++' (or 'cphy++ optimize') when you have a "final" version of a program, ready to be distributed.
