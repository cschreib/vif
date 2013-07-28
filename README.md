phy++
=====
This is a set of library and tools build to provide user friendly data manipulation as offered in
interpreted languages like IDL [1] or its open source clone GDL [2], but with the added benefit of
C++, i.e. increased robusteness and speed. The library offers the following features:
 - multidimensional arrays with overloaded mathematical operators and mathematical functions,
   providing easy manipulation of tabulated data
 - a whole library of general purpose functions to modify and analyse these data
 - a FITS module allowing read/write operations on FITS images and tables
 - an ASCII module allowing read/write operations on ASCII tables
 - an astrophysics module providing tools such as cosmological calculations, PSF fitting, unit
   conversions, catalog cross-matching, ...
 - a full featured plot module to output graphs in EPS format

Currently, the only drawback of using this library instead of IDL/GDL is that, at the time of
writing, there is no stable and full featured C++ interpreter (although cling [3] is a promising
candidate). This means that it is impossible to do "real-time" testing of code, and that any 
program must be recompiled whenever a change is made. This issue can be partly solved by making use
of command line arguments to quicly modify the behavior of a program.

The library itself relies on standard C++11 (using a few features from the upcoming C++14), and a
small tool ('refgen') is used to parse the source code to generate necessary data for compile-time
reflection (which is not yet part of the C++ standard).

[1]: http://www.exelisvis.com/ProductsServices/IDL.aspx
[2]: http://gnudatalanguage.sourceforge.net/
[3]: http://root.cern.ch/drupal/content/cling


Using the 'phypp' library
========================
1) create a new C++ file (foobar.cpp for example)
2) include the phy++ header: #include <phypp.hpp>
3) write some code in the main function: int main(int argc, char* argv[]) { ... return 0; }
4) invoque the phy++ compiler from the directory where the source code is located: phy++ foobar
5) the program will compile and run automatically (if no compiler error occured)


The different phy++ compilers
=============================
- phy++: standard compiler. Should compile fast but will not fully optimize the executable.
- ophy++: optimized compiled. Can take twice more time to compile, but will generate faster code.
- gphy++: debug compiler. As with 'phy++', but compiles with debug information and fires 'gdb'.
- gophy++: debug compiler. As with 'ophy++', but compiles with debug information and fires 'gdb'.
          Note that, since optimizations can remove intermediary steps/variables, debugging can
          be troublesome (some variables or functions do disappear). Only use this option when
          needed.

As a rule of thumb, prefer using 'phy++' for testing purposes (when you compile often). Then, if 
the execution time gets larger, try with 'ophy++'. In any case, prefer using 'ophy++' when you 
have a "final" version of a program, ready to be distributed.

