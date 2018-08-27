# What is vif ?

vif is a C++ library built to provide user friendly data manipulation as offered in interpreted languages like [python], [IDL], [GDL], but with the added benefit of C++, i.e. increased robustness and speed. The library offers the following features:

 - Multidimensional arrays with standard mathematical operators and
   functions, providing easy manipulation of images and tabulated data.
 - A whole library of general purpose functions to modify and analyze these data (reverse, shuffle, reduce, ...).
 - A number of optional modules:
    - a FITS module allowing read/write operations on FITS images and tables,
    - an ASCII module allowing read/write operations on ASCII tables,
    - an astrophysics module providing tools such as cosmological calculations,
   PSF fitting, unit conversions, catalog cross-matching, image manipulation, ...

It relies only on standard C++11, and a couple of well-known libraries (see INSTALL). Some features can optionally use compile-time reflection, which is not [yet] part of the C++ standard, so a small tool ('refgen') is used to fill the gap. These features are non-essential though, and reflection can be completely disabled if not needed.

The library was designed following three core principles:
 1. "Safety first".
 2. "You do not pay for what you do not use".
 3. "Minimize burden to the user".

The first principle requires that the most easily accessible interface is the safest to use, meaning that the library will perform, by default, some basic checks on the data before processing it (bounds checking, etc.). The second principle requires the existence of mechanisms to bypass these safety checks, such that it is always possible to write optimal code when checks are not needed. The third and last principle requires that both interfaces are as simple to use as possible, to minimize the amount of code the user has to write to reach their goal.

Thus far, vif has been a one-man project. Development started in 2013 during my PhD, and has continued until today. The library is now mature and fully functional; all that remains to be done before a version 1.0 is a final review of the code before API freeze, and wrapping up the documentation (see [roadmap]).

[python]: https://www.python.org/
[IDL]: http://www.exelisvis.com/ProductsServices/IDL.aspx
[GDL]: http://gnudatalanguage.sourceforge.net/
[cling]: http://root.cern.ch/drupal/content/cling
[yet]: https://groups.google.com/a/isocpp.org/forum/#!forum/reflection
[roadmap]: https://github.com/cschreib/vif/projects/1


# Installing the library

Detailed instructions are provided in the ``INSTALL.md`` file. The library is header-only, so installing it is very easy, the main difficulty will be to install the (optional) dependencies.


# Using the library

 - create a new C++ file (foobar.cpp for example)
 - write some code:

   ```cpp
       // Include the vif headers
       #include <vif.hpp>

       // Declare the main function, entry point of the program
       int vif_main(int argc, char* argv[]) {
           // Read some arguments from the command line
           vif::uint_t n = 10;
           vif::read_args(argc, argv, arg_list(n));

           // Create a vector containing 1, 2, ..., n-1
           vif::vec1d v = vif::dindgen(n);

           // Print the square root of each element
           vif::print(sqrt(v));

           return 0;
       }
   ```

 - invoke your compiler from the directory where the source code is
   located:

   ```bash
      g++ -std=c++11 foobar.cpp -o foobar
   ```

 - once the program is compiled, you can run it

   ```bash
      ./foobar n=20
   ```

 - alternatively, you may also compile the program using ``cvif``, a handy shortcut
   script that does the same thing as above (plus it also enables useful warnings,
   and takes care of specifying include paths and linking dependencies):

   ```bash
      cvif foobar.cpp
      ./foobar n=20

      # bonus: compile with full optimizations
      cvif optimize foobar.cpp
      ./foobar n=20
   ```
