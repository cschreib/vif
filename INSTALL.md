Pre-requisites
==============

1) You must have a recent compiler to compile C++ code. Either ``clang`` (>= 3.3)
   or ``gcc`` (>= 4.7) will do.

2) You should have ``cmake`` (>= 2.6) installed (see Appendix below for instructions on how
   to install it).

3) Installing the ``cfitsio``, ``wcslib``, ``lapack``, ``gsl``, ``fftw-3``, ``tcmalloc``,
   ``libdwarf``, ``libunwind``, and ``libelf`` libraries is not mandatory, but recommended
   to use all the features of the vif library (see Appendix below for instructions on how
   to install these libraries). The ``libclang`` is needed if you want to use reflection.
   In Ubuntu or other Debian-based distributions:

   ```bash
   sudo apt-get install libcfitsio3-dev wcslib-dev liblapack-dev libgsl0-dev libfftw3-dev libgoogle-perftools-dev libdwarf-dev libunwind8-dev libelf-dev
   ```

   Here are the features that will not be available if you do not install some of these libraries:
    - code backtraces in case of errors (requires ``libunwind``, ``libelf``, ``libdwarf``).
    - Fourier transforms (requires ``fftw``).
    - Exotic math functions like incomplete gamma functions (requires ``gsl``).
    - Eigenvalue/eigenvector decomposition and faster matrix inversion (requires ``lapack``).
    - FITS file input/output (requires ``cfitsio``).
    - WCS coordinate system conversions (requires ``cfitsio`` and ``wcslib``).

4) Installing the GNU debugger ``gdb`` is recommended but not mandatory.


Install the library with CMake
==============================

1) Make a new directory called ``build`` inside the main vif directory and
   navigate to it using a terminal.

2) In this directory, call ``cmake ../``.
   If you have properly installed all the required libraries, CMake will find
   them and configure the project accordingly.

   If you have installed some dependencies in non-standard directories, the
   CMake script will not be able to find them. Therefore you need to specify
   the install directory manually for each of these libraries. To do so, just
   add the following command line arguments to the ``cmake ../`` command:

    - for cfitsio: ``-DCFITSIO_ROOT_DIR=...``
    - for wcslib: ``-DWCSLIB_ROOT_DIR=...``
    - for lapack: ``-DLAPACK_ROOT_DIR=...``
    - for libclang: ``-DCLANG_ROOT=...``
    - for gsl: ``-DGSL_DIR=...``
    - for fftw: ``-DFFTW_ROOT=...``
    - for tcmalloc: ``-DTCMALLOC_ROOT_DIR=...``

   The ``...`` have to be replaced by the path to where the library was
   installed, i.e. it must contain the ``include`` and ``lib`` folders in which
   vif will look for headers and shared libraries, respectively.

   For example, if cfitsio was installed in ``/opt/local/share/cfitsio/``, then
   you will have to call ``cmake`` using:

   ```bash
   cmake ../ -DCFITSIO_ROOT_DIR=/opt/local/share/cfitsio/
   ```

3) Build the library by calling ``make`` (will do nothing if reflection is
   disabled and/or libclang was not found).

4) Install the library by calling ``sudo make install``.

5) As the install script suggests, you should then source the ``.vifrc`` file
   that was just created in your home directory, so that the vif configuration
   is loaded in all your terminals. This can be done by adding a line in your
   ``.bashrc`` file containing ``source ~/.vifrc``.

6) You may then build the vif tools. Inside the ``build`` directory, create & navigate to
   the ``tools`` directory. Note: if you are still in the same terminal session, you
   will need to either close it and open a new one, or source the vif configuration
   using ``source ~/.vifrc``

7) Invoque cmake again, this time pointing to the tools directory:
   ``cmake ../../tools/``. Like in point 2) above, you may need to manually
   specify the location of some dependencies if they were not installed in
   standard paths.

8) Build and install, using ``make`` then ``sudo make install``.


Veryfing installation
=====================

1) Open a new terminal so that the vif configuration is loaded.

2) Go to the ``test`` subfolder.

3) Type ``vif unit_test``.

4) This should print a series of tests, eventually saying that all tests were
   passed. If this is not the case, then either the library has not been
   properly installed, or your compiler is either old or broken.


Appendix: installing dependencies (ubuntu and derivatives)
==========================================================

Installing the required libraries in ubuntu-based systems is easy. You first
need to run the following command:

   ```bash
  sudo apt-get install cmake
   ```

If you want, you can also install the optional (but recommanded) libraries and
tools:

   ```bash
  sudo apt-get install libcfitsio3-dev liblapack-dev wcslib-dev libgsl0-dev libgoogle-perftools-dev libfftw3-dev gdb libdwarf-dev libunwind8-dev libelf-dev
   ```

Lastly, if you want to use reflection, you will need to get a recent version
of libclang.

If running ubuntu Trusty Tarh (14.04) and above:

   ```bash
  sudo apt-get install libclang-3.5-dev
   ```

You may have to increase the version number if your distribution has moved to a newer version of the library. You may install any other version provided it is >= 3.3, however note that the
official Ubuntu packages for versions 3.3 to 3.4 are sometimes broken.


Appendix: installing dependencies (from source)
===============================================

Alternatively, you may also build all the dependencies from source. This will
take longer, and may require that you install sub-dependencies, but should work
if everything else fails.

Since the build system of all these libraries can vary, please refer to each
libraries' documentation for precise build instructions. Note that I will not
offer support on how to compile and install these libraries, you will have to
contact the libraries' respective authors.

cfitsio:
http://heasarc.nasa.gov/fitsio/fitsio.html

wcslib:
http://www.atnf.csiro.au/people/mcalabre/WCS/

lapack:
http://www.netlib.org/lapack/

libclang (llvm toolchain):
http://llvm.org/releases/download.html

gsl:
http://www.gnu.org/software/gsl/

tcmalloc (google perf tools):
http://code.google.com/p/gperftools/

fftw3:
http://www.fftw.org/download.html

libdwarf:
https://github.com/tomhughes/libdwarf

libelf:
https://directory.fsf.org/wiki/Libelf

libunwind:
https://www.nongnu.org/libunwind/

google-perftools:
https://github.com/gperftools/gperftools
