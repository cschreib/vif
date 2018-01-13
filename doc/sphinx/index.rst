phy++ documentation
===================

phy++ is a set of library and tools built to provide user-friendly vector data manipulation, as offered in interpreted languages like IDL_, its open source clone GDL_, or Python and numpy_, but with the added benefit of C++: increased robustness, and optimal speed.

The library can be split into two components:

* The "core" library
* The "support" library

The core library introduces the ``vec`` type (a "data vector"), which is the most important data type in phy++, while the support library provides functions and other tools to manipulate these vectors and perform common tasks. You can think of phy++ as a separate language inside C++, where the core library defines this language, and the support library is the the "standard" library where all the useful functions are stored.

Below is a code sample written in phy++ that illustrates the most basic functionalities.

.. code-block:: c++

    vec2f img = fits::read("img.fits"); // read a FITS image
    img -= median(img);                 // subtract the median of the whole image
    float imax = max(img);              // find the maximum of the image
    vec1u ids = where(img > 0.5*imax);  // find pixels at least half as bright
    float sum = total(img[ids]);        // compute the sum of these pixels
    img[ids] = log(img[ids]/sum);       // modify these pixels with a logarithm
    fits::write("new.fits", img);       // save the modified image to a FITS file

.. _IDL: http://www.exelisvis.com/ProductsServices/IDL.aspx
.. _GDL: http://gnudatalanguage.sourceforge.net/
.. _numpy: http://www.numpy.org/

.. _core-library:

.. toctree::
   :maxdepth: 2
   :caption: The core library

   overview
   vector
   indexing
   view
   generic-function-guidelines

.. _support-library:

.. toctree::
   :maxdepth: 2
   :caption: The support library

   generic
   argv
   string
   os
   io
   print
   time
   math
   thread
   image
   astro
