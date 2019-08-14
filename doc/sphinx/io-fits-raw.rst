.. _FITS files:

C interface
===========

Defined in header ``<vif/io/fits.hpp>``.

cfitstio_status
---------------

.. code-block:: c++

    int fits::file_base::cfitstio_status() const noexcept;

This function returns the current CFITSIO error code. Only useful for debugging purposes. If no file is currently open, it will return zero.

**Example:**

.. code-block:: c++

    fits::input_image img("my_image.fits");
    img.cfitsio_status(); // most likely 0


cfitsio_ptr
-----------

.. code-block:: c++

    fitsfile*       fits::file_base::cfitsio_ptr()       noexcept;
    const fitsfile* fits::file_base::cfitsio_ptr() const noexcept;

These functions returns the underlying CFITSIO file pointer. This is useful if you need to perform an operation that is not available as part of the C++ interface. It is safe to perform any operation with this pointer and then fall back to the C++ interface, however if you do so you must call the ``update_state()`` function before using any function of the C++ interface.

If no file is currently open, it will return a null pointer.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::input_image img("my_image.fits");
    // Get the underlying CFITSIO pointer
    fitsptr* fptr = img.cfitsio_ptr();
    // Use the pointer with the raw C interface
    // ...
    // Update the internal state
    img.update_internal_state();
    // Continue using the C++ interface


update_internal_state
---------------------

.. code-block:: c++

    void fits::file_base::update_internal_state();

This function is called internally by ``open()`` and ``reach_hdu()``, and is used to update the internal state of the C++ wrapper based on the current content of the file. You only need to use this function if you perform operations on the file using the raw CFITSIO interface. See ``cfitsio_ptr()`` for more information. Will throw an exception if no file is currently open.
