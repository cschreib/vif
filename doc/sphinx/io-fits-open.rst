.. _FITS files:

Open/close files
================

Defined in header ``<vif/io/fits.hpp>``.

open
----

.. code-block:: c++

    void fits::file_base::open(std::string f);

This function will open a new FITS file for read or write access (depending on the actual class that is used).

**Example:**

.. code-block:: c++

    // Open a file directly in the constructor
    fits::input_image img1("my_image.fits");

    // Open a file later using open()
    fits::input_image img2;
    img2.open("my_image.fits");

If there is already a file open when ``open()`` is called, that file is closed before the new file is opened.

When requesting only read access (i.e., with the classes ``fits::input_image`` or ``fits::input_table``), an exception will be raised if the file does not exist or cannot be accessed with read permission. When requesting only write access (i.e., with the classes ``fits::output_image`` or ``fits::output_table``), a new file will be created regardless of whether a file with the provided name already exists or not, and an exception will be raised if the file cannot be created. When requesting read/write access (i.e., with the classes ``fits::image`` or ``fits::table``), an exception will be raised if the file does not exist or cannot be accessed with read/write permission.

This function is partially exception-safe: if a file was previously open before the call and an exception is raised, that file will be closed and no data will be lost. Aside from this minor point, the instance can be used safely after recovering from an exception raised by ``open()``.

It is possible to open the same file multiple times as different objects, but this is not safe when performing write operations. It is, however, perfectly safe to read data from the same file through two objects:

.. code-block:: c++

    // Open the same file twice for reading data
    fits::input_image img1, img2;
    img1.open("my_image.fits");
    img2.open("my_image.fits");
    // Perform read operations (safe)
    vec2d image1, image2;
    img1.read(image1);
    img1.read(image2);


close
-----

.. code-block:: c++

    void fits::file_base::close() noexcept;

This function will close the currently opened FITS file (if any). If data was written to the file, it will be force-flushed to the disk to ensure no data is lost before the file is closed.

This function is called automatically in the destructor, so you do not need to call it explicitly unless you want to close the file before the end of the object's lifetime.

If the file cannot be properly closed for any reason, this function will not raise an exception and simply consider the file as closed.

**Example:**

.. code-block:: c++

    // Open a file
    fits::input_image img("my_image.fits");
    // Perform some operations
    // ...
    // Close the file early
    img.close();
    // A new file must now be opened before doing further operations


is_open
-------

.. code-block:: c++

    bool fits::file_base::is_open() const noexcept;

This function checks if a file is currently open.

**Example:**

.. code-block:: c++

    // Create a FITS image object with no opened file yet
    fits::input_image img;
    img.is_open(); // false
    // Open a file
    img.open("my_image.fits");
    img.is_open(); // true


filename
--------

.. code-block:: c++

    const std::string& fits::file_base::filename() const noexcept;

This function returns the name of the currently opened file (or blank if no file is opened).

**Example:**

.. code-block:: c++

    fits::input_image img("my_image.fits");
    img.filename(); // "my_image.fits"


flush
-----

.. code-block:: c++

    void fits::output_file_base::flush();

This function will perform any pending write operation to the disk and only return when all the data has been written. It will perform a full update of the file, including binary data and header data. Only available for output files. Will throw an exception if no file is currently open.

Indeed, as with any disk write operation in the C++ standard library, CFITSIO write operations use a write buffer which is only written to the disk occasionally, rather than on any write operation. This is done for performance reasons. The downside of this approach is that the data is not always immediately written to the disk, even after a call to ``write()`` has returned. This usually is not an issue, except when one wants to access the content of the file while it is being written, or if the program crashed while data was not yet written to the file.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    // Write some data
    img.write(data);
    // Force writing data to disk now
    img.flush();


flush_buffer
------------

.. code-block:: c++

    void fits::output_file_base::flush_buffer();

This function will perform any pending write operation to the disk and only return when all the data has been written. Contrary to ``flush()``, it will only flush the binary data, and not the header data. This will be faster but less complete; only use this if you know the header data is likely to already be up-to-date. See ``flush()`` for more information. Only available for output files. Will throw an exception if no file is currently open.
