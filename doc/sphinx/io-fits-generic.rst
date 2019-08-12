.. _FITS files:

FITS files
==========

Defined in header ``<vif/io/fits.hpp>``.

A few words on the FITS format
------------------------------

The FITS (Flexible Image Transport System) format is a general purpose file format originally developed for astrophysics data. In particular, FITS files can store images with integer or floating point pixel values, image cubes (with more than two dimensions), but also binary data tables with an arbitrary number of columns and rows. Using a meta-data system (keywords in a header), FITS files usually carry a number of important additional informations about their content. For example, for images files, this can be the mapping between image pixels and sky coordinates (WCS coordinates), or the physical unit of the pixel values.

A single file can contain any number of Header-Data Units (HDUs), which can be seen as independent "extensions" of the file. Each such extension can contain either image data or table data, and has its own header and keywords for storing meta-data. The first HDU of a FITS file has a special status and is called the "primary HDU"; it can only contain image data. For this reason FITS tables always have an empty primary HDU, and the table data located in the first extension.

Writing tabulated data in a binary FITS file is a space-efficient and fast way to store and read non-image data, vastly superior to using human-readable ASCII tables. FITS tables come in two fashions: row-oriented and column-oriented tables. In row-oriented tables, all the data about one row (e.g., about one galaxy in the table) is stored contiguously on disk. This means that it is very fast to retrieve all the information about a given object. In column-oriented tables however, a whole column is stored contiguously in memory. This means that it is very fast to read a given column for all the objects in the table. This distinction is analogous to the dilemma of choosing between a structure-of-array (column-oriented) or an array-of-structures (row-oriented).

Since vif vectors are also contiguous in memory and are used to store data from a given column, the column-oriented format is the most efficient, and is therefore the default format in vif. An additional benefit of this format is that it allows storing columns of different lengths, which is particularly useful to carry meta-data that would be hard to store in FITS keywords. The column-oriented format is not well known, but most softwares and libraries do support it. Topcat_ does, and in IDL column-oriented FITS files are supported naturally by the mrdfits_ and mwrfits_ procedures. But since row-oriented files are nevertheless very common, vif is capable of reading and writing in both formats.

.. _Topcat: http://www.star.bris.ac.uk/~mbt/topcat/
.. _mrdfits: https://www.harrisgeospatial.com/docs/mrdfits.html
.. _mwrfits: https://www.harrisgeospatial.com/docs/mwrfits.html


The FITS implementation in vif
------------------------------

Because FITS files can have a complex structure, vif offers *three* different interfaces to interact with them, which offer different levels of abstraction, capabilities, and ease of use.

The lowest level interface is the raw C interface of the CFITSIO_ library. This interface allows the finest control over the FITS file structure, however it is more cumbersome to use and more error prone. Unless you need to do something very specific that the other interfaces cannot achieve, it is not recommended to use the CFITSIO interface directly. If you do need to use it, please refer to the official CFITSIO documentation.

The second interface already has a significantly higher abstraction level. It is available through the classes ``fits::input_table``,  ``fits::output_table``, ``fits::table``, ``fits::input_image``, ``fits::output_image``, and ``fits::image``. This is an object-oriented interface, where an instance of one of the classes above represents an open FITS file, through which you can navigate to read and/or write data. This is versatile enough to allow you to create multiple table or image extensions, modify data inside the file (i.e., one pixel, or one cell in a table), and edit individual header keywords.

The last and simplest interface is provided through the free functions ``fits::read_table()``, ``fits::write_table()``, ``fits::update_table()``, ``fits::read_image()``, ``fits::write_image()``, and ``fits::update_image()``. These allow you to read/write/update data from a single extension in a FITS file, all with a single function call. They are most convenient when dealing with simple FITS files, however they are less powerful than the object-oriented interface described above.

The object-oriented interface is implemented directly on top of the CFITSIO C API, and the simple free-function interface is implemented on top of the object-oriented interface. These interfaces are not based on the official CFITSIO C++ wrapper, CCFITS_, mostly because its level of abstraction is too high such that it was not possible to implement all the required features with it. This official wrapper could have simply been adopted as the default FITS implementation in vif, however its design choices conflicted too blatantly with the vif "mindset" (different choice of container types, class hierarchy too deep, and notion of HDU too central).

.. _CFITSIO: https://heasarc.gsfc.nasa.gov/fitsio/
.. _CCFITS: https://heasarc.gsfc.nasa.gov/fitsio/CCfits/


General information on the FITS classes
---------------------------------------

The class hierarchy is roughly modeled around the ``std::iostream`` interface:

- ``impl::fits_impl::file_base``. File manipulation (open/close), access extensions, read header data.
    - ``fits::input_image``. Read image data.
    - ``fits::input_table``. Read table data.
    - ``impl::fits_impl::output_file_base``. Edit extensions, edit header keywords.
        - ``fits::output_image``. Write image data.
        - ``fits::output_table``. Write table data.
            - ``fits::input_file``. Read image and table data.
            - ``fits::output_file``. Write image and table data.
            - ``fits::image``. Read/write/update image data.
            - ``fits::table``. Read/write/update table data.
                - ``fits::file``. Read/write/update image and table data.

A file can be opened directly by providing the file name in the constructor, or using the ``open()`` member function. Memory and other resources are freed automatically in the destructor of the object, however it is possible to close the file early if needed using the ``close()`` member function. When the file is open, classes ``fits::input_image``, ``fits::output_image``, and ``fits::image`` will automatically go to the first extension of the file containing image data, likewise with ``fits::input_table``, ``fits::output_table``, and ``fits::table`` and table data.

Any invalid operation will raise an exception of the type ``fits::exception`` (with ``fits::exception::what()`` providing a text message describing the error). It is safe to recover from exceptions raised when attempting invalid data read operations (incorrect image format, unknown table columns, ...). However, write operations are not exception-safe; if they happen, the only safe course of action is to close the file and consider its content as corrupted.

If an instance of any of these classes is shared by multiple execution threads, all operations (even simple queries about the state of the file) must be protected by a mutex lock. In contrast, multiple instances can be used in multiple threads without locking as long as:

- each instance is used by a single thread only,
- all instances are pointing to different files, or instances pointing to the same file are all performing read operations only.


file_base::open
---------------

.. code-block:: c++

    void file_base::open(std::string f);

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


file_base::close
----------------

.. code-block:: c++

    void file_base::close() noexcept;

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


file_base::is_open
------------------

.. code-block:: c++

    bool file_base::is_open() const noexcept;

This function checks if a file is currently open.

**Example:**

.. code-block:: c++

    // Create a FITS image object with no opened file yet
    fits::input_image img;
    img.is_open(); // false
    // Open a file
    img.open("my_image.fits");
    img.is_open(); // true


file_base::filename
-------------------

.. code-block:: c++

    const std::string& file_base::filename() const noexcept;

This function returns the name of the currently opened file (or blank if no file is opened).

**Example:**

.. code-block:: c++

    fits::input_image img("my_image.fits");
    img.filename(); // "my_image.fits"


file_base::cfitstio_status
--------------------------

.. code-block:: c++

    int file_base::cfitstio_status() const noexcept;

This function returns the current CFITSIO error code. Only useful for debugging purposes. If no file is currently open, it will return zero.

**Example:**

.. code-block:: c++

    fits::input_image img("my_image.fits");
    img.cfitsio_status(); // most likely 0


file_base::cfitsio_ptr
----------------------

.. code-block:: c++

    fitsfile* file_base::cfitsio_ptr() noexcept;
    const fitsfile* file_base::cfitsio_ptr() const noexcept;

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


file_base::update_internal_state
--------------------------------

.. code-block:: c++

    void file_base::update_internal_state();

This function is called internally by ``open()`` and ``reach_hdu()``, and is used to update the internal state of the C++ wrapper based on the current content of the file. You only need to use this function if you perform operations on the file using the raw CFITSIO interface. See ``cfitsio_ptr()`` for more information. Will throw an exception if no file is currently open.


output_file_base::flush
-----------------------

.. code-block:: c++

    void output_file_base::flush();

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


output_file_base::flush_buffer
------------------------------

.. code-block:: c++

    void output_file_base::flush_buffer();

This function will perform any pending write operation to the disk and only return when all the data has been written. Contrary to ``flush()``, it will only flush the binary data, and not the header data. This will be faster but less complete; only use this if you know the header data is likely to already be up-to-date. See ``flush()`` for more information. Only available for output files. Will throw an exception if no file is currently open.


file_base::hdu_count
--------------------

.. code-block:: c++

    uint_t file_base::hdu_count() const;

This function returns the number of HDU (or extensions) currently present in the file. This includes the "primary HDU" (extension with ID ``0``), and therefore should always be larger or equal to one. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.hdu_count(); // 1 (only the primary HDU)
    // Reach some other HDU
    img.reach_hdu(1);
    img.hdu_count(); // 2


file_base::current_hdu
----------------------

.. code-block:: c++

    uint_t file_base::current_hdu() const;

This function returns the ID of the current HDU (or extensions). The "primary HDU" has ID of ``0``, and every following HDU has its ID incremented by one. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.current_hdu(); // 0 (the primary HDU)
    // Reach some other HDU
    img.reach_hdu(1);
    img.current_hdu(); // 1


file_base::hdu_type
-------------------

.. code-block:: c++

    fits::hdu_type file_base::hdu_type() const;

This function attempts to identify the content in the current HDU, determining whether it is an image (``fits::image_hdu``), a table (``fits::table_hdu``), or an empty HDU (``fits::empty_hdu``). If it could not decide, it returns ``fits::null_hdu``. The function will throw an exception if the header contains keywords with invalid values, or if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.hdu_type(); // fits::empty_hdu (the primary HDU is initially empty)
    // Write some data
    img.write(data);
    img.hdu_type(); // fits::image_hdu


file_base::reach_hdu
--------------------

.. code-block:: c++

    void file_base::reach_hdu(uint_t hdu);

This function attempts to reach the requested HDU to start reading/writing data from/to it. If this HDU does not exist and the file was opened only with read access, the function will throw an exception. If the file was opened with write access, the function will insert as many empty HDUs as required so that the requested HDU exists, and then reach it for read/write operations. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing; we start at the primary HDU (ID 0)
    fits::output_image img("my_image.fits");
    // Reach some other HDU
    img.reach_hdu(2);
    // Write data there
    vec2d data(10,10);
    img.write(data);
    // The file now contains:
    //  - an empty primary HDU (ID 0)
    //  - an empty first extension (ID 1)
    //  - the image data in the second extension (ID 2)


output_file_base::remove_hdu
----------------------------

.. code-block:: c++

    void file_base::remove_hdu();

This function removes the current HDU from the file. If other HDUs existed after the current HDU, their IDs are decreased by one, to fill the gap. This function will throw an exception when attempting to remove the primary HDU, as by definition it cannot be removed. Will throw an exception if no file is currently open. Only available for output files.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing; we start at the primary HDU (ID 0)
    fits::output_image img("my_image.fits");
    // Reach some other HDU
    img.reach_hdu(2);
    // Write some data
    vec2d data(10,10);
    img.write(data)
    // The file now contains:
    //  - an empty primary HDU (ID 0)
    //  - an empty first extension (ID 1)
    //  - the image data in the second extension (ID 2)

    // Move to the HDU 1
    img.reach_hdu(1);
    // Remove it
    img.remove_hdu();
    // The file now contains:
    //  - an empty primary HDU (ID 0)
    //  - the image data in the first extension (ID 1)


file_base::axis_count
---------------------

.. code-block:: c++

    uint_t file_base::axis_count() const;

This function returns the number of axes of the data located in the current HDU. For image data, this is the number of axes (1 for 1D data, 2 for images, 3 for cubes, etc.). For table data and empty HDUs, the function returns zero. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.axis_count(); // 0 (the primary HDU is initially empty)
    // Write some data
    vec2d data(10,10);
    img.write(data);
    img.axis_count(); // 2


file_base::image_dims
---------------------

.. code-block:: c++

    vec1u file_base::image_dims() const;

This function returns the dimensions of the image in the current HDU. If the current HDU is empty or contains table data, this returns an empty vector. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.image_dims(); // {} (the primary HDU is initially empty)
    // Write some data
    vec2d data(8,10);
    img.write(data);
    img.image_dims(); // {8,10}


file_base::has_keyword
----------------------

.. code-block:: c++

    bool file_base::has_keyword(std::string name) const;

This function checks if a given keyword exists in the header of the current HDU. This check is not case-sensitive, and the function automatically supports long keyword names specified with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::input_image img("my_image.fits");
    img.has_keyword("BUNIT"); // does this image have a unit?


file_base::read_keyword
-----------------------

.. code-block:: c++

    template<typename T>
    bool file_base::read_keyword(std::string name, T& value) const;

This function checks if a given keyword exists in the header of the current HDU, and if the keyword exits, attempts to read its value and store it into the variable ``value``. This check is not case-sensitive, and the function automatically supports long keyword names specified with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. If any of these steps fail, the content of ``value`` is unchanged and the function returns ``false``. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::input_image img("my_image.fits");
    std::string unit;
    if (img.read_keyword("BUNIT", unit)) {
        // We know the unit of the image
    }
    double frequency;
    if (img.read_keyword("FREQ", frequency)) {
        // We know the frequency at which the image was obtained
    }


output_file_base::write_keyword, output_file_base::add_keyword
--------------------------------------------------------------

.. code-block:: c++

    template<typename T>
    void output_file_base::write_keyword(std::string name, const T& value); // [1]
    template<typename T>
    void output_file_base::add_keyword(std::string name, const T& value); // [2]

These functions write the given keyword into the header of the current HDU, setting its value to the provided ``value``. If a keyword with this name already exist, function [1] will update its value, while function [2] will simply ignore it and add a new keyword with the same name at the end of the header (it is indeed possible to have multiple keywords with the same name). If the keyword name is longer than 8 characters, CFITSIO will automatically write the keyword with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::output_image img("my_image.fits");
    vec2d data(10,10);
    img.write(data);
    img.write_keyword("BUNIT", "W/m2/sr"); // write a string
    img.write_keyword("FREQ", 1.4e9);      // write a number


output_file_base::remove_keyword
--------------------------------

.. code-block:: c++

    void output_file_base::remove_keyword(std::string name);

This function will remove the first keyword in the header whose name matches the provided string. No error is generated if no such keyword exists. If the keyword name is longer than 8 characters, CFITSIO will automatically write the keyword with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::output_image img("my_image.fits");
    vec2d data(10,10);
    img.write(data);
    img.write_keyword("BUNIT", "W/m2/sr"); // write a string
    img.remove_keyword("BUNIT");           // we changed our mind, remove it
