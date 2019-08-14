.. _FITS files:

Managing HDUs
=============

Defined in header ``<vif/io/fits.hpp>``.

hdu_count
---------

.. code-block:: c++

    uint_t fits::file_base::hdu_count() const;

This function returns the number of HDUs (or extensions) currently present in the file. This includes the "primary HDU" (extension with ID ``0``), and therefore should always be larger or equal to one. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.hdu_count(); // 1 (only the primary HDU)
    // Reach some other HDU
    img.reach_hdu(1);
    img.hdu_count(); // 2


current_hdu
-----------

.. code-block:: c++

    uint_t fits::file_base::current_hdu() const;

This function returns the ID of the current HDU (or extension). The "primary HDU" has ID of ``0``, and every following HDU has its ID incremented by one. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.current_hdu(); // 0 (the primary HDU)
    // Reach some other HDU
    img.reach_hdu(1);
    img.current_hdu(); // 1


hdu_type
--------

.. code-block:: c++

    fits::hdu_type fits::file_base::hdu_type() const;

This function attempts to identify the content in the current HDU, determining whether it is an image (``fits::image_hdu``), a table (``fits::table_hdu``), or an empty HDU (``fits::empty_hdu``). If it could not decide, it returns ``fits::null_hdu``. The function will throw an exception if the header contains keywords with invalid values, or if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image for writing
    fits::output_image img("my_image.fits");
    img.hdu_type(); // fits::empty_hdu (the primary HDU is initially empty)
    // Write some data
    img.write(data);
    img.hdu_type(); // fits::image_hdu


reach_hdu
---------

.. code-block:: c++

    void fits::file_base::reach_hdu(uint_t hdu);

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


remove_hdu
----------

.. code-block:: c++

    void fits::file_base::remove_hdu();

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


axis_count
----------

.. code-block:: c++

    uint_t fits::file_base::axis_count() const;

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


image_dims
----------

.. code-block:: c++

    vec1u fits::file_base::image_dims() const;

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


has_keyword
-----------

.. code-block:: c++

    bool fits::file_base::has_keyword(std::string name) const;

This function checks if a given keyword exists in the header of the current HDU. This check is not case-sensitive, and the function automatically supports long keyword names specified with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::input_image img("my_image.fits");
    img.has_keyword("BUNIT"); // does this image have a unit?


read_keyword
------------

.. code-block:: c++

    template<typename T>
    bool fits::file_base::read_keyword(std::string name, T& value) const;

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


write_keyword, add_keyword
--------------------------

.. code-block:: c++

    template<typename T>
    void fits::output_file_base::write_keyword(std::string name, const T& value); // [1]
    template<typename T>
    void fits::output_file_base::add_keyword(std::string name, const T& value); // [2]

These functions write the given keyword into the header of the current HDU, setting its value to the provided ``value``. If a keyword with this name already exist, function [1] will update its value, while function [2] will simply ignore it and add a new keyword with the same name at the end of the header (it is indeed possible to have multiple keywords with the same name). If the keyword name is longer than 8 characters, CFITSIO will automatically write the keyword with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::output_image img("my_image.fits");
    vec2d data(10,10);
    img.write(data);
    img.write_keyword("BUNIT", "W/m2/sr"); // write a string
    img.write_keyword("FREQ", 1.4e9);      // write a number


remove_keyword
--------------

.. code-block:: c++

    void fits::output_file_base::remove_keyword(std::string name);

This function will remove the first keyword in the header whose name matches the provided string. No error is generated if no such keyword exists. If the keyword name is longer than 8 characters, CFITSIO will automatically write the keyword with the ``HIERARCH`` convention; it is not necessary to specify the ``HIERARCH`` explicitly. Will throw an exception if no file is currently open.

**Example:**

.. code-block:: c++

    // Open a FITS image
    fits::output_image img("my_image.fits");
    vec2d data(10,10);
    img.write(data);
    img.write_keyword("BUNIT", "W/m2/sr"); // write a string
    img.remove_keyword("BUNIT");           // we changed our mind, remove it
