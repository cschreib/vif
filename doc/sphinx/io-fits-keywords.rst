.. _FITS files:

Header & keywords
=================

Defined in header ``<vif/io/fits.hpp>``.


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
