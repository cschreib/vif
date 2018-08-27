Processes
=========

Defined in header ``<vif/utility/os.hpp>``.

spawn
-----

.. code-block:: c++

    bool spawn(const std::string& cmd);

This function executes the shell command passed in argument (on UNIX systems, using ``/bin/sh``). It puts execution of the current program on hold until the shell command terminates. This can be used to launch another application, for example an image viewer that the user can use to inspect some temporary calculations.

If you wish to run this command in parallel with the current program, see ``fork()``.

**Example:**

.. code-block:: c++

    // Calculate some things to make an image
    vec2d img = /* ... */;

    // Write this image to the disk in FITS format
    fits::write("tmp.fits", img);

    // Call the FITS viewer DS9 on this data and let user inspect it.
    // The program will be paused until the user closes DS9.
    spawn("ds9 tmp.fits");

    // Now execution comes back to us.
    // We can, for example, ask the user if the data was satisfactory.


fork
-----

.. code-block:: c++

    bool fork(const std::string& cmd);

This function creates a new ``child`` process to execute the shell command passed in argument (on UNIX systems, using ``/bin/sh``). The child process runs in the background, while execution continues in the main process. Note that any child process will be automatically killed when the main process terminates: child processes cannot survive longer than the main process.

If you wish simply to run a command and wait for it to finish, see ``spawn()``.

**Example:**

.. code-block:: c++

    // Calculate some things to make an image
    vec2d img = /* ... */;

    // Write this image to the disk in FITS format
    fits::write("tmp.fits", img);

    // Call the FITS viewer DS9 on this data to let user inspect it.
    // The DS9 window will open in the background.
    fork("ds9 tmp.fits");

    // Execution comes back to us immediately, so the user can keep
    // inspecting the image while we do some more calculations.
    // Just be careful not to open too many windows in this way!
