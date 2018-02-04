Error checking
==============

In ``phypp/core/error.hpp``.

phypp_check
-----------

.. code-block:: c++

    template<typename ... Args>
    void phypp_check(bool b, Args&& ... args);
