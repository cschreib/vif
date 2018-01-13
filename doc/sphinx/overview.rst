Overview of the core library
============================

At the core of the phy++ library is the *vector* class, ``vec``. This is basically an enhanced ``std::vector``, and it therefore shares most of its features and strengths. On top of the ``std::vector`` interface, the phy++ vectors have extra functionalities to simplify data analysis and calculations, including overloaded mathematical operators, multi-dimensional indexing, and the ability to create "views" to access and edit subsets of a given vector.

Here we will first describe the properties of the vector class, and then describe the vector views. Lastly, a guide for writing "generic" functions that work with any vector type is provided.
