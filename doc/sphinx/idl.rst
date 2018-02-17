IDL equivalents
===============

.. raw:: html

   <style>
   table.docutils {
       width: 100%;
       table-layout: fixed;
       border: none;
   }

   table.docutils th {
       text-align: center;
   }

   table.docutils .line-block {
       margin-left: 0;
       margin-bottom: 0;
   }

   table.docutils code.literal {
       color: initial;
   }

   code.docutils {
       background: initial;
       border: none;
   }

   * {
       border: none;
   }

   .rst-content table.docutils td {
       border-bottom: none;
       border-left: none;
   }
   </style>

The interface of phy++ was designed to facilitate the migration from IDL, an interpreted language that, wile dated, is still commonly used in astronomy. Although IDL as a language suffers from a number of design issues, there is much good in its interface and API that one may wish to emulate in C++. But not all of it.

This page lists common language constructs in IDL and their C++ equivalent with phy++. The following table is inspired from the `xtensor documentation <https://xtensor.readthedocs.io/en/latest/numpy.html>`_.

Main differences to keep in mind
--------------------------------

* C++ is a statically typed language: a variable, once created, can never change its type.
* C++ is a row-major language, IDL is a column-major language: for the same layout in memory (or on the disk), the dimensions of a vector in C++ are reversed compared to IDL. In particular, an image is accessed as ``img[x,y]`` in IDL, and ``img(y,x)`` in C++.


Creating, accessing, modifying vectors
--------------------------------------

+------------------------------------------------+------------------------------------------------+
|             IDL                                |               C++ 11 - phy++                   |
+================================================+================================================+
| | ``v = fltarr(10)``                           | | ``vec1f v(10);``                             |
| | ``v = fltarr(20)``                           | | ``v.resize(20);``                            |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(5)``                            | | ``vec1i v(10);``                             |
| | ``d = double(v)``                            | | ``vec1d d = v;``                             |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(5)``                            | No equivalent. Types in C++ are *static*,      |
| | ``v = double(v)``                            | cannot change ``int`` to ``double``.           |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(6)``                            | | ``vec1i v(6);``                              |
| | ``d = reform(v, 3, 2)``                      | | ``vec2i d = reform(v, 2, 3);``               |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(2, 3)``                         | | ``vec2i v(3, 2);``                           |
| | ``v = reform(v, 3, 2)``                      | | ``v = reform(v, 2, 3);``                     |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(6)``                            | No equivalent. The number of dimensions of a   |
| | ``v = reform(v, 3, 2)``                      | vector is part of its type, and cannot change. |
+------------------------------------------------+------------------------------------------------+
| | ``v = [1,2,5,7]``                            | | ``vec1i v = {1,2,5,7};``                     |
| | ``v = [1,2,3]``                              | | ``v = {1,2,3};``                             |
+------------------------------------------------+------------------------------------------------+
| ``n_elements(v)``                              | ``v.size();``                                  |
+------------------------------------------------+------------------------------------------------+
| ``v = dindgen(5)``                             | ``vec1d v = dindgen(5);``                      |
+------------------------------------------------+------------------------------------------------+
| | ``v = indgen(2,3)``                          | | ``vec2i v = indgen(3,2);``                   |
| | ``v[0] = 1``                                 | | ``v[0] = 1;``                                |
| | ``v[0,2] = 2``                               | | ``v(2,0) = 2;``                              |
| | ``v[0,*] = [2,5,6]``                         | | ``v(_,0) = {2,5,6};``                        |
| | ``v[0,*:1] = [5,6]``                         | | ``v(_-1,0) = {5,6};``                        |
| | ``v[0,1:*] = [5,6]``                         | | ``v(1-_,0) = {5,6};``                        |
| | ``v[0,1:2] = [5,6]``                         | | ``v(1-_-2,0) = {5,6};``                      |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(5)``                            | | ``vec1i v(5);``                              |
| | ``w = intarr(5)``                            | | ``vec1i w(5);``                              |
| | ``id = [1,3,4]``                             | | ``vec1u id = {1,3,4};``                      |
| | ``v[id] = 1``                                | | ``v[id] = 1;``                               |
| | ``v[id] = [-1,0,1]``                         | | ``v[id] = {-1,0,1};``                        |
| | ``w[id] = v[id]``                            | | ``w[id] = v[id];``                           |
+------------------------------------------------+------------------------------------------------+
| | ``v = intarr(5)``                            | | ``vec1i v(5);``                              |
| | ``v[0] = [1,2]`` (optimized assignment)      | | ``v[0-_-1] = {1,2};`` (need explicit range)  |
+------------------------------------------------+------------------------------------------------+


Vector operations
-----------------

+------------------------------------------------+------------------------------------------------+
|             IDL                                |               C++ 11 - phy++                   |
+================================================+================================================+
| | ``v = indgen(5)``                            | | ``vec1i v = indgen(5);``                     |
| | ``w = indgen(5)``                            | | ``vec1i w = indgen(5);``                     |
| |                                              | | ``vec1i x;``                                 |
| | ``x = v + w``                                | | ``x = v + w;``                               |
| | ``x = v - w``                                | | ``x = v - w;``                               |
| | ``x = v * w``                                | | ``x = v * w;``                               |
| | ``x = v / w``                                | | ``x = v / w;``                               |
| | ``x = v ^ w``                                | | ``x = pow(v, w);``                           |
| | ``x = v mod w``                              | | ``x = v % w;``                               |
| | ``x = v gt w``                               | | ``x = v > w;``                               |
| | ``x = v ge w``                               | | ``x = v >= w;``                              |
| | ``x = v lt w``                               | | ``x = v < w;``                               |
| | ``x = v le w``                               | | ``x = v <= w;``                              |
| | ``x = v and w``                              | | ``x = v && w;``                              |
| | ``x = v or w``                               | | ``x = v || w;``                              |
| | ``x = v > w``                                | | ``x = max(v, w);``                           |
| | ``x = v < w``                                | | ``x = min(v, w);``                           |
+------------------------------------------------+------------------------------------------------+
| | ``v = indgen(5,5)+1``                        | | ``vec2i v = indgen(5,5)+1;``                 |
| | ``w = indgen(5,5)+0``                        | | ``vec2i w = indgen(5,5)+0;``                 |
| | ``x = v # w``                                | | ``vec2i x = matrix::product(w, v);``         |
+------------------------------------------------+------------------------------------------------+
| ``x = v ## w``                                 | No direct equivalent. Do the operation         |
|                                                | explicitly with indices in a loop.             |
+------------------------------------------------+------------------------------------------------+


Control flow
------------

+------------------------------------------------+------------------------------------------------+
|             IDL                                |               C++ 11 - phy++                   |
+================================================+================================================+
| | ``if x lt y then begin``                     | | ``if (x < y) {``                             |
| | ...                                          | | ...                                          |
| | ``endif else begin``                         | | ``} else {``                                 |
| | ...                                          | | ...                                          |
| | ``endelse``                                  | | ``}``                                        |
+------------------------------------------------+------------------------------------------------+


Finding values
--------------

+------------------------------------------------+------------------------------------------------+
|             IDL                                |               C++ 11 - phy++                   |
+================================================+================================================+
| | ``v = [1,2,3,4,5]``                          | | ``vec1f v = {1,2,3,4,5};``                   |
| | ``id = where(v gt 3, cnt)``                  | | ``vec1u id = where(v > 3);``                 |
| | ``if cnt ne 0 then v[id] = 0``               | | ``v[id] = 0;``                               |
+------------------------------------------------+------------------------------------------------+
