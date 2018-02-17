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

   table.docutils code.docutils {
       background: initial;
       border: none;
   }

   table.docutils * {
       border: none;
   }

   .rst-content table.docutils td {
       border-bottom: none;
       white-space: normal;
   }
   </style>

The interface of phy++ was designed to facilitate the migration from IDL, an interpreted language that, wile dated, is still commonly used in astronomy. Although IDL as a language suffers from a number of design issues, there is much good in its interface and API that one may wish to emulate in C++. But not all of it.

This page lists common language constructs in IDL and their C++ equivalent with phy++. The following table is inspired from the `xtensor documentation <https://xtensor.readthedocs.io/en/latest/numpy.html>`_.

Main differences to keep in mind
--------------------------------

* C++ is a statically typed language: a variable, once created, can never change its type.
* C++ is a row-major language, IDL is a column-major language: for the same layout in memory (or on the disk), the dimensions of a vector in C++ are reversed compared to IDL. In particular, an image is accessed as ``img[x,y]`` in IDL, and ``img(y,x)`` in C++.
* C++ is case-sensitive, so ``a`` and ``A`` are different objects. Keywords and functions in C++ are always lowercase by convention.
* C++'s loops *much* faster than IDL's, so there is no need to avoid them.
* C++ has no procedures, but functions are allowed to have no return values.
* C++ does not support keywords for functions, only normal arguments. A structure with named member values can be used instead.
* phy++ vectors can be empty, while IDL vectors cannot. An empty vector has zero elements, and it is possible to do operations with an empty vector (which have no cost and do nothing) if all the other vectors involved, if any, are also empty.


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
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endif else begin``                         | | ``} else {``                                 |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endelse``                                  | | ``}``                                        |
+------------------------------------------------+------------------------------------------------+
| | ``for i=0, n-1 do begin``                    | | ``for (uint_t i : range(n)) {``              |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``break``                                    | | ``break;``                                   |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``continue``                                 | | ``continue;``                                |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endfor``                                   | | ``}``                                        |
+------------------------------------------------+------------------------------------------------+
| | ``array = ['foo','bar','blob']``             | | ``vec1s array = {"foo","bar","blob"};``      |
| | ``foreach val, array do begin``              | | ``for (std::string val : array) {``          |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endforeach``                               | | ``}``                                        |
+------------------------------------------------+------------------------------------------------+
| | ``while a gt b do begin``                    | | ``while (a > b) {``                          |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endfor``                                   | | ``}``                                        |
+------------------------------------------------+------------------------------------------------+
| | ``repeat begin``                             | | ``do {``                                     |
| | ``// ...``                                   | | ``// ...``                                   |
| | ``endrep until a gt b``                      | | ``} while (a > b);``                         |
+------------------------------------------------+------------------------------------------------+
| | ``switch i of``                              | | ``switch (i) {``                             |
| | ``1: print, 'one'``                          | | ``case 1: print("one");``                    |
| | ``2: print, 'two'``                          | | ``case 2: print("two");``                    |
| | ``3: print, 'three'``                        | | ``case 3: print("three");``                  |
| | ``4: begin``                                 | | ``case 4:``                                  |
| | ``    print, 'four'``                        | | ``    print("four");                         |
| | ``    break``                                | | ``    break;``                               |
| | ``   end``                                   | |                                              |
| | ``else: print, 'other'``                     | | ``default: print("other");``                 |
| | ``endswitch``                                | | ``}``                                        |
| |                                              | | Note: only works with integers, no strings.  |
+------------------------------------------------+------------------------------------------------+
| | ``case i of``                                | | No direct equivalent. Use ``switch()`` and   |
| | ``// ...``                                   | | be sure to call ``break;`` at the end of     |
| | ``endcase``                                  | | each case.                                   |
+------------------------------------------------+------------------------------------------------+


Finding values
--------------

+------------------------------------------------+------------------------------------------------+
|             IDL                                |               C++ 11 - phy++                   |
+================================================+================================================+
| | ``v = [1,2,3,4,5]``                          | | ``vec1f v = {1,2,3,4,5};``                   |
| | ``id = where(v gt 3, cnt)``                  | | ``vec1u id = where(v > 3);``                 |
| | ``if cnt ne 0 then v[id] = 0``               | | ``v[id] = 0;``                               |
| |                                              | | Note: empty vectors are allowed in phy++,    |
| |                                              | | so the check for ``cnt`` is not needed.      |
+------------------------------------------------+------------------------------------------------+
