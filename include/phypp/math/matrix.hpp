#ifndef PHYPP_MATH_MATRIX_HPP
#define PHYPP_MATH_MATRIX_HPP

#ifndef NO_LAPACK
#include "phypp/math/lapack.hpp"
#endif
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/range.hpp"
#include "phypp/math/base.hpp"

#define PHYPP_INCLUDING_MATH_MATRIX_BITS
#include "phypp/math/bits/matrix_types.hpp"
#include "phypp/math/bits/matrix_functions.hpp"
#undef PHYPP_INCLUDING_MATH_MATRIX_BITS

#endif
