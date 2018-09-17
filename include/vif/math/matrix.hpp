#ifndef VIF_MATH_MATRIX_HPP
#define VIF_MATH_MATRIX_HPP

#ifndef NO_LAPACK
#include "vif/math/lapack.hpp"
#endif
#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/math/base.hpp"

#define VIF_INCLUDING_MATH_MATRIX_BITS
#include "vif/math/bits/matrix-types.hpp"
#include "vif/math/bits/matrix-functions.hpp"
#undef VIF_INCLUDING_MATH_MATRIX_BITS

#endif
