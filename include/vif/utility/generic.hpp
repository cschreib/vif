#ifndef VIF_UTILITY_GENERIC_HPP
#define VIF_UTILITY_GENERIC_HPP

#include <numeric>

#include "vif/core/vec.hpp"
#include "vif/core/meta.hpp"
#include "vif/core/range.hpp"
#include "vif/core/error.hpp"
#include "vif/core/string_conversion.hpp"

#define VIF_INCLUDING_GENERIC_BITS
#include "vif/utility/bits/generic-sequences.hpp"
#include "vif/utility/bits/generic-indices.hpp"
#include "vif/utility/bits/generic-dims.hpp"
#include "vif/utility/bits/generic-find.hpp"
#include "vif/utility/bits/generic-rearrange.hpp"
#undef VIF_INCLUDING_GENERIC_BITS

#endif
