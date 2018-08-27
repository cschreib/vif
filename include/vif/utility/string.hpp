#ifndef VIF_UTILITY_STRING_HPP
#define VIF_UTILITY_STRING_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/core/string_conversion.hpp"
#include "vif/utility/generic.hpp"

#define VIF_INCLUDING_STRING_BITS
#include "vif/utility/bits/string-base.hpp"
#include "vif/utility/bits/string-find_replace.hpp"
#include "vif/utility/bits/string-split.hpp"
#include "vif/utility/bits/string-format.hpp"
#include "vif/utility/bits/string-regex.hpp"
#include "vif/utility/bits/string-hash.hpp"
#undef VIF_INCLUDING_STRING_BITS

#endif

