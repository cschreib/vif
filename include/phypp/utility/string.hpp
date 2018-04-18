#ifndef PHYPP_UTILITY_STRING_HPP
#define PHYPP_UTILITY_STRING_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/string_conversion.hpp"
#include "phypp/utility/generic.hpp"

#define PHYPP_INCLUDING_STRING_BITS
#include "phypp/utility/bits/string-base.hpp"
#include "phypp/utility/bits/string-find_replace.hpp"
#include "phypp/utility/bits/string-split.hpp"
#include "phypp/utility/bits/string-format.hpp"
#include "phypp/utility/bits/string-regex.hpp"
#include "phypp/utility/bits/string-hash.hpp"
#undef PHYPP_INCLUDING_STRING_BITS

#endif

