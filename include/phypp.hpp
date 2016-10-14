// Core headers, containing essential features from the library
#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/print.hpp"
#include "phypp/core/error.hpp"

// New entry point of the program, necessary to get stack trace in case of errors
#include "phypp/core/main.hpp"

// Generic utility functions
#include "phypp/utility/string.hpp"
#include "phypp/utility/argv.hpp"
#include "phypp/utility/time.hpp"
#include "phypp/utility/thread.hpp"
#include "phypp/utility/generic.hpp"

// Reflection tools
#include "phypp/reflex/reflex.hpp"
#include "phypp/reflex/reflex_helpers.hpp"

// Input/output (IO) functions for reading and writing files
#include "phypp/io/filesystem.hpp"
#include "phypp/io/ascii.hpp"
#include "phypp/io/fits.hpp"

// Mathematics functions
#include "phypp/math/base.hpp"
#include "phypp/math/interpolate.hpp"
#include "phypp/math/reduce.hpp"
#include "phypp/math/histogram.hpp"
#include "phypp/math/random.hpp"
#include "phypp/math/matrix.hpp"
#include "phypp/math/linfit.hpp"
#include "phypp/math/convex_hull.hpp"
#include "phypp/math/complex.hpp"
#include "phypp/math/fourier.hpp"
#include "phypp/math/transform.hpp"

// Astronomy functions
#include "phypp/astro/image.hpp"
#include "phypp/astro/astro.hpp"
#include "phypp/astro/wcs.hpp"

// Bring common functions into the global namespace
using namespace phypp;
using namespace phypp::astro;

