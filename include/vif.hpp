// Core headers, containing essential features from the library
#include "vif/core/vec.hpp"
#include "vif/core/range.hpp"
#include "vif/core/print.hpp"
#include "vif/core/error.hpp"

// New entry point of the program, necessary to get stack trace in case of errors
#include "vif/core/main.hpp"

// Generic utility functions
#include "vif/utility/string.hpp"
#include "vif/utility/argv.hpp"
#include "vif/utility/os.hpp"
#include "vif/utility/time.hpp"
#include "vif/utility/thread.hpp"
#include "vif/utility/generic.hpp"

// Reflection tools
#include "vif/reflex/reflex.hpp"
#include "vif/reflex/reflex_helpers.hpp"

// Input/output (IO) functions for reading and writing files
#include "vif/io/filesystem.hpp"
#include "vif/io/ascii.hpp"
#include "vif/io/fits.hpp"

// Mathematics functions
#include "vif/math/base.hpp"
#include "vif/math/interpolate.hpp"
#include "vif/math/reduce.hpp"
#include "vif/math/histogram.hpp"
#include "vif/math/random.hpp"
#include "vif/math/matrix.hpp"
#include "vif/math/linfit.hpp"
#include "vif/math/convex_hull.hpp"
#include "vif/math/complex.hpp"
#include "vif/math/fourier.hpp"
#include "vif/math/transform.hpp"

// Astronomy functions
#include "vif/astro/image.hpp"
#include "vif/astro/astro.hpp"
#include "vif/astro/wcs.hpp"

// Bring common functions into the global namespace
using namespace vif;
using namespace vif::astro;

