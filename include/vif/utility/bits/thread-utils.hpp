#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
    using pool_t = std::vector<thread_t>;

    inline pool_t pool(uint_t n) {
        return pool_t(n);
    }

    inline void sleep_for(double duration) {
        std::this_thread::sleep_for(std::chrono::microseconds(uint_t(duration*1e6)));
    }
}
}
