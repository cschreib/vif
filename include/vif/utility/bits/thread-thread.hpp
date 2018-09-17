#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
    struct thread_t {
        std::unique_ptr<std::thread> impl;

        template<typename F, typename ... Args>
        void start(F&& f, Args&& ... args) {
            impl = std::unique_ptr<std::thread>(new std::thread(
                std::forward<F>(f), std::forward<Args>(args)...
            ));
        }

        void join() {
            impl->join();
        }
    };
}
}
