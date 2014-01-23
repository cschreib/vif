#ifndef THREAD_HPP
#define THREAD_HPP

#include "phypp/vec.hpp"
#include <thread>
#include <atomic>

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

    using pool_t = std::vector<thread_t>;

    pool_t pool(uint_t n) {
        return pool_t(n);
    }

    void sleep_for(double duration) {
        std::this_thread::sleep_for(std::chrono::microseconds(uint_t(duration*1e6)));
    }
}

#endif
