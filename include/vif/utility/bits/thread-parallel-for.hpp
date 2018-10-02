#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
    struct parallel_for {
        std::vector<thread_t> workers;
        bool verbose = false;
        uint_t progress_step = 1;
        double update_rate = 0.1;

        parallel_for() = default;

        explicit parallel_for(uint_t nthread) {
            workers.resize(nthread);
        }

        template<typename F>
        void execute(const F& f, uint_t ifirst, uint_t ilast) {
            uint_t n = ilast - ifirst;

            if (workers.empty()) {
                // Single-threaded execution
                auto pg = progress_start(n);
                for (uint_t i : range(ifirst, ilast)) {
                    f(i);
                    if (verbose) {
                        progress(pg, progress_step);
                    }
                }
            } else {
                // Multi-threaded execution
                std::atomic<uint_t> iter(0);
                uint_t di = n/workers.size() + 1;
                uint_t i0 = ifirst;
                uint_t i1 = ifirst + di;
                bool local_verbose = verbose;
                for (auto& t : workers) {
                    t.start([&f,&iter,i0,i1,local_verbose]() {
                        for (uint_t i : range(i0, i1)) {
                            f(i);
                            if (local_verbose) {
                                ++iter;
                            }
                        }
                    });

                    i0 += di;
                    i1 += di;
                    if (i1 > n) {
                        i1 = n;
                    }
                }

                // Wait for all threads to finish before returning
                if (verbose) {
                    auto pg = progress_start(n);
                    while (iter < n) {
                        thread::sleep_for(update_rate);
                        if (verbose) print_progress(pg, iter);
                    }
                }

                for (auto& t : workers) {
                    t.join();
                }
            }
        }

        template<typename F>
        void execute(const F& f, uint_t i1) {
            execute(f, 0, i1);
        }

        uint_t size() const {
            return workers.size();
        }
    };
}
}
