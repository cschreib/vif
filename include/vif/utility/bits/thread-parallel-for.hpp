#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
    struct parallel_for {
        // Setup
        bool verbose = false;
        uint_t progress_step = 1;
        double update_rate = 0.1;
        uint_t chunk_size = 0;

    private :

        struct worker_t {
            thread_t thread;
            parallel_for& pfor;

            explicit worker_t(parallel_for& p) : pfor(p) {}

            worker_t(worker_t&& w) : pfor(w.pfor) {}

            worker_t(const worker_t& w) = delete;

            template<typename F>
            void start(const F& f, bool verbose) {
                thread.start([this,&f,verbose]() {
                    uint_t i0, i1;
                    while (pfor.query_chunk(i0, i1)) {
                        for (uint_t i : range(i0, i1)) {
                            f(i);
                            if (verbose) {
                                ++pfor.iter;
                            }
                        }
                    }
                });
            }

            void join() {
                thread.join();
            }
        };

        // Workers
        std::vector<parallel_for::worker_t> workers;

        // Internal
        std::mutex query_mutex;
        std::atomic<uint_t> iter;
        uint_t n, i0, i1, di;

        bool query_chunk(uint_t& oi0, uint_t& oi1) {
            std::unique_lock<std::mutex> l(query_mutex);

            if (i0 == n) {
                return false;
            }

            oi0 = i0;
            oi1 = i1;

            i0 = i1;
            i1 += di;
            if (i1 > n) {
                i1 = n;
            }

            return true;
        }

    public :

        parallel_for() = default;
        parallel_for(const parallel_for&) = delete;
        parallel_for(parallel_for&&) = delete;

        explicit parallel_for(uint_t nthread) {
            for (uint_t i : range(nthread)) {
                workers.emplace_back(*this);
            }
        }

        template<typename F>
        void execute(const F& f, uint_t ifirst, uint_t ilast) {
            n = ilast - ifirst;

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
                uint_t nchunk = (chunk_size == 0 ?
                    workers.size() : std::max(workers.size(), n/chunk_size));

                // Setup chunks
                di = n/nchunk + 1;
                i0 = ifirst;
                i1 = ifirst + di;
                iter = 0;

                // Launch threads
                for (auto& t : workers) {
                    t.start(f, verbose);
                }

                progress_t pg;
                if (verbose) {
                    pg = progress_start(n);
                }

                // Wait for all threads to finish before returning
                if (verbose) {
                    while (iter < n) {
                        thread::sleep_for(update_rate);
                        print_progress(pg, iter);
                    }
                }

                // Join all
                for (auto& t : workers) {
                    t.join();
                }

                if (verbose) {
                    print_progress(pg, iter);
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
