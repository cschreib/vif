#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
    class worker {
        std::atomic<bool>                      stop_on_empty_;
        lock_free_queue<std::function<void()>> jobs_;
        std::unique_ptr<std::thread>           thread_;

        void consume_loop_() {
            while (!stop_on_empty_ || !jobs_.empty()) {
                std::function<void()> job;
                while (jobs_.pop(job)) {
                    job();
                }

                sleep_for(0.05);
            }
        }

    public :

        worker() : stop_on_empty_(false) {}

        ~worker() {
            wait();
        }

        void push(std::function<void()> job) {
            if (!thread_) {
                thread_ = std::unique_ptr<std::thread>(new std::thread(
                    &worker::consume_loop_, this
                ));
            }

            jobs_.push(std::move(job));
        }

        void wait() {
            if (thread_ && thread_->joinable()) {
                stop_on_empty_ = true;
                thread_->join();
                thread_ = nullptr;
                stop_on_empty_ = false;
            }
        }
    };
}
}
