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

    /// Thread-safe and lock-free FIFO queue.
    /// Single Producer, Single Consumer (SPSC).
    /** Note: implementation is from:
        http://www.drdobbs.com/parallel/writing-lock-free-code-a-corrected-queue/210604448
    **/
    template<typename T>
    class lock_free_queue {
        struct node {
            node() : next(nullptr) {}
            template<typename U>
            node(U&& t) : data(std::forward<U>(t)), next(nullptr) {}

            T     data;
            node* next;
        };

        // Modified and read by 'procuder' only
        node* first_;
        // Modified and read by 'consumer', read by 'producer'
        std::atomic<node*> dummy_;
        // Modified and read by 'procuder', read by 'consumer'
        std::atomic<node*> last_;

    public :
        lock_free_queue() {
            // Always keep a dummy separator between head and tail
            first_ = last_ = dummy_ = new node();
        }

        ~lock_free_queue() {
            // Clear the whole queue
            while (first_ != nullptr) {
                node* temp = first_;
                first_ = temp->next;
                delete temp;
            }
        }

        lock_free_queue(const lock_free_queue& q) = delete;
        lock_free_queue& operator = (const lock_free_queue& q) = delete;

        /// Push a new element at the back of the queue.
        /** Called by the 'producer' thread only.
        **/
        template<typename U>
        void push(U&& t) {
            // Add the new item to the queue
            (*last_).next = new node(std::forward<U>(t));
            last_ = (*last_).next;

            // Clear consumed items
            while (first_ != dummy_) {
                node* temp = first_;
                first_ = temp->next;
                delete temp;
            }
        }

        /// Pop an element from the front of the queue.
        /** Called by the 'consumer' thread only.
        **/
        bool pop(T& t) {
            // Return false if queue is empty
            if (dummy_ != last_) {
                // Consume the value
                t = std::move((*dummy_).next->data);
                dummy_ = (*dummy_).next;
                return true;
            } else {
                return false;
            }
        }

        /// Compute the current number of elements in the queue.
        /** Called by the 'consumer' thread only.
            Note that the true size may actually be larger than the returned value if the
            'producer' thread pushes new elements while the size is computed.
        **/
        std::size_t size() const {
            node* tmp = dummy_.load();
            node* end = last_.load();
            std::size_t n = 0;

            while (tmp != end) {
                tmp = tmp->next;
                ++n;
            }

            return n;
        }

        /// Check if this queue is empty.
        /** Called by the 'consumer' thread only.
        **/
        bool empty() const {
            return dummy_ == last_;
        }

        /// Delete all elements from the queue.
        /** This method should not be used in concurrent situations
        **/
        void clear() {
            while (first_ != dummy_) {
                node* temp = first_;
                first_ = temp->next;
                delete temp;
            }

            first_ = dummy_.load();

            while (first_->next != nullptr) {
                node* temp = first_;
                first_ = temp->next;
                delete temp;
            }

            last_ = dummy_ = first_;
        }
    };

    class worker {
        std::atomic<bool>                      stop_on_empty_;
        lock_free_queue<std::function<void()>> jobs_;
        std::unique_ptr<std::thread>           thread_;

        void consume_loop_() {
            while (!stop_on_empty_) {
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

#endif
