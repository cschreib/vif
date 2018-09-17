#ifndef VIF_INCLUDING_THREAD_BITS
#error this file is not meant to be included separately, include "vif/utilty/thread.hpp" instead
#endif

namespace vif {
namespace thread {
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

        // Modified and read by 'producer' only
        node* first_;
        // Modified and read by 'consumer', read by 'producer'
        std::atomic<node*> dummy_;
        // Modified and read by 'producer', read by 'consumer'
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
}
}
