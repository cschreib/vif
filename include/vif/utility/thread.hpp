#ifndef VIF_UTILITY_THREAD_HPP
#define VIF_UTILITY_THREAD_HPP

#include <thread>
#include <atomic>
#include <mutex>
#include <functional>
#include "vif/core/vec.hpp"

#define VIF_INCLUDING_THREAD_BITS
#include "vif/core/bits/thread-thread.hpp"
#include "vif/core/bits/thread-utils.hpp"
#include "vif/core/bits/thread-queue.hpp"
#include "vif/core/bits/thread-worker.hpp"
#include "vif/core/bits/thread-worker-pool.hpp"
#undef VIF_INCLUDING_THREAD_BITS

#endif
