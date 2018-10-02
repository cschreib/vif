#ifndef VIF_UTILITY_THREAD_HPP
#define VIF_UTILITY_THREAD_HPP

#include <thread>
#include <atomic>
#include <mutex>
#include <functional>
#include "vif/core/vec.hpp"

#define VIF_INCLUDING_THREAD_BITS
#include "vif/utility/bits/thread-thread.hpp"
#include "vif/utility/bits/thread-utils.hpp"
#include "vif/utility/bits/thread-queue.hpp"
#include "vif/utility/bits/thread-worker.hpp"
#include "vif/utility/bits/thread-worker-pool.hpp"
#include "vif/utility/bits/thread-parallel-for.hpp"
#undef VIF_INCLUDING_THREAD_BITS

#endif
