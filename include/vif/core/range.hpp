#ifndef VIF_CORE_RANGE_HPP
#define VIF_CORE_RANGE_HPP

#include <utility>
#include "vif/core/error.hpp"
#include "vif/core/typedefs.hpp"

namespace vif {
namespace impl {
    namespace range_impl {
        template<typename T>
        struct range_t;

        template<typename T>
        struct has_size_t {
            template <typename U> static std::true_type dummy(typename std::decay<
                decltype(std::declval<U&>().size())>::type*);
            template <typename U> static std::false_type dummy(...);
            using type = decltype(dummy<T>(nullptr));
        };

        template<typename T>
        using has_size = typename has_size_t<T>::type;

        template<typename T, bool integral>
        struct dtype_impl;

        template <typename T>
        struct dtype_impl<T,true> {
            using type = typename std::make_signed<T>::type;
        };

        template <typename T>
        struct dtype_impl<T,false> {
            using type = T;
        };

        template<typename T>
        using dtype = typename dtype_impl<T, std::is_integral<T>::value>::type;

        template <typename T>
        dtype<T> make_step(dtype<T> e, dtype<T> b, dtype<T> n) {
            return (e-b)/n;
        }

        template<typename T>
        struct iterator_t;

        template<typename T>
        struct iterator_t<range_t<T>> {
            const range_t<T>& range;
            uint_t i;

            iterator_t& operator++ (int) {
                ++i;
                return *this;
            }

            iterator_t operator++ () {
                return iterator_t{range, i++};
            }

            T operator * () const {
                return range.b + range.d*i;
            }

            bool operator == (const iterator_t& iter) const {
                return iter.i == i && &iter.range == &range;
            }

            bool operator != (const iterator_t& iter) const {
                return iter.i != i || &iter.range != &range;
            }
        };

        template<typename T>
        struct range_t {
            T b, e;
            dtype<T> d;
            uint_t n;

            explicit range_t(T b_, T e_, uint_t n_) : b(b_), e(e_),
                d(n_ == 0 ? 0 : make_step<T>(e_, b_, n_)), n(n_) {}

            using iterator = iterator_t<range_t>;

            iterator begin() const { return iterator{*this, 0}; }
            iterator end() const { return iterator{*this, n}; }
        };

        // Full range v(_)
        struct full_range_t {};

        // Left range variable v(_-last)
        struct left_range_t {
            explicit left_range_t(uint_t l) : last(l) {}
            uint_t last;
        };

        // Right range variable v(first-_)
        struct right_range_t {
            explicit right_range_t(uint_t f) : first(f) {}
            uint_t first;
        };

        // Left-right range variable v(first-_-last)
        struct left_right_range_t {
            explicit left_right_range_t(uint_t f, uint_t l) : first(f), last(l) {}
            uint_t first, last;
        };

        template<typename T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
        left_range_t operator - (full_range_t, T last) {
            return left_range_t{static_cast<uint_t>(last)};
        }

        template<typename T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
        left_right_range_t operator - (T first, left_range_t left) {
            return left_right_range_t{static_cast<uint_t>(first), left.last};
        }

        template<typename T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
        right_range_t operator - (T first, full_range_t) {
            return right_range_t{static_cast<uint_t>(first)};
        }

        template<typename T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
        left_right_range_t operator - (right_range_t right, T last) {
            return left_right_range_t{right.first, static_cast<uint_t>(last)};
        }

        template<typename T>
        void check_lower_bounds(const T& right, uint_t size) {
            vif_check(right.first < size, "lower bound of range goes past the "
                "size of the data (", right.first, " vs. ", size, ")");
        }

        template<typename T>
        void check_upper_bounds(const T& left, uint_t size) {
            vif_check(left.last < size, "upper bound of range goes past the "
                "size of the data (", left.last, " vs. ", size, ")");
        }

        inline void check_bounds(full_range_t, uint_t size) {}

        inline void check_bounds(const left_range_t& left, uint_t size) {
            check_upper_bounds(left, size);
        }

        inline void check_bounds(const right_range_t& right, uint_t size) {
            check_lower_bounds(right, size);
        }

        inline void check_bounds(const left_right_range_t& rng, uint_t size) {
            check_upper_bounds(rng, size);
            check_lower_bounds(rng, size);
        }

        inline uint_t range_size(full_range_t, uint_t size) {
            return size;
        }
        inline uint_t range_size(left_range_t rng, uint_t size) {
            return rng.last+1;
        }
        inline uint_t range_size(right_range_t rng, uint_t size) {
            return size - rng.first;
        }
        inline uint_t range_size(left_right_range_t rng, uint_t size) {
            return (rng.last+1) - rng.first;
        }
    }

    using placeholder_t = impl::range_impl::full_range_t;
}

    // Full range variable v(_)
    static impl::range_impl::full_range_t _;

    // Define a range from start 'i', end 'e' (exclusive), and number of steps 'n'
    template<typename T, typename U = T>
    impl::range_impl::range_t<T> range(T i, U e, uint_t n) {
        return impl::range_impl::range_t<T>(i, e, n);
    }

    // Define a range from start 'i', end 'e' (exclusive) in integer steps
    template<typename T, typename U = T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
    impl::range_impl::range_t<T> range(T i, U e) {
        return range(i, e, std::abs(impl::range_impl::dtype<T>(e)-impl::range_impl::dtype<T>(i)));
    }

    // Define a range from '0' to 'n' (exclusive) in integer steps
    template<typename T, typename enable = typename std::enable_if<std::is_integral<T>::value>::type>
    impl::range_impl::range_t<T> range(T n) {
        return range(T(0), n);
    }

    // Define a range from '0' to 'n.size()' (exclusive) in integer steps
    template<typename T, typename enable = typename std::enable_if<impl::range_impl::has_size<T>::value>::type>
    auto range(const T& n) -> impl::range_impl::range_t<decltype(n.size())> {
        return range(n.size());
    }

    inline uint_t range_begin(impl::range_impl::full_range_t, uint_t size) {
        return 0u;
    }

    inline uint_t range_begin(impl::range_impl::left_range_t, uint_t size) {
        return 0u;
    }

    inline uint_t range_begin(impl::range_impl::right_range_t r, uint_t size) {
        return r.first;
    }

    inline uint_t range_begin(impl::range_impl::left_right_range_t r, uint_t size) {
        return r.first;
    }

    inline uint_t range_end(impl::range_impl::full_range_t, uint_t size) {
        return size;
    }

    inline uint_t range_end(impl::range_impl::left_range_t r, uint_t size) {
        return r.last+1;
    }

    inline uint_t range_end(impl::range_impl::right_range_t, uint_t size) {
        return size;
    }

    inline uint_t range_end(impl::range_impl::left_right_range_t r, uint_t size) {
        return r.last+1;
    }

namespace meta {
    template<typename T>
    struct is_range : std::integral_constant<bool,
        std::is_same<T,impl::range_impl::full_range_t>::value ||
        std::is_same<T,impl::range_impl::left_range_t>::value ||
        std::is_same<T,impl::range_impl::right_range_t>::value ||
        std::is_same<T,impl::range_impl::left_right_range_t>::value> {};
}
}

#endif
