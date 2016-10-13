#ifndef PHYPP_CORE_RANGE_HPP
#define PHYPP_CORE_RANGE_HPP

#include "phypp/core/error.hpp"
#include "phypp/core/typedefs.hpp"

template<typename T>
struct range_t;

template<typename T>
struct has_size_t {
    template <typename U> static std::true_type dummy(typename std::decay<
        decltype(std::declval<U&>().size())>::type*);
    template <typename U> static std::false_type dummy(...);
    using type = decltype(dummy<T>(0));
};

template<typename T>
using has_size = typename has_size_t<T>::type;

template<typename T>
struct range_iterator_t;

template<typename T>
struct range_iterator_t<range_t<T>> {
    const range_t<T>& range;
    std::size_t i;

    range_iterator_t& operator++ (int) {
        ++i;
        return *this;
    }

    range_iterator_t operator++ () {
        return range_iterator_t{range, i++};
    }

    T operator * () const {
        return range.b + range.d*i;
    }

    bool operator == (const range_iterator_t& iter) const {
        return iter.i == i && &iter.range == &range;
    }

    bool operator != (const range_iterator_t& iter) const {
        return iter.i != i || &iter.range != &range;
    }
};

template<typename T>
using range_dtype = typename std::conditional<
    std::is_unsigned<T>::value,
    typename std::make_signed<T>::type,
    T
>::type;

template<typename T>
struct range_t {
    T b, e;
    range_dtype<T> d;
    std::size_t n;

    range_t(T b_, T e_, std::size_t n_) : b(b_), e(e_),
        d(n_ == 0 ? 0 : (range_dtype<T>(e_)-range_dtype<T>(b_))/range_dtype<T>(n_)), n(n_) {}

    using iterator = range_iterator_t<range_t>;

    iterator begin() const { return iterator{*this, 0}; }
    iterator end() const { return iterator{*this, n}; }
};

template<typename T, typename U = T>
range_t<T> range(T i, U e, std::size_t n) {
    return range_t<T>(i, e, n);
}

template<typename T, typename U = T>
range_t<T> range(T i, U e) {
    return range(i, e, std::abs(range_dtype<T>(e)-range_dtype<T>(i)));
}

template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
range_t<T> range(T n) {
    return range(T(0), n);
}

template<typename T, typename enable = typename std::enable_if<has_size<T>::value>::type>
range_t<std::size_t> range(const T& n) {
    return range(n.size());
}

// Full range variable v(_)
static struct full_range_t {} _;

using placeholder_t = full_range_t;

// Left range variable v(_-last)
struct left_range_t {
    uint_t last;
};

// Right range variable v(first-_)
struct right_range_t {
    uint_t first;
};

// Left-right range variable v(first-_-last)
struct left_right_range_t {
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

namespace range_impl {
    template<typename T>
    void check_lower_bounds(const T& right, uint_t size) {
        phypp_check(right.first < size, "lower bound of range goes past the "
            "size of this vector (", right.first, " vs. ", size, ")");
    }

    template<typename T>
    void check_upper_bounds(const T& left, uint_t size) {
        phypp_check(left.last < size, "upper bound of range goes past the "
            "size of this vector (", left.last, " vs. ", size, ")");
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

    inline std::size_t range_size(full_range_t, std::size_t size) {
        return size;
    }
    inline std::size_t range_size(left_range_t rng, std::size_t size) {
        return rng.last+1;
    }
    inline std::size_t range_size(right_range_t rng, std::size_t size) {
        return size - rng.first;
    }
    inline std::size_t range_size(left_right_range_t rng, std::size_t size) {
        return (rng.last+1) - rng.first;
    }
}

inline range_t<std::size_t> range(full_range_t, std::size_t size) {
    return range(size);
}

inline range_t<std::size_t> range(const left_range_t& rng, std::size_t size) {
    range_impl::check_bounds(rng, size);
    return range(rng.last+1);
}

inline range_t<std::size_t> range(const right_range_t& rng, std::size_t size) {
    range_impl::check_bounds(rng, size);
    return range(rng.first, size);
}

inline range_t<std::size_t> range(const left_right_range_t& rng, std::size_t size) {
    range_impl::check_bounds(rng, size);
    return range(rng.first, rng.last+1);
}

template<typename T>
struct is_range : std::integral_constant<bool,
    std::is_same<T,full_range_t>::value || std::is_same<T,left_range_t>::value ||
    std::is_same<T,right_range_t>::value || std::is_same<T,left_right_range_t>::value> {};

#endif
