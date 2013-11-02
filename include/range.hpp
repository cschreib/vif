#ifndef RANGE_HPP
#define RANGE_HPP

#include <vec.hpp>

template<typename T>
struct range_t;

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
        return iter.i == i && iter.range == range;
    }

    bool operator != (const range_iterator_t& iter) const {
        return iter.i != i || &iter.range != &range;
    }
};

template<typename T>
struct range_t {
    T b, e, d;
    std::size_t n;

    range_t(T b_, T e_, std::size_t n_) : b(b_), e(e_), d((e_-b_)/n_), n(n_) {}

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
    return range(i, e, e-i);
}

template<typename T>
range_t<T> range(T n) {
    return range(T(0), n);
}

template<typename T, typename enable = typename std::enable_if<is_vec<T>::value>::type>
range_t<T> range(T n) {
    return range(n.size());
}

#endif
