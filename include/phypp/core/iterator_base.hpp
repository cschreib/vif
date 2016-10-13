#ifndef PHYPP_CORE_ITERATOR_BASE_HPP
#define PHYPP_CORE_ITERATOR_BASE_HPP

template<typename T, typename C, typename P>
class iterator_base;

template<typename T, typename C, typename P>
class const_iterator_base {
    friend C;
    template<typename N, typename K, typename V>
    friend class const_reverse_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    const_iterator_base(const T& j) : i(j) {}

public :
    const_iterator_base() {}
    const_iterator_base(const iterator_base<T,C,P>& iter) : i(iter.i) {}

    const_iterator_base operator ++ (int) {
        const_iterator_base iter = *this;
        ++i;
        return iter;
    }

    const_iterator_base& operator ++ () {
        ++i; return *this;
    }

    const_iterator_base operator -- (int) {
        const_iterator_base iter = *this;
        --i;
        return iter;
    }

    const_iterator_base& operator -- () {
        --i; return *this;
    }

    void operator += (int_t n) {
        i += n;
    }

    void operator -= (int_t n) {
        i -= n;
    }

    const_iterator_base operator + (int_t n) const {
        const_iterator_base iter = *this;
        iter.i += n;
        return iter;
    }

    const_iterator_base operator - (int_t n) const {
        const_iterator_base iter = *this;
        iter.i -= n;
        return iter;
    }

    std::ptrdiff_t operator - (const const_iterator_base& iter) const {
        return i - iter.i;
    }

    bool operator == (const const_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const const_iterator_base& iter) const {
        return i != iter.i;
    }

    bool operator < (const const_iterator_base& iter) const {
        return i < iter.i;
    }

    bool operator <= (const const_iterator_base& iter) const {
        return i <= iter.i;
    }

    bool operator > (const const_iterator_base& iter) const {
        return i > iter.i;
    }

    bool operator >= (const const_iterator_base& iter) const {
        return i >= iter.i;
    }

    auto operator * () -> decltype(P::get_obj(i)) {
        return P::get_obj(i);
    }

    auto operator -> () -> decltype(P::get_ptr(i)) {
        return P::get_ptr(i);
    }
};

template<typename T, typename C, typename P>
class iterator_base {
    friend C;
    template<typename N, typename K, typename V>
    friend class const_iterator_base;
    template<typename N, typename K, typename V>
    friend class reverse_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    iterator_base(const T& j) : i(j) {}

public :
    iterator_base() {}

    iterator_base operator ++ (int) {
        iterator_base iter = *this;
        ++i;
        return iter;
    }

    iterator_base& operator ++ () {
        ++i; return *this;
    }

    iterator_base operator -- (int) {
        iterator_base iter = *this;
        --i;
        return iter;
    }

    iterator_base& operator -- () {
        --i; return *this;
    }

    void operator += (int_t n) {
        i += n;
    }

    void operator -= (int_t n) {
        i -= n;
    }

    iterator_base operator + (int_t n) const {
        iterator_base iter = *this;
        iter.i += n;
        return iter;
    }

    iterator_base operator - (int_t n) const {
        iterator_base iter = *this;
        iter.i -= n;
        return iter;
    }

    std::ptrdiff_t operator - (const iterator_base& iter) const {
        return i - iter.i;
    }

    bool operator == (const iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const iterator_base& iter) const {
        return i != iter.i;
    }

    bool operator < (const iterator_base& iter) const {
        return i < iter.i;
    }

    bool operator <= (const iterator_base& iter) const {
        return i <= iter.i;
    }

    bool operator > (const iterator_base& iter) const {
        return i > iter.i;
    }

    bool operator >= (const iterator_base& iter) const {
        return i >= iter.i;
    }

    auto operator * () -> decltype(P::get_obj(i)) {
        return P::get_obj(i);
    }

    auto operator -> () -> decltype(P::get_ptr(i)) {
        return P::get_ptr(i);
    }
};

template<typename T, typename C, typename P>
class reverse_iterator_base;

template<typename T, typename C, typename P>
class const_reverse_iterator_base {
    friend C;

    T i;

    T stdbase() {
        return i;
    }

    const_reverse_iterator_base(const T& j) : i(j) {}

public :
    const_reverse_iterator_base() {}
    const_reverse_iterator_base(const reverse_iterator_base<T,C,P>& iter) : i(iter.i) {}

    const_reverse_iterator_base operator ++ (int) {
        const_reverse_iterator_base iter = *this;
        ++i;
        return iter;
    }

    const_reverse_iterator_base& operator ++ () {
        ++i; return *this;
    }

    const_reverse_iterator_base operator -- (int) {
        const_reverse_iterator_base iter = *this;
        --i;
        return iter;
    }

    const_reverse_iterator_base& operator -- () {
        --i; return *this;
    }

    void operator += (int_t n) {
        i += n;
    }

    void operator -= (int_t n) {
        i -= n;
    }

    const_reverse_iterator_base operator + (int_t n) const {
        const_reverse_iterator_base iter = *this;
        iter.i += n;
        return iter;
    }

    const_reverse_iterator_base operator - (int_t n) const {
        const_reverse_iterator_base iter = *this;
        iter.i -= n;
        return iter;
    }

    std::ptrdiff_t operator - (const const_reverse_iterator_base& iter) const {
        return i - iter.i;
    }

    bool operator == (const const_reverse_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const const_reverse_iterator_base& iter) const {
        return i != iter.i;
    }

    bool operator < (const const_reverse_iterator_base& iter) const {
        return i < iter.i;
    }

    bool operator <= (const const_reverse_iterator_base& iter) const {
        return i <= iter.i;
    }

    bool operator > (const const_reverse_iterator_base& iter) const {
        return i > iter.i;
    }

    bool operator >= (const const_reverse_iterator_base& iter) const {
        return i >= iter.i;
    }

    auto operator * () -> decltype(P::get_obj(i)) {
        return P::get_obj(i);
    }

    auto operator -> () -> decltype(P::get_ptr(i)) {
        return P::get_ptr(i);
    }

    const_iterator_base<T,C,P> base() {
        return const_iterator_base<T,C,P>(i.base());
    }
};

template<typename T, typename C, typename P>
class reverse_iterator_base {
    friend C;
    template<typename N, typename K, typename V>
    friend class const_reverse_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    reverse_iterator_base(const T& j) : i(j) {}

public :
    reverse_iterator_base() {}

    reverse_iterator_base operator ++ (int) {
        reverse_iterator_base iter = *this;
        ++i;
        return iter;
    }

    reverse_iterator_base& operator ++ () {
        ++i; return *this;
    }

    reverse_iterator_base operator -- (int) {
        reverse_iterator_base iter = *this;
        --i;
        return iter;
    }

    reverse_iterator_base& operator -- () {
        --i; return *this;
    }

    void operator += (int_t n) {
        i += n;
    }

    void operator -= (int_t n) {
        i -= n;
    }

    reverse_iterator_base operator + (int_t n) const {
        reverse_iterator_base iter = *this;
        iter.i += n;
        return iter;
    }

    reverse_iterator_base operator - (int_t n) const {
        reverse_iterator_base iter = *this;
        iter.i -= n;
        return iter;
    }

    std::ptrdiff_t operator - (const reverse_iterator_base& iter) const {
        return i - iter.i;
    }

    bool operator == (const reverse_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const reverse_iterator_base& iter) const {
        return i != iter.i;
    }

    bool operator < (const reverse_iterator_base& iter) const {
        return i < iter.i;
    }

    bool operator <= (const reverse_iterator_base& iter) const {
        return i <= iter.i;
    }

    bool operator > (const reverse_iterator_base& iter) const {
        return i > iter.i;
    }

    bool operator >= (const reverse_iterator_base& iter) const {
        return i >= iter.i;
    }

    auto operator * () -> decltype(P::get_obj(i)) {
        return P::get_obj(i);
    }

    auto operator -> () -> decltype(P::get_ptr(i)) {
        return P::get_ptr(i);
    }

    iterator_base<T,C,P> base() {
        return iterator_base<T,C,P>(i.base());
    }
};

#endif
