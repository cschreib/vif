#ifndef ITERATOR_HPP
#define ITERATOR_HPP

template<typename T, typename C>
class ptr_iterator_base;

template<typename T, typename C>
class const_ptr_iterator_base {
    friend C;
    template<typename N, typename K>
    friend class const_reverse_ptr_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    const_ptr_iterator_base(const T& j) : i(j) {}

public :
    const_ptr_iterator_base() {}
    const_ptr_iterator_base(const ptr_iterator_base<T,C>& iter) : i(iter.i) {}

    const_ptr_iterator_base& operator ++ (int) {
        ++i; return *this;
    }

    const_ptr_iterator_base operator ++ () {
        ++i; return *this;
    }

    bool operator == (const const_ptr_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const const_ptr_iterator_base& iter) const {
        return i != iter.i;
    }

    auto operator - (const const_ptr_iterator_base& iter) const -> decltype(i - iter.i) {
        return i - iter.i;
    }

    auto operator * () -> decltype(**i) {
        return **i;
    }

    auto operator -> () -> decltype(*i) {
        return *i;
    }
};

template<typename T, typename C>
class ptr_iterator_base {
    friend C;
    template<typename N, typename K>
    friend class const_ptr_iterator_base;
    template<typename N, typename K>
    friend class reverse_ptr_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    ptr_iterator_base(const T& j) : i(j) {}

public :
    ptr_iterator_base() {}

    ptr_iterator_base& operator ++ (int) {
        ++i; return *this;
    }

    ptr_iterator_base operator ++ () {
        ++i; return *this;
    }

    bool operator == (const ptr_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const ptr_iterator_base& iter) const {
        return i != iter.i;
    }

    auto operator - (const ptr_iterator_base& iter) const -> decltype(i - iter.i) {
        return i - iter.i;
    }

    auto operator * () -> decltype(**i) {
        return **i;
    }

    auto operator -> () -> decltype(*i) {
        return *i;
    }
};

template<typename T, typename C>
class reverse_ptr_iterator_base;

template<typename T, typename C>
class const_reverse_ptr_iterator_base {
    friend C;

    T i;

    T stdbase() {
        return i;
    }

    const_reverse_ptr_iterator_base(const T& j) : i(j) {}

public :
    const_reverse_ptr_iterator_base() {}
    const_reverse_ptr_iterator_base(const reverse_ptr_iterator_base<T,C>& iter) : i(iter.i) {}

    const_reverse_ptr_iterator_base& operator ++ (int) {
        ++i; return *this;
    }

    const_reverse_ptr_iterator_base operator ++ () {
        ++i; return *this;
    }

    bool operator == (const const_reverse_ptr_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const const_reverse_ptr_iterator_base& iter) const {
        return i != iter.i;
    }

    auto operator - (const const_reverse_ptr_iterator_base& iter) const -> decltype(i - iter.i) {
        return i - iter.i;
    }

    auto operator * () -> decltype(**i) {
        return **i;
    }

    auto operator -> () -> decltype(*i) {
        return *i;
    }

    const_ptr_iterator_base<T,C> base() {
        return const_ptr_iterator_base<T,C>(i.base());
    }
};

template<typename T, typename C>
class reverse_ptr_iterator_base {
    friend C;
    template<typename N, typename K>
    friend class const_reverse_ptr_iterator_base;

    T i;

    T stdbase() {
        return i;
    }

    reverse_ptr_iterator_base(const T& j) : i(j) {}

public :
    reverse_ptr_iterator_base() {}

    reverse_ptr_iterator_base& operator ++ (int) {
        ++i; return *this;
    }

    reverse_ptr_iterator_base operator ++ () {
        ++i; return *this;
    }

    bool operator == (const reverse_ptr_iterator_base& iter) const {
        return i == iter.i;
    }

    bool operator != (const reverse_ptr_iterator_base& iter) const {
        return i != iter.i;
    }

    auto operator - (const reverse_ptr_iterator_base& iter) const -> decltype(i - iter.i) {
        return i - iter.i;
    }

    auto operator * () -> decltype(**i) {
        return **i;
    }

    auto operator -> () -> decltype(*i) {
        return *i;
    }

    ptr_iterator_base<T,C> base() {
        return ptr_iterator_base<T,C>(i.base());
    }
};

#endif
