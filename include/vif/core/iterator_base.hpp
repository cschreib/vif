#ifndef VIF_CORE_ITERATOR_ADAPTOR_HPP
#define VIF_CORE_ITERATOR_ADAPTOR_HPP

namespace vif {
namespace impl {
    template<typename T, typename C, typename P>
    struct iterator_adaptor;

    template<typename T, typename C, typename P>
    struct const_iterator_adaptor {
        P policy;
        T i;

        T stdbase() {
            return i;
        }

        const_iterator_adaptor() = default;
        const_iterator_adaptor(const T& j, P p = P()) : policy(std::move(p)), i(j) {
            policy.initialize(i);
        }
        const_iterator_adaptor(const iterator_adaptor<T,C,P>& iter) : policy(iter.policy), i(iter.i) {}

        const_iterator_adaptor operator ++ (int) {
            const_iterator_adaptor iter = *this;
            policy.increment(i);
            return iter;
        }

        const_iterator_adaptor& operator ++ () {
            policy.increment(i); return *this;
        }

        const_iterator_adaptor operator -- (int) {
            const_iterator_adaptor iter = *this;
            policy.decrement(i);
            return iter;
        }

        const_iterator_adaptor& operator -- () {
            policy.decrement(i); return *this;
        }

        void operator += (int_t n) {
            policy.advance(i, n);
        }

        void operator -= (int_t n) {
            policy.advance(i, -n);
        }

        const_iterator_adaptor operator + (int_t n) const {
            const_iterator_adaptor iter = *this;
            iter += n;
            return iter;
        }

        const_iterator_adaptor operator - (int_t n) const {
            const_iterator_adaptor iter = *this;
            iter -= n;
            return iter;
        }

        std::ptrdiff_t operator - (const const_iterator_adaptor& iter) const {
            return policy.difference(i, iter.i);
        }

        bool operator == (const const_iterator_adaptor& iter) const {
            return i == iter.i;
        }

        bool operator != (const const_iterator_adaptor& iter) const {
            return i != iter.i;
        }

        bool operator < (const const_iterator_adaptor& iter) const {
            return i < iter.i;
        }

        bool operator <= (const const_iterator_adaptor& iter) const {
            return i <= iter.i;
        }

        bool operator > (const const_iterator_adaptor& iter) const {
            return i > iter.i;
        }

        bool operator >= (const const_iterator_adaptor& iter) const {
            return i >= iter.i;
        }

        auto operator * () -> decltype(policy.get_obj(i)) {
            return policy.get_obj(i);
        }

        auto operator -> () -> decltype(policy.get_ptr(i)) {
            return policy.get_ptr(i);
        }
    };

    template<typename T, typename C, typename P>
    struct iterator_adaptor {
        P policy;
        T i;

        T stdbase() {
            return i;
        }

        iterator_adaptor() = default;
        iterator_adaptor(const T& j, P p = P()) : policy(std::move(p)), i(j) {
            policy.initialize(i);
        }

        iterator_adaptor operator ++ (int) {
            iterator_adaptor iter = *this;
            policy.increment(i);
            return iter;
        }

        iterator_adaptor& operator ++ () {
            policy.increment(i); return *this;
        }

        iterator_adaptor operator -- (int) {
            iterator_adaptor iter = *this;
            policy.decrement(i);
            return iter;
        }

        iterator_adaptor& operator -- () {
            policy.decrement(i); return *this;
        }

        void operator += (int_t n) {
            policy.advance(i, n);
        }

        void operator -= (int_t n) {
            policy.advance(i, -n);
        }

        iterator_adaptor operator + (int_t n) const {
            iterator_adaptor iter = *this;
            iter += n;
            return iter;
        }

        iterator_adaptor operator - (int_t n) const {
            iterator_adaptor iter = *this;
            iter -= n;
            return iter;
        }

        std::ptrdiff_t operator - (const iterator_adaptor& iter) const {
            return policy.difference(i, iter.i);
        }

        bool operator == (const iterator_adaptor& iter) const {
            return i == iter.i;
        }

        bool operator != (const iterator_adaptor& iter) const {
            return i != iter.i;
        }

        bool operator < (const iterator_adaptor& iter) const {
            return i < iter.i;
        }

        bool operator <= (const iterator_adaptor& iter) const {
            return i <= iter.i;
        }

        bool operator > (const iterator_adaptor& iter) const {
            return i > iter.i;
        }

        bool operator >= (const iterator_adaptor& iter) const {
            return i >= iter.i;
        }

        auto operator * () -> decltype(policy.get_obj(i)) {
            return policy.get_obj(i);
        }

        auto operator -> () -> decltype(policy.get_ptr(i)) {
            return policy.get_ptr(i);
        }
    };

    template<typename T, typename C, typename P>
    struct reverse_iterator_adaptor;

    template<typename T, typename C, typename P>
    struct const_reverse_iterator_adaptor {
        P policy;
        T i;

        T stdbase() {
            return i;
        }

        const_reverse_iterator_adaptor() = default;
        const_reverse_iterator_adaptor(const T& j, P p = P()) : policy(std::move(p)), i(j) {
            policy.initialize(i);
        }
        const_reverse_iterator_adaptor(const reverse_iterator_adaptor<T,C,P>& iter) : policy(iter.policy), i(iter.i) {}

        const_reverse_iterator_adaptor operator ++ (int) {
            const_reverse_iterator_adaptor iter = *this;
            policy.increment(i);
            return iter;
        }

        const_reverse_iterator_adaptor& operator ++ () {
            policy.increment(i); return *this;
        }

        const_reverse_iterator_adaptor operator -- (int) {
            const_reverse_iterator_adaptor iter = *this;
            policy.decrement(i);
            return iter;
        }

        const_reverse_iterator_adaptor& operator -- () {
            policy.decrement(i); return *this;
        }

        void operator += (int_t n) {
            policy.advance(i, n);
        }

        void operator -= (int_t n) {
            policy.advance(i, -n);
        }

        const_reverse_iterator_adaptor operator + (int_t n) const {
            const_reverse_iterator_adaptor iter = *this;
            iter += n;
            return iter;
        }

        const_reverse_iterator_adaptor operator - (int_t n) const {
            const_reverse_iterator_adaptor iter = *this;
            iter -= n;
            return iter;
        }

        std::ptrdiff_t operator - (const const_reverse_iterator_adaptor& iter) const {
            return policy.difference(i, iter.i);
        }

        bool operator == (const const_reverse_iterator_adaptor& iter) const {
            return i == iter.i;
        }

        bool operator != (const const_reverse_iterator_adaptor& iter) const {
            return i != iter.i;
        }

        bool operator < (const const_reverse_iterator_adaptor& iter) const {
            return i < iter.i;
        }

        bool operator <= (const const_reverse_iterator_adaptor& iter) const {
            return i <= iter.i;
        }

        bool operator > (const const_reverse_iterator_adaptor& iter) const {
            return i > iter.i;
        }

        bool operator >= (const const_reverse_iterator_adaptor& iter) const {
            return i >= iter.i;
        }

        auto operator * () -> decltype(policy.get_obj(i)) {
            return policy.get_obj(i);
        }

        auto operator -> () -> decltype(policy.get_ptr(i)) {
            return policy.get_ptr(i);
        }

        const_iterator_adaptor<T,C,P> base() {
            return const_iterator_adaptor<T,C,P>(i.base());
        }
    };

    template<typename T, typename C, typename P>
    struct reverse_iterator_adaptor {
        P policy;
        T i;

        T stdbase() {
            return i;
        }

        reverse_iterator_adaptor() = default;
        reverse_iterator_adaptor(const T& j, P p = P()) : policy(std::move(p)), i(j) {
            policy.initialize(i);
        }

        reverse_iterator_adaptor operator ++ (int) {
            reverse_iterator_adaptor iter = *this;
            policy.increment(i);
            return iter;
        }

        reverse_iterator_adaptor& operator ++ () {
            policy.increment(i); return *this;
        }

        reverse_iterator_adaptor operator -- (int) {
            reverse_iterator_adaptor iter = *this;
            policy.decrement(i);
            return iter;
        }

        reverse_iterator_adaptor& operator -- () {
            policy.decrement(i); return *this;
        }

        void operator += (int_t n) {
            policy.advance(i, n);
        }

        void operator -= (int_t n) {
            policy.advance(i, -n);
        }

        reverse_iterator_adaptor operator + (int_t n) const {
            reverse_iterator_adaptor iter = *this;
            iter += n;
            return iter;
        }

        reverse_iterator_adaptor operator - (int_t n) const {
            reverse_iterator_adaptor iter = *this;
            iter -= n;
            return iter;
        }

        std::ptrdiff_t operator - (const reverse_iterator_adaptor& iter) const {
            return policy.difference(i, iter.i);
        }

        bool operator == (const reverse_iterator_adaptor& iter) const {
            return i == iter.i;
        }

        bool operator != (const reverse_iterator_adaptor& iter) const {
            return i != iter.i;
        }

        bool operator < (const reverse_iterator_adaptor& iter) const {
            return i < iter.i;
        }

        bool operator <= (const reverse_iterator_adaptor& iter) const {
            return i <= iter.i;
        }

        bool operator > (const reverse_iterator_adaptor& iter) const {
            return i > iter.i;
        }

        bool operator >= (const reverse_iterator_adaptor& iter) const {
            return i >= iter.i;
        }

        auto operator * () -> decltype(policy.get_obj(i)) {
            return policy.get_obj(i);
        }

        auto operator -> () -> decltype(policy.get_ptr(i)) {
            return policy.get_ptr(i);
        }

        iterator_adaptor<T,C,P> base() {
            return iterator_adaptor<T,C,P>(i.base());
        }
    };
}
}

#endif
