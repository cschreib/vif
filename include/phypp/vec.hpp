#ifndef PHYPP_VEC_HPP
#define PHYPP_VEC_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <initializer_list>
#include "phypp/range.hpp"
#include "phypp/variadic.hpp"
#include "phypp/error.hpp"

// Tag type to mark initialization of a reference vector.
static struct vec_ref_tag_t {} vec_ref_tag;
// Tag type to permit copy initialization without data copy.
static struct vec_nocopy_tag_t {} vec_nocopy_tag;

// Tag a value to indicate that is has to repeated N times
template<std::size_t N, typename T>
struct repeated_value {
    T value;
};

template<std::size_t N, typename T>
repeated_value<N,T> repeat(T t) {
    return repeated_value<N,T>{t};
}

// Load a bunch of helper code to make this file more readable
#include "phypp/vec_helpers.hpp"
#include "phypp/vec_iterator.hpp"
#include "phypp/vec_access.hpp"
#include "phypp/vec_initializer_list.hpp"

////////////////////////////////////////////
//            Generic vector              //
////////////////////////////////////////////

template<std::size_t Dim, typename Type>
struct vec;

// Disable dimension zero
template<typename Type>
struct vec<0,Type> {};

// The vector
template<std::size_t Dim, typename Type>
struct vec {
    static_assert(!std::is_pointer<Type>::value, "library bug: pointer specialization failed");

    using effective_type = vec;
    using rtype = rtype_t<Type>;
    using dtype = dtype_t<Type>;
    using drtype = dtype_t<Type>;
    using vtype = std::vector<dtype>;
    using dim_type = std::array<std::size_t, Dim>;
    struct comparator {
        constexpr bool operator() (const dtype& t1, const dtype& t2) const {
            return t1 < t2;
        }
        template<typename U>
        constexpr bool operator() (const dtype& t1, const U& t2) const {
            return t1 < t2;
        }
        template<typename U>
        constexpr bool operator() (const U& t1, const dtype& t2) const {
            return t1 < t2;
        }
    };

    vtype    data;
    dim_type dims = {{0}};

    // Default constructor
    vec() : safe(*this) {}

    // Copy constructor
    vec(const vec& v) : data(v.data), dims(v.dims), safe(*this) {}
    // Copy constructor (only dims)
    vec(vec_nocopy_tag_t, const vec& v) : dims(v.dims), safe(*this) {}

    // Move constructor
    vec(vec&& v) : data(std::move(v.data)), dims(v.dims), safe(*this) {
        for (uint_t i : range(Dim)) {
            v.dims[i] = 0;
        }
    }

    // Dimension constructor
    template<typename ... Args, typename enable =
        typename std::enable_if<is_dim_list<Args...>::value>::type>
    explicit vec(Args&& ... d) : safe(*this) {
        static_assert(dim_total<Args...>::value == Dim, "dimension list does not match "
            "the dimensions of this vector");

        set_array(dims, std::forward<Args>(d)...);
        resize();
    }

    // Initializer list constructor
    vec(nested_initializer_list<Dim,dtype_t<Type>> il) : safe(*this) {
        vec_ilist::helper<Dim, Type>::fill(*this, il);
    }

    // Implicit conversion
    template<typename T, typename std::enable_if<!vec_only_explicit_convertible<rtype_t<T>,Type>::value, bool>::type = false>
    vec(const vec<Dim,T>& v) : dims(v.dims), safe(*this) {
        static_assert(vec_implicit_convertible<rtype_t<T>,Type>::value,
            "could not construct vector from non-convertible type");
        data.resize(v.data.size());
        for (uint_t i : range(v)) {
            data[i] = v.safe[i];
        }
    }

    // Explicit conversion
    template<typename T, typename std::enable_if<vec_only_explicit_convertible<rtype_t<T>,Type>::value, bool>::type = false>
    explicit vec(const vec<Dim,T>& v) : dims(v.dims), safe(*this) {
        static_assert(vec_explicit_convertible<rtype_t<T>,Type>::value,
            "could not construct vector from non-convertible type");
        data.resize(v.data.size());
        for (uint_t i : range(v)) {
            data[i] = v.safe[i];
        }
    }

    vec& operator = (nested_initializer_list<Dim,dtype_t<Type>> il) {
        vec_ilist::helper<Dim, Type>::fill(*this, il);
        return *this;
    }

    vec& operator = (const vec& v) {
        data = v.data;
        dims = v.dims;
        return *this;
    }

    vec& operator = (vec&& v) {
        data = std::move(v.data);
        dims = std::move(v.dims);
        return *this;
    }

    template<typename T>
    vec& operator = (const vec<Dim,T>& v) {
        static_assert(vec_implicit_convertible<T,Type>::value, "could not assign vectors of "
            "non-implicitly-convertible types");

        if (v.view_same(*this)) {
            std::vector<dtype> t(v.data.size());
            for (uint_t i : range(v)) {
                t[i] = v.safe[i];
            }
            dims = v.dims;
            data.resize(v.data.size());
            for (uint_t i : range(v)) {
                data[i] = t[i];
            }
        } else {
            dims = v.dims;
            data.resize(v.data.size());
            for (uint_t i : range(v)) {
                data[i] = v.safe[i];
            }
        }

        return *this;
    }

    template<std::size_t D, typename T>
    bool view_same(const vec<D,T>&) const {
        return false;
    }

    template<std::size_t D>
    bool view_same(const vec<D,Type*>& v) const {
        return v.view_same(*this);
    }

    bool empty() const {
        return data.empty();
    }

    std::size_t size() const {
        return data.size();
    }

    void resize() {
        std::size_t size = 1;
        for (uint_t i : range(Dim)) {
            size *= dims[i];
        }

        data.resize(size);
    }

    template<typename ... Args>
    void resize(Args&& ... d) {
        set_array(dims, d...);
        resize();
    }

    void clear() {
        data.clear();
        for (uint_t i : range(Dim)) {
            dims[i] = 0;
        }
    }

    const Type& back() const {
        static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
        return reinterpret_cast<const Type&>(data.back());
    }

    Type& back() {
        static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
        return const_cast<Type&>(const_cast<const vec&>(*this).back());
    }

    const Type& front() const {
        static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
        return reinterpret_cast<const Type&>(data.front());
    }

    Type& front() {
        static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
        return const_cast<Type&>(const_cast<const vec&>(*this).front());
    }

    void push_back(const Type& t) {
        static_assert(Dim == 1, "cannot call push_back(Type) on multidimensional vectors");
        data.push_back(t);
        ++dims[0];
    }

    void push_back(Type&& t) {
        static_assert(Dim == 1, "cannot call push_back(Type) on multidimensional vectors");
        data.push_back(std::move(t));
        ++dims[0];
    }

    template<typename T = Type, typename enable =
        typename std::enable_if<std::is_convertible<T,Type>::value>::type>
    void push_back(const vec<Dim-1,T>& t) {
        static_assert(Dim > 1, "cannot call push_back(vec<D-1>) on monodimensional vectors");
        if (empty()) {
            dims[0] = 1;
            for (uint_t i = 1; i < Dim; ++i) {
                dims[i] = t.dims[i-1];
            }

            data.resize(t.data.size());
            std::copy(t.data.begin(), t.data.end(), data.begin());
        } else {
            for (uint_t i = 0; i < Dim-1; ++i) {
                phypp_check(dims[i+1] == t.dims[i],
                    "push_back: incompatible dimensions (", dims, " vs ", t.dims, ")");
            }

            data.insert(data.end(), t.data.begin(), t.data.end());
            ++dims[0];
        }
    }

    template<typename T = Type, typename enable =
        typename std::enable_if<std::is_convertible<T,Type>::value>::type>
    void push_back(const vec<Dim-1,T*>& t) {
        static_assert(Dim > 1, "cannot call push_back(vec<D-1>) on monodimensional vectors");
        push_back(t.concretise());
    }

    void reserve(uint_t n) {
        data.reserve(n);
    }

    const vec<Dim,Type>& concretise() const {
        return *this;
    }

    template<typename T>
    uint_t to_idx_(T ui, cte_t<false>) const {
        phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
            data.size(), ")");
        return ui;
    }

    template<typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += data.size();
        phypp_check(i >= 0, "operator[]: index out of bounds (", i+data.size(), " vs. ",
            data.size(), ")");
        uint_t ui(i);
        phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
            data.size(), ")");
        return ui;
    }

    template<std::size_t D, typename T>
    uint_t to_idx_(T ui, cte_t<false>) const {
        phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
            dims[D], ")");
        return ui;
    }

    template<std::size_t D, typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += dims[D];
        phypp_check(i >= 0, "operator(): index out of bounds (", i+data.size(), " vs. ",
            dims[D], ")");
        uint_t ui(i);
        phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
            dims[D], ")");
        return ui;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    uint_t to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    uint_t to_idx(T ix) const {
        return to_idx_<D>(ix, cte_t<std::is_signed<T>::value>());
    }

    uint_t pitch(uint_t i) const {
        uint_t p = 1;
        for (uint_t j = i+1; j < Dim; ++j) {
            p *= dims[j];
        }
        return p;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    Type& operator [] (T i) {
        return const_cast<Type&>(const_cast<const vec&>(*this)[i]);
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    const Type& operator [] (T i) const {
        return reinterpret_cast<const Type&>(data[to_idx(i)]);
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,Type*> operator [] (const vec<1,T>& i) {
        vec<1,Type*> v(vec_ref_tag, *this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = ptr<Type>(data[to_idx(i.safe[j])]);
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,Type*> operator [] (const vec<1,T*>& i) {
        vec<1,Type*> v(vec_ref_tag, *this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = ptr<Type>(data[to_idx(i.safe[j])]);
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,const Type*> operator [] (const vec<1,T>& i) const {
        vec<1,const Type*> v(vec_ref_tag, *this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = ptr<Type>(data[to_idx(i.safe[j])]);
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
        vec<1,const Type*> v(vec_ref_tag, *this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = ptr<Type>(data[to_idx(i.safe[j])]);
        }
        return v;
    }

    vec<1,Type*> operator [] (full_range_t rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const left_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const right_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const left_right_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (full_range_t rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const left_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const right_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const left_right_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) ->
        typename vec_access::helper<true, false, Dim, Type, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, false, Dim, Type, Args...>::access(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const ->
        typename vec_access::helper<true, true, Dim, Type, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, true, Dim, Type, Args...>::access(*this, i...);
    }

    struct safe_proxy {
        vec& parent;

        safe_proxy(vec& p) : parent(p) {}
        safe_proxy(const vec&) = delete;
        safe_proxy(vec&&) = delete;
        safe_proxy& operator=(const vec&) = delete;
        safe_proxy& operator=(vec&&) = delete;

        Type& operator [] (uint_t i) {
            return const_cast<Type&>(const_cast<const safe_proxy&>(*this)[i]);
        }

        const Type& operator [] (uint_t i) const {
            return reinterpret_cast<const Type&>(parent.data[i]);
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T>& i) {
            vec<1,Type*> v(vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = ptr<Type>(parent.data[i.safe[j]]);
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T*>& i) {
            vec<1,Type*> v(vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = ptr<Type>(parent.data[i.safe[j]]);
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T>& i) const {
            vec<1,const Type*> v(vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = ptr<Type>(parent.data[i.safe[j]]);
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
            vec<1,const Type*> v(vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = ptr<Type>(parent.data[i.safe[j]]);
            }
            return v;
        }

        vec<1,Type*> operator [] (full_range_t rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const left_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const right_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const left_right_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (full_range_t rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const left_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const right_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const left_right_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) ->
            typename vec_access::helper<false, false, Dim, Type, Args...>::type {
            static_assert(vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return vec_access::helper<false, false, Dim, Type, Args...>::access(parent, i...);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) const ->
            typename vec_access::helper<false, true, Dim, Type, Args...>::type {
            static_assert(vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return vec_access::helper<false, true, Dim, Type, Args...>::access(parent, i...);
        }
    } safe;

    vec operator + () const {
        return *this;
    }

    vec operator - () const {
        vec v = *this;
        for (auto& t : v) {
            t = -t;
        }

        return v;
    }

    #define OPERATOR(op) \
        template<typename U> \
        vec& operator op (const vec<Dim,U>& u) { \
            phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                "' ("+strn(dims)+" vs "+strn(u.dims)+")"); \
            if (u.view_same(*this)) { \
                std::vector<dtype> t; t.resize(data.size()); \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    t[i] = u.safe[i]; \
                } \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    data[i] op t[i]; \
                } \
            } else { \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    data[i] op u.safe[i]; \
                } \
            } \
            return *this; \
        } \
        template<typename U> \
        vec& operator op (U u) { \
            for (auto& v : data) { \
                v op u; \
            } \
            return *this; \
        }

    OPERATOR(*=)
    OPERATOR(/=)
    OPERATOR(+=)
    OPERATOR(-=)

    #undef OPERATOR

    using iterator = typename vec_iterator_type<vec>::iterator;
    using const_iterator = typename vec_iterator_type<vec>::const_iterator;

    iterator begin() {
        return data.begin();
    }

    iterator end() {
        return data.end();
    }

    const_iterator begin() const {
        return data.begin();
    }

    const_iterator end() const {
        return data.end();
    }
};


////////////////////////////////////////////
//              Vector view               //
////////////////////////////////////////////

template<std::size_t Dim, typename Type>
struct vec<Dim,Type*> {
    using rtype = rtype_t<Type*>;
    using effective_type = vec<Dim,rtype>;
    using dtype = Type;
    using drtype = rtype;
    using vtype = std::vector<dtype*>;
    using dim_type = std::array<std::size_t, Dim>;
    struct comparator {
        bool operator() (const dtype* t1, const dtype* t2) {
            return *t1 < *t2;
        }
        template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
        bool operator() (const dtype* t1, const U& t2) {
            return *t1 < t2;
        }
        template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
        bool operator() (const U& t1, const dtype* t2) {
            return t1 < *t2;
        }

        constexpr bool operator() (const dtype& t1, const dtype& t2) const {
            return t1 < t2;
        }
        template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
        constexpr bool operator() (const dtype& t1, const U& t2) const {
            return t1 < t2;
        }
        template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
        constexpr bool operator() (const U& t1, const dtype& t2) const {
            return t1 < t2;
        }
    };

    void*    parent = nullptr;
    vtype    data;
    dim_type dims = {{0}};

    vec() = delete;
    vec(const vec& v) : parent(v.parent), data(v.data), dims(v.dims), safe(*this) {}
    vec(vec_nocopy_tag_t, const vec& v) : parent(v.parent), dims(v.dims), safe(*this) {}
    vec(vec&& v) : parent(v.parent), data(std::move(v.data)), dims(std::move(v.dims)), safe(*this) {}

    template<std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_same<rtype_t<T>, rtype>::value>::type>
    explicit vec(vec_ref_tag_t, const vec<D,T>& p) : safe(*this) {
        parent = static_cast<void*>(const_cast<vec<D,T>*>(&p));
    }

    explicit vec(vec_ref_tag_t, void* p) : safe(*this) {
        parent = p;
    }

    vec& operator = (const Type& t) {
        for (uint_t i = 0; i < data.size(); ++i) {
            *data[i] = t;
        }
        return *this;
    }

    vec& operator = (nested_initializer_list<Dim,dtype_t<Type>> il) {
        vec_ilist::helper<Dim, Type*>::fill(*this, il);
        return *this;
    }

    template<typename T>
    void assign_(const vec<Dim,T>& v) {
        if (view_same(v)) {
            // Make a copy to prevent aliasing
            std::vector<dtype> t(v.data.size());
            for (uint_t i : range(v)) {
                t[i] = v.safe[i];
            }

            // Actual assignment
            for (uint_t i : range(v)) {
                *data[i] = t[i];
            }
        } else {
            // No aliasing possible, assign directly
            for (uint_t i : range(v)) {
                *data[i] = v.safe[i];
            }
        }
    }

    template<typename T>
    vec& operator = (const vec<Dim,T>& v) {
        static_assert(vec_implicit_convertible<T,Type>::value, "could not assign vectors of "
            "non-implicitly-convertible types");
        phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
            v.dims, " to ", dims, ")");

        assign_(v);

        return *this;
    }

    vec& operator = (const vec& v) {
        phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
            v.dims, " to ", dims, ")");

        assign_(v);

        return *this;
    }

    template<std::size_t D, typename T>
    bool view_same(const vec<D,T>& v) const {
        return false;
    }

    template<std::size_t D>
    bool view_same(const vec<D,Type*>& v) const {
        return v.parent == parent;
    }

    template<std::size_t D>
    bool view_same(const vec<D,Type>& v) const {
        return static_cast<void*>(const_cast<vec<D,Type>*>(&v)) == parent;
    }

    bool empty() const {
        return data.empty();
    }

    std::size_t size() const {
        return data.size();
    }

    void resize() {
        std::size_t size = 1;
        for (uint_t i = 0; i < Dim; ++i) {
            size *= dims[i];
        }

        data.resize(size);
    }

    effective_type concretise() const {
        return *this;
    }

    const Type& back() const {
        static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
        return reinterpret_cast<const Type&>(*data.back());
    }

    Type& back() {
        static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
        return const_cast<Type&>(const_cast<const vec&>(*this).back());
    }

    const Type& front() const {
        static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
        return reinterpret_cast<const Type&>(*data.front());
    }

    Type& front() {
        static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
        return const_cast<Type&>(const_cast<const vec&>(*this).front());
    }

    template<typename T>
    T to_idx_(T ui, cte_t<false>) const {
        phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
            data.size(), ")");
        return ui;
    }

    template<typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += data.size();
        phypp_check(i >= 0, "operator[]: index out of bounds (", i+data.size(), " vs. ",
            data.size(), ")");
        uint_t ui(i);
        phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
            data.size(), ")");
        return ui;
    }

    template<std::size_t D, typename T>
    T to_idx_(T ui, cte_t<false>) const {
        phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
            dims[D], ")");
        return ui;
    }

    template<std::size_t D, typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += dims[D];
        phypp_check(i >= 0, "operator(): index out of bounds (", i+data.size(), " vs. ",
            dims[D], ")");
        uint_t ui(i);
        phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
            dims[D], ")");
        return ui;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_<D>(ix, cte_t<std::is_signed<T>::value>());
    }

    uint_t pitch(uint_t i) const {
        uint_t p = 1;
        for (uint_t j = i+1; j < Dim; ++j) {
            p *= dims[j];
        }
        return p;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    Type& operator [] (T i) {
        return const_cast<Type&>(const_cast<const vec&>(*this)[i]);
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    const Type& operator [] (T i) const {
        return *data[to_idx(i)];
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,Type*> operator [] (const vec<1,T>& i) {
        vec<1,Type*> v(vec_ref_tag, parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i.safe[j])];
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,Type*> operator [] (const vec<1,T*>& i) {
        vec<1,Type*> v(vec_ref_tag, parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i.safe[j])];
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,const Type*> operator [] (const vec<1,T>& i) const {
        vec<1,const Type*> v(vec_ref_tag, parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i.safe[j])];
        }
        return v;
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_integral<T>::value>::type>
    vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
        vec<1,const Type*> v(vec_ref_tag, parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i.safe[j])];
        }
        return v;
    }

    vec<1,Type*> operator [] (full_range_t rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const left_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const right_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,Type*> operator [] (const left_right_range_t& rng) {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (full_range_t rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const left_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const right_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    vec<1,const Type*> operator [] (const left_right_range_t& rng) const {
        return vec_access::bracket_access(*this, rng);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) ->
        typename vec_access::helper<true, false, Dim, Type*, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, false, Dim, Type*, Args...>::access(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const ->
        typename vec_access::helper<true, true, Dim, Type*, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, true, Dim, Type*, Args...>::access(*this, i...);
    }

    struct safe_proxy {
        vec& parent;

        safe_proxy(vec& p) : parent(p) {}
        safe_proxy(const vec&) = delete;
        safe_proxy(vec&&) = delete;
        safe_proxy& operator=(const vec&) = delete;
        safe_proxy& operator=(vec&&) = delete;

        Type& operator [] (uint_t i) {
            return const_cast<Type&>(const_cast<const safe_proxy&>(*this)[i]);
        }

        const Type& operator [] (uint_t i) const {
            return *parent.data[i];
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T>& i) {
            vec<1,Type*> v(vec_ref_tag, parent.parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = parent.data[i.safe[j]];
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T*>& i) {
            vec<1,Type*> v(vec_ref_tag, parent.parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = parent.data[i.safe[j]];
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T>& i) const {
            vec<1,const Type*> v(vec_ref_tag, parent.parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = parent.data[i.safe[j]];
            }
            return v;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_unsigned<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
            vec<1,const Type*> v(vec_ref_tag, parent.parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j = 0; j < i.data.size(); ++j) {
                v.data[j] = parent.data[i.safe[j]];
            }
            return v;
        }

        vec<1,Type*> operator [] (full_range_t rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const left_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const right_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,Type*> operator [] (const left_right_range_t& rng) {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (full_range_t rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const left_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const right_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        vec<1,const Type*> operator [] (const left_right_range_t& rng) const {
            return vec_access::bracket_access(parent, rng);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) ->
            typename vec_access::helper<false, false, Dim, Type*, Args...>::type {
            static_assert(vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return vec_access::helper<false, false, Dim, Type*, Args...>::access(parent, i...);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) const ->
            typename vec_access::helper<false, true, Dim, Type*, Args...>::type {
            static_assert(vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return vec_access::helper<false, true, Dim, Type*, Args...>::access(parent, i...);
        }
    } safe;

    effective_type operator + () const {
        return *this;
    }

    effective_type operator - () const {
        effective_type v = *this;
        for (auto& t : v) {
            t = -t;
        }

        return v;
    }

    #define OPERATOR(op) \
        template<typename U> \
        vec& operator op (const vec<Dim,U>& u) { \
            phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                "' ("+strn(dims)+" vs "+strn(u.dims)+")"); \
            if (u.view_same(*this)) { \
                std::vector<dtype> t; t.resize(data.size()); \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    t[i] = u.safe[i]; \
                } \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    *data[i] op t[i]; \
                } \
            } else { \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    *data[i] op u.safe[i]; \
                } \
            } \
            return *this; \
        } \
        \
        template<typename U> \
        vec& operator op (U u) { \
            for (auto& v : data) { \
                *v op u; \
            } \
            return *this; \
        }

    OPERATOR(*=)
    OPERATOR(/=)
    OPERATOR(%=)
    OPERATOR(+=)
    OPERATOR(-=)

    #undef OPERATOR

    using base_iterator = typename vtype::iterator;
    using base_const_iterator = typename vtype::const_iterator;

    using iterator = ptr_iterator_base<base_iterator, vec>;
    using const_iterator = const_ptr_iterator_base<base_const_iterator, vec>;

    iterator begin() {
        return data.begin();
    }

    iterator end() {
        return data.end();
    }

    const_iterator begin() const {
        return data.begin();
    }

    const_iterator end() const {
        return data.end();
    }

};


////////////////////////////////////////////
//            Type shortcuts              //
////////////////////////////////////////////

// Shortcuts for the most used types: vec1f, vec2d, ... up to vec6.
#define MAKE_TYPEDEFS(N) \
    using vec##N##f = vec<N, float>; \
    using vec##N##d = vec<N, double>; \
    using vec##N##i = vec<N, int_t>; \
    using vec##N##u = vec<N, std::size_t>; \
    using vec##N##s = vec<N, std::string>; \
    using vec##N##c = vec<N, char>; \
    using vec##N##b = vec<N, bool>;

MAKE_TYPEDEFS(1)
MAKE_TYPEDEFS(2)
MAKE_TYPEDEFS(3)
MAKE_TYPEDEFS(4)
MAKE_TYPEDEFS(5)
MAKE_TYPEDEFS(6)

#undef MAKE_TYPEDEFS


////////////////////////////////////////////
//               Operators                //
////////////////////////////////////////////

// Mathematical operators
struct op_mul_t;
struct op_div_t;
struct op_mod_t;
struct op_add_t;
struct op_sub_t;

struct op_node_t {
    op_mul_t operator * (op_node_t);
    op_div_t operator / (op_node_t);
    op_mod_t operator % (op_node_t);
    op_add_t operator + (op_node_t);
    op_sub_t operator - (op_node_t);
};

#define OP_TYPE(op) decltype(op_node_t{} op op_node_t{})

template<typename T>
using math_bake_type = typename std::decay<typename std::remove_pointer<
    data_type_t<typename std::remove_pointer<T>::type>>::type>::type;

template<typename OP, typename T, typename U>
struct op_res_t;

template<typename T, typename U>
struct op_res_t<op_mul_t, T, U> {
    using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() *
        std::declval<math_bake_type<U>>())>::type;
};

template<typename T, typename U>
struct op_res_t<op_div_t, T, U> {
    using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() /
        std::declval<math_bake_type<U>>())>::type;
};

template<typename T, typename U>
struct op_res_t<op_mod_t, T, U> {
    using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() %
        std::declval<math_bake_type<U>>())>::type;
};

template<typename T, typename U>
struct op_res_t<op_add_t, T, U> {
    using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() +
        std::declval<math_bake_type<U>>())>::type;
};

template<typename T, typename U>
struct op_res_t<op_sub_t, T, U> {
    using type = typename std::decay<decltype(std::declval<math_bake_type<T>>() -
        std::declval<math_bake_type<U>>())>::type;
};

template<typename T, typename U>
using res_mul_t = typename op_res_t<op_mul_t, T, U>::type;
template<typename T, typename U>
using res_div_t = typename op_res_t<op_div_t, T, U>::type;
template<typename T, typename U>
using res_mod_t = typename op_res_t<op_mod_t, T, U>::type;
template<typename T, typename U>
using res_add_t = typename op_res_t<op_add_t, T, U>::type;
template<typename T, typename U>
using res_sub_t = typename op_res_t<op_sub_t, T, U>::type;

template<typename T>
const T& get_element_(const T& t, uint_t i) {
    return t;
}

template<std::size_t Dim, typename T>
auto get_element_(const vec<Dim,T>& t, uint_t i) -> decltype(t.safe[i]) {
    return t.safe[i];
}

#define VECTORIZE(op, sop) \
    template<std::size_t Dim, typename T, typename U> \
    vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec<Dim,T>& v, const vec<Dim,U>& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i : range(v)) { \
            tv.data.push_back(get_element_(v, i) op get_element_(u, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U> \
    vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec<Dim,T>& v, const U& u) { \
        vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i : range(v)) { \
            tv.data.push_back(get_element_(v, i) op u); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
    vec<Dim,T> operator op (vec<Dim,T>&& v, const vec<Dim,U>& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        for (uint_t i : range(v)) { \
            v.data[i] sop get_element_(u, i); \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
    vec<Dim,T> operator op (vec<Dim,T>&& v, const U& u) { \
        for (auto& t : v) { \
            t sop u; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        !is_vec<U>::value>::type> \
    vec<Dim,typename op_res_t<OP_TYPE(op),U,T>::type> operator op (const U& u, const vec<Dim,T>& v) { \
        vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i : range(v)) { \
            tv.data.push_back(u op get_element_(v, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
    vec<Dim,T> operator op (const vec<Dim,U>& u, vec<Dim,T>&& v) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        for (uint_t i : range(v)) { \
            v.data[i] = get_element_(u, i) op v.data[i]; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
    vec<Dim,T> operator op (const U& u, vec<Dim,T>&& v) { \
        for (auto& t : v) { \
            t = u op t; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value && \
        !std::is_pointer<U>::value>::type> \
    vec<Dim,T> operator op (vec<Dim,T>&& v, vec<Dim,U>&& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        for (uint_t i : range(v)) { \
            v.data[i] sop u.data[i]; \
        } \
        return std::move(v); \
    }

VECTORIZE(*, *=)
VECTORIZE(+, +=)
VECTORIZE(/, /=)
VECTORIZE(%, %=)
VECTORIZE(-, -=)

#undef VECTORIZE

// Logical operators
#define VECTORIZE(op) \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!is_vec<U>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v, const U& u) { \
        vec<Dim,bool> tv(v.dims); \
        for (uint_t i : range(v)) { \
            tv.safe[i] = (v.safe[i] op u); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!is_vec<U>::value>::type> \
    vec<Dim,bool> operator op (const U& u, const vec<Dim,T>& v) { \
        vec<Dim,bool> tv(v.dims); \
        for (uint_t i : range(v)) { \
            tv.safe[i] = (u op v.safe[i]); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v, const vec<Dim,U>& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        vec<Dim,bool> tv(v.dims); \
        for (uint_t i : range(v)) { \
            tv.safe[i] = (v.safe[i] op u.safe[i]); \
        } \
        return tv; \
    }

VECTORIZE(==)
VECTORIZE(!=)
VECTORIZE(<)
VECTORIZE(<=)
VECTORIZE(>)
VECTORIZE(>=)

#undef VECTORIZE

template<typename T>
using is_bool_t = std::is_same<typename std::decay<typename std::remove_pointer<T>::type>::type, bool>;

#define VECTORIZE(op) \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        is_bool_t<T>::value && is_bool_t<U>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v1, const vec<Dim,U>& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        vec<Dim,bool> tv = v1; \
        for (uint_t i : range(v1)) { \
            tv.safe[i] = tv.safe[i] op v2.safe[i]; \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename U, typename enable = typename std::enable_if< \
        is_bool_t<U>::value>::type> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, const vec<Dim,U>& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i : range(v1)) { \
            v1.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v1, vec<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i : range(v2)) { \
            v2.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v2); \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, vec<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i : range(v1)) { \
            v1.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v1, bool b) { \
        vec<Dim,bool> tv = v1; \
        for (uint_t i : range(v1)) { \
            tv.safe[i] = tv.safe[i] op b; \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (bool b, const vec<Dim,T>& v2) { \
        vec<Dim,bool> tv = v2; \
        for (uint_t i : range(v2)) { \
            tv.safe[i] = b op tv.safe[i]; \
        } \
        return tv; \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, bool b) { \
        for (uint_t i : range(v1)) { \
            v1.safe[i] = v1.safe[i] op b; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (bool b, vec<Dim,bool>&& v2) { \
        for (uint_t i : range(v2)) { \
            v2.safe[i] = b op v2.safe[i]; \
        } \
        return std::move(v2); \
    }

VECTORIZE(&&)
VECTORIZE(||)

#undef VECTORIZE

template<std::size_t Dim, typename T, typename enable = typename std::enable_if<
    is_bool_t<T>::value>::type>
vec<Dim,bool> operator ! (const vec<Dim,T>& v) {
    vec<Dim,bool> tv = v;
    for (uint_t i : range(v)) {
        tv.safe[i] = !tv.safe[i];
    }

    return tv;
}

template<std::size_t Dim>
vec<Dim,bool> operator ! (vec<Dim,bool>&& v) {
    for (uint_t i : range(v)) {
        v.safe[i] = !v.safe[i];
    }

    return std::move(v);
}

// Print a vector into a stream.
template<typename O, std::size_t Dim, typename Type, typename enable = typename std::enable_if<!std::is_same<Type, bool>::value>::type>
O& operator << (O& o, const vec<Dim,Type>& v) {
    o << '{';
    for (uint_t i : range(v)) {
        if (i != 0) o << ", ";
        o << v.safe[i];
    }
    o << '}';

    return o;
}

template<typename O, std::size_t Dim>
O& operator << (O& o, const vec<Dim,bool>& v) {
    o << '{';
    for (uint_t i : range(v)) {
        if (i != 0) o << ", ";
        o << bool(v.safe[i]);
    }
    o << '}';

    return o;
}

#endif

