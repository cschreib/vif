#ifndef VEC_HPP
#define VEC_HPP

#include <tuple>
#include <vector>
#include <bitset>
#include <string>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include "color.hpp"
#include "variadic.hpp"
#include "print.hpp"
#include "iterator.hpp"

using int_t = std::ptrdiff_t;
using uint_t = std::size_t;
static const uint_t npos = uint_t(-1);

// Generic vector type
template<std::size_t Dim, typename Type>
struct vec_t;

// Helper to get the vector internal storage type.
// This is used to avoid std::vector<bool>.
template<typename T>
struct dtype_ {
    using type = T;
};
template<>
struct dtype_<bool> {
    using type = char;
};

template<typename T>
using dtype_t = typename dtype_<T>::type;

template<typename T>
struct vb_ref {
    static T& get(T& t) { return t; }
    static const T& get(const T& t) { return t; }
};
template<typename T>
struct vb_ref<const T> {
    static const T& get(const T& t) { return t; }
};
template<>
struct vb_ref<bool> {
    static bool& get(dtype_t<bool>& t) { return *(bool*)&t; }
    static const bool& get(const dtype_t<bool>& t) { return *(const bool*)&t; }
};
template<>
struct vb_ref<const bool> {
    static const bool& get(const dtype_t<bool>& t) { return *(const bool*)&t; }
};

template<typename T>
using rtype_t = typename std::remove_cv<typename std::remove_pointer<T>::type>::type;

// Helpers to reference/dereference variables only when necessary.
template<typename T>
T& dref(T& t) {
    return t;
}

template<typename T>
T& dref(T* t) {
    return *t;
}

template<typename T>
T* ref(T& t) {
    return &t;
}

template<typename T>
T* ref(T* t) {
    return t;
}

// Helper to check if a given type is a generic vector.
template<typename T>
struct is_vec : public std::false_type {};

template<std::size_t Dim, typename Type>
struct is_vec<vec_t<Dim,Type>> : public std::true_type {};

template<std::size_t Dim, typename Type>
struct is_vec<vec_t<Dim,Type>&> : public std::true_type {};

template<std::size_t Dim, typename Type>
struct is_vec<const vec_t<Dim,Type>&> : public std::true_type {};


// Return the data type of the provided type
template<typename T>
struct data_type_t_ {
    using type = T;
};

template<std::size_t Dim, typename T>
struct data_type_t_<vec_t<Dim,T>> {
    using type = T;
};

template<typename T>
using data_type_t = typename data_type_t_<T>::type;

// Helper to cound the number of generic vector in an argument list.
template<typename T, typename ... Args>
struct count_vec {
    static const std::size_t value = is_vec<T>::value + count_vec<Args...>::value;
};

template<typename T>
struct count_vec<T> {
    static const std::size_t value = is_vec<T>::value;
};

template<typename T>
struct make_vtype {
    using type = typename std::conditional<
        std::is_pointer<T>::value,
        T,
        typename std::remove_cv<T>::type
    >::type;
};

// Helper to build the result of v[ix(1,2,3)], i.e. when all indices are scalars.
// The result is a scalar.
template<std::size_t Dim, typename Type, typename ... Args>
struct make_vrtype_i {
    using type = dtype_t<typename std::remove_pointer<Type>::type>&;
    using vrtype = vec_t<Dim, typename make_vtype<Type>::type>;
    using vtype = typename std::conditional<std::is_const<Type>::value, const vrtype, vrtype>::type;

    static uint_t make_(vtype& vi, uint_t idx) {
        return idx;
    }

    template<typename T, typename ... Args2>
    static uint_t make_(vtype& v, uint_t idx, const T& ix, const Args2& ... i) {
        uint_t pitch = 1;
        for (uint_t j = Dim-sizeof...(Args2); j < Dim; ++j) {
            pitch *= v.dims[j];
        }

        idx += v.template to_idx<Dim-1-sizeof...(Args2)>(ix)*pitch;

        return make_(v, idx, i...);
    }

    static type make(vtype& v, const Args& ... i) {
        return reinterpret_cast<type>(dref(v.data[make_(v, 0, i...)]));
    }

    static type make(vtype& v, const std::tuple<Args...>& i) {
        return make_tuple_1(v, i, gen_seq_t<sizeof...(Args)>());
    }

    template<std::size_t ... S>
    static type make_tuple_1(vtype& v, const std::tuple<Args...>& i, seq_t<S...>) {
        return make(v, std::get<S>(i)...);
    }
};

template<typename Type, typename ... Args>
struct make_vrtype_i<1, Type, Args...> {
    using type = dtype_t<typename std::remove_pointer<Type>::type>&;
    using vrtype = vec_t<1, typename make_vtype<Type>::type>;
    using vtype = typename std::conditional<std::is_const<Type>::value, const vrtype, vrtype>::type;

    static type make(vtype& v, const std::tuple<Args...>& i) {
        return dref(v.data[v.to_idx(std::get<0>(i))]);
    }

    template<typename T>
    static type make(vtype& v, const T& i) {
        return dref(v.data[v.to_idx(i)]);
    }
};

template<std::size_t Dim, typename Type>
const vec_t<Dim,Type>& get_parent(const vec_t<Dim,Type>& v) {
    return v;
}

template<std::size_t Dim, typename Type>
void* get_parent(const vec_t<Dim,Type*>& v) {
    return v.parent;
}

// Helper to build the result of v[ix(_, rx(1,2), 5)], i.e. when at least one index is not scalar.
// The result is another array.
template<std::size_t Dim, typename Type, typename ... Args>
struct make_vrtype_v {
    static const std::size_t ODim = count_vec<decay_t<Args>...>::value + count_same<placeholder_t, decay_t<Args>...>::value;
    using type = vec_t<ODim, typename std::remove_pointer<Type>::type*>;
    using vrtype = vec_t<Dim, typename make_vtype<Type>::type>;
    using vtype = typename std::conditional<std::is_const<Type>::value, const vrtype, vrtype>::type;

    template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value && !std::is_same<placeholder_t,T>::value>::type>
    static void resize_(vtype& v, type& t, cte_t<ODim>, const T& ix) {
        t.resize();
    }

    static void resize_(vtype& v, type& t, cte_t<ODim-1>, placeholder_t) {
        t.dims[ODim-1] = v.dims[Dim-1];
        t.resize();
    }

    template<typename T>
    static void resize_(vtype& v, type& t, cte_t<ODim-1>, const vec_t<1,T>& r) {
        t.dims[ODim-1] = r.data.size();
        t.resize();
    }

    template<std::size_t M, typename T, typename ... Args2, typename enable = typename std::enable_if<!is_vec<T>::value && !std::is_same<placeholder_t,T>::value>::type>
    static void resize_(vtype& v, type& t, cte_t<M>, const T& ix, const Args2&... i) {
        resize_(v, t, cte_t<M>(), i...);
    }

    template<std::size_t M, typename ... Args2>
    static void resize_(vtype& v, type& t, cte_t<M>, placeholder_t, const Args2&... i) {
        t.dims[M] = v.dims[Dim-1-sizeof...(Args2)];
        resize_(v, t, cte_t<M+1>(), i...);
    }

    template<std::size_t M, typename T, typename ... Args2>
    static void resize_(vtype& v, type& t, cte_t<M>, const vec_t<1,T>& r, const Args2&... i) {
        t.dims[M] = r.data.size();
        resize_(v, t, cte_t<M+1>(), i...);
    }


    template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value && !std::is_same<placeholder_t,T>::value>::type>
    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, const T& ix) {
        t.data[itx] = ref(v.data[ivx+v.template to_idx<Dim-1>(ix)]);
        ++itx;
    }

    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, placeholder_t) {
        for (uint_t j = 0; j < v.dims[Dim-1]; ++j) {
            t.data[itx] = ref(v.data[ivx+j]);
            ++itx;
        }
    }

    template<typename T>
    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, const vec_t<1,T>& r) {
        for (uint_t j = 0; j < r.data.size(); ++j) {
            t.data[itx] = ref(v.data[ivx+v.template to_idx<Dim-1>(dref(r.data[j]))]);
            ++itx;
        }
    }

    template<typename T, typename ... Args2, typename enable = typename std::enable_if<!is_vec<T>::value && !std::is_same<placeholder_t,T>::value>::type>
    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, const T& ix, const Args2&... i) {
        make_(v, ivx +
            v.template to_idx<Dim-1-sizeof...(Args2)>(ix)*pitch[Dim-1-sizeof...(Args2)],
            pitch, t, itx, i...
        );
    }

    template<typename ... Args2>
    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, placeholder_t, const Args2&... i) {
        for (uint_t j = 0; j < v.dims[Dim-1-sizeof...(Args2)]; ++j) {
            make_(v, ivx + j*pitch[Dim-1-sizeof...(Args2)], pitch, t, itx, i...);
        }
    }

    template<typename T, typename ... Args2>
    static void make_(vtype& v, uint_t ivx, const std::array<uint_t, Dim>& pitch, type& t, uint_t& itx, const vec_t<1,T>& r, const Args2&... i) {
        for (uint_t j = 0; j < r.data.size(); ++j) {
            make_(v, ivx +
                v.template to_idx<Dim-1-sizeof...(Args2)>(dref(r.data[j]))*
                pitch[Dim-1-sizeof...(Args2)], pitch, t, itx, i...
            );
        }
    }

    static type make(vtype& v, const Args&... i) {
        type t(get_parent(v));
        resize_(v, t, cte_t<0>(), i...);
        t.resize();

        std::array<uint_t, Dim> pitch;
        for (uint_t j = 0; j < Dim; ++j) {
            pitch[j] = 1;
            for (uint_t k = j+1; k < Dim; ++k) {
                pitch[j] *= v.dims[k];
            }
        }

        uint_t idx = 0;
        make_(v, 0, pitch, t, idx, i...);
        return t;
    }

    static type make(vtype& v, const std::tuple<Args...>& i) {
        return make_tuple_(v, i, gen_seq_t<sizeof...(Args)>());
    }

    template<std::size_t ... S>
    static type make_tuple_(vtype& v, const std::tuple<Args...>& i, seq_t<S...>) {
        return make(v, std::get<S>(i)...);
    }
};

template<std::size_t Dim, typename Type, typename ... Args>
using make_vrtype = typename std::conditional<
    (count_vec<decay_t<Args>...>::value + count_same<placeholder_t, decay_t<Args>...>::value) == 0,
    make_vrtype_i<Dim, Type, Args...>,
    make_vrtype_v<Dim, Type, Args...>
>::type;

// Helper to intialize a generic vector from an initializer_list.
template<std::size_t N, std::size_t Dim, typename Type>
struct make_ilist_type {
    using type = std::initializer_list<typename make_ilist_type<N-1, Dim, Type>::type>;
};

template<std::size_t Dim, typename Type>
struct make_ilist_type<0, Dim, Type> {
    using type = std::initializer_list<Type>;
};

template<std::size_t Dim, typename Type>
struct ilist_t {
    using dtype = dtype_t<Type>;
    using type = typename make_ilist_type<Dim-1, Dim, dtype>::type;

    static void resize_(vec_t<Dim,Type>& v, const std::initializer_list<dtype>& il, cte_t<Dim-1>) {
        v.dims[Dim-1] = il.size();
        v.resize();
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void resize_(vec_t<Dim,Type>& v, const typename make_ilist_type<Dim-N-1, Dim-N, dtype>::type& il, cte_t<N>) {
        v.dims[N] = il.size();
        resize_(v, *il.begin(), cte_t<N+1>());
    }

    static void fill_(vec_t<Dim,Type>& v, const std::initializer_list<dtype>& il, uint_t& idx, cte_t<Dim-1>) {
        assert(il.size() == v.dims[Dim-1]);
        for (auto& t : il) {
            v.data[idx] = t;
            ++idx;
        }
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void fill_(vec_t<Dim,Type>& v, const typename make_ilist_type<Dim-N-1, Dim-N, dtype>::type& il, uint_t& idx, cte_t<N>) {
        assert(il.size() == v.dims[N]);
        for (auto& t : il) {
            fill_(v, t, idx, cte_t<N+1>());
        }
    }

    static void fill(vec_t<Dim,Type>& v, const type& il) {
        resize_(v, il, cte_t<0>());
        uint_t idx = 0;
        fill_(v, il, idx, cte_t<0>());
    }
};

template<std::size_t Dim, typename Type>
struct ilist_t<Dim, Type*> {
    using dtype = dtype_t<typename std::remove_cv<Type>::type>;
    using type = typename make_ilist_type<Dim-1, Dim, dtype>::type;

    static void fill_(vec_t<Dim,Type*>& v, const std::initializer_list<dtype>& il, uint_t& idx, cte_t<Dim-1>) {
        assert(il.size() == v.dims[Dim-1]);
        for (auto& t : il) {
            *v.data[idx] = t;
            ++idx;
        }
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void fill_(vec_t<Dim,Type*>& v, const typename make_ilist_type<Dim-N-1, Dim-N, dtype>::type& il, uint_t& idx, cte_t<N>) {
        assert(il.size() == v.dims[N]);
        for (auto& t : il) {
            fill_(v, t, idx, cte_t<N+1>());
        }
    }

    static void fill(vec_t<Dim,Type*>& v, const type& il) {
        uint_t idx = 0;
        fill_(v, il, idx, cte_t<0>());
    }
};

static bool do_print = false;

// The generic vector itself.
template<std::size_t Dim, typename Type>
struct vec_t {
    using effective_type = vec_t;
    using rtype = rtype_t<Type>;
    using dtype = dtype_t<Type>;
    using drtype = dtype_t<Type>;
    using vtype = std::vector<dtype>;
    using comparator = std::less<dtype>;

    vtype                        data;
    std::array<std::size_t, Dim> dims = {{0}};

    vec_t() = default;
    vec_t(const vec_t&) = default;
    vec_t(vec_t&&) = default;

    template<typename ... T>
    explicit vec_t(T ... d) {
        set_array(dims, d...);
        resize();
    }

    explicit vec_t(const std::array<std::size_t, Dim> d) : dims(d) {
        resize();
    }

    template<std::size_t N,
        typename enable = typename std::enable_if<std::is_same<Type, std::string>::value>::type>
    vec_t(const char t[N]) {
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 1;
        }
        resize();
        data[0] = t;
    }

    template<typename enable = typename std::enable_if<std::is_same<Type, std::string>::value>::type>
    vec_t(const char* t) {
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 1;
        }
        resize();
        data[0] = t;
    }

    vec_t(typename ilist_t<Dim, Type>::type t) {
        ilist_t<Dim, Type>::fill(*this, t);
    }

    template<typename T>
    vec_t(const vec_t<Dim,T>& v) : dims(v.dims) {
        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = dref(v.data[i]);
        }
    }

    vec_t& operator = (const Type& t) {
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 1;
        }
        resize();
        data[0] = t;
        return *this;
    }

    vec_t& operator = (typename ilist_t<Dim, Type>::type t) {
        ilist_t<Dim, Type>::fill(*this, t);
        return *this;
    }

    vec_t& operator = (const vec_t&) = default;
    vec_t& operator = (vec_t&&) = default;

    vec_t& operator = (const vec_t<Dim,Type*>& v) {
        if (v.is_same(*this)) {
            std::vector<dtype> t; t.resize(v.data.size());
            for (uint_t i = 0; i < v.data.size(); ++i) {
                t[i] = *v.data[i];
            }
            dims = v.dims;
            data.resize(v.data.size());
            for (uint_t i = 0; i < v.data.size(); ++i) {
                data[i] = t[i];
            }
        } else {
            dims = v.dims;
            data.resize(v.data.size());
            for (uint_t i = 0; i < v.data.size(); ++i) {
                data[i] = *v.data[i];
            }
        }

        return *this;
    }

    template<std::size_t D, typename T>
    bool is_same(const vec_t<D,T>&) const {
        return false;
    }

    template<std::size_t D>
    bool is_same(const vec_t<D,Type*>& v) const {
        return v.is_same(*this);
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

    void clear() {
        data.clear();
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 0;
        }
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


    template<typename enable = typename std::enable_if<(Dim > 1)>::type>
    void push_back(const vec_t<Dim-1,Type>& t) {
        if (empty()) {
            dims[0] = 1;
            for (uint_t i = 1; i < Dim; ++i) {
                dims[i] = t.dims[i-1];
            }

            data = t.data;
        } else {
            for (uint_t i = 0; i < Dim-1; ++i) {
                phypp_check(dims[i+1] == t.dims[i], "push_back: incompatible dimensions (", dims, " vs ", t.dims, ")");
            }

            data.insert(data.end(), t.data.begin(), t.data.end());
            ++dims[0];
        }
    }

    template<typename T, typename enable = typename std::enable_if<(Dim > 1)>::type>
    void push_back(const vec_t<Dim-1,T*>& t) {
        push_back(t.concretise());
    }

    void reserve(uint_t n) {
        data.reserve(n);
    }

    const vec_t<Dim,Type>& concretise() const {
        return *this;
    }

    vec_t<1, uint_t> ids(uint_t i) const {
        vec_t<1, uint_t> v(Dim);
        for (uint_t j = 0; j < Dim; ++j) {
            v[Dim-1-j] = i % dims[Dim-1-j];
            i /= dims[Dim-1-j];
        }

        return v;
    }

    template<typename T>
    T to_idx_(T ui, cte_t<false>) const {
        assert(ui < data.size());
        return ui;
    }

    template<typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += data.size();
        assert(i >= 0);
        uint_t ui(i);
        assert(ui < data.size());
        return ui;
    }

    template<std::size_t D, typename T>
    T to_idx_(T ui, cte_t<false>) const {
        assert(ui < dims[D]);
        return ui;
    }

    template<std::size_t D, typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += dims[D];
        assert(i >= 0);
        uint_t ui(i);
        assert(ui < dims[D]);
        return ui;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_<D>(ix, cte_t<std::is_signed<T>::value>());
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    Type& operator [] (T i) {
        return const_cast<Type&>(const_cast<const vec_t&>(*this)[i]);
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    const Type& operator [] (T i) const {
        return vb_ref<Type>::get(data[to_idx(i)]);
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,Type*> operator [] (const vec_t<1,T>& i) {
        vec_t<1,Type*> v(*this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = &vb_ref<Type>::get(data[to_idx(i[j])]);
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,Type*> operator [] (const vec_t<1,T*>& i) {
        vec_t<1,Type*> v(*this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = &vb_ref<Type>::get(data[to_idx(i[j])]);
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,const Type*> operator [] (const vec_t<1,T>& i) const {
        vec_t<1,const Type*> v(*this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = &vb_ref<Type>::get(data[to_idx(i[j])]);
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,const Type*> operator [] (const vec_t<1,T*>& i) const {
        vec_t<1,const Type*> v(*this);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = &vb_ref<Type>::get(data[to_idx(i[j])]);
        }
        return v;
    }

    vec_t<1,Type*> operator [] (placeholder_t) {
        vec_t<1,Type*> v(*this);
        v.data.resize(data.size());
        v.dims[0] = data.size();
        for (uint_t i = 0; i < data.size(); ++i) {
            v.data[i] = &vb_ref<Type>::get(data[i]);
        }
        return v;
    }

    vec_t<Dim,const Type*> operator [] (placeholder_t) const {
        vec_t<1,const Type*> v(*this);
        v.data.resize(data.size());
        v.dims[0] = data.size();
        for (uint_t i = 0; i < data.size(); ++i) {
            v.data[i] = &vb_ref<Type>::get(data[i]);
        }
        return v;
    }

    template<typename ... Args>
    auto operator [] (const std::tuple<Args...>& i) -> typename make_vrtype<Dim, Type, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, Type, Args...>::make(*this, i);
    }

    template<typename ... Args>
    auto operator [] (const std::tuple<Args...>& i) const -> typename make_vrtype<Dim, const Type, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, const Type, Args...>::make(*this, i);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) -> typename make_vrtype<Dim, Type, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, Type, Args...>::make(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const -> typename make_vrtype<Dim, const Type, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, const Type, Args...>::make(*this, i...);
    }

    vec_t operator + () const {
        return *this;
    }

    vec_t operator - () const {
        vec_t v = *this;
        for (auto& t : v) {
            t = -t;
        }

        return v;
    }

    #define OPERATOR(op) \
        template<typename U> \
        vec_t& operator op (const vec_t<Dim,U>& u) { \
            assert(data.size() == u.data.size()); \
            if (u.is_same(*this)) { \
                std::vector<dtype> t; t.resize(data.size()); \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    t[i] = dref(u.data[i]); \
                } \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    data[i] op t[i]; \
                } \
            } else { \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    data[i] op dref(u.data[i]); \
                } \
            } \
            return *this; \
        } \
        template<typename U> \
        vec_t& operator op (U u) { \
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

    using iterator = typename vtype::iterator;
    using const_iterator = typename vtype::const_iterator;

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

// "View" into a generic vector, allowing modification of indexed arrays: v[rx(1,3)] = indgen(3)
template<std::size_t Dim, typename Type>
struct vec_t<Dim,Type*> {
    using rtype = rtype_t<Type*>;
    using effective_type = vec_t<Dim,rtype>;
    using dtype = Type;
    using drtype = rtype;
    using vtype = std::vector<dtype*>;
    struct comparator {
        bool operator() (const dtype* t1, const dtype* t2) {
            return *t1 < *t2;
        }
        bool operator() (const dtype* t1, const dtype& t2) {
            return *t1 < t2;
        }
        bool operator() (const dtype& t1, const dtype* t2) {
            return t1 < *t2;
        }
        bool operator() (const dtype& t1, const dtype& t2) {
            return t1 < t2;
        }
    };

    void*                       parent = nullptr;
    vtype                       data;
    std::array<std::size_t,Dim> dims = {{0}};

    vec_t() = delete;
    vec_t(const vec_t&) = default;
    vec_t(vec_t&&) = default;

    template<std::size_t D>
    explicit vec_t(const vec_t<D,rtype>& p) {
        parent = static_cast<void*>(const_cast<vec_t<D,rtype>*>(&p));
    }

    explicit vec_t(void* p) {
        parent = p;
    }

    vec_t& operator = (const Type& t) {
        for (uint_t i = 0; i < data.size(); ++i) {
            *data[i] = t;
        }
        return *this;
    }

    vec_t& operator = (typename ilist_t<Dim, Type*>::type t) {
        ilist_t<Dim, Type*>::fill(*this, t);
        return *this;
    }

    vec_t& operator = (const vec_t<Dim,Type*>& v) {
        assert(data.size() == v.data.size());
        // Make a copy to prevent aliasing
        std::vector<dtype> t; t.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            t[i] = *v.data[i];
        }
        for (uint_t i = 0; i < v.data.size(); ++i) {
            *data[i] = t[i];
        }
        return *this;
    }

    template<typename T>
    vec_t& operator = (const vec_t<Dim,T>& v) {
        assert(data.size() == v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            *data[i] = dref(v.data[i]);
        }
        return *this;
    }

    template<std::size_t D, typename T>
    bool is_same(const vec_t<D,T>& v) const {
        return static_cast<void*>(const_cast<vec_t<D,T>*>(&v)) == parent;
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

    template<typename T>
    T to_idx_(T ui, cte_t<false>) const {
        assert(ui < data.size());
        return ui;
    }

    template<typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += data.size();
        assert(i >= 0);
        uint_t ui(i);
        assert(ui < data.size());
        return ui;
    }

    template<std::size_t D, typename T>
    T to_idx_(T ui, cte_t<false>) const {
        assert(ui < dims[D]);
        return ui;
    }

    template<std::size_t D, typename T>
    uint_t to_idx_(T i, cte_t<true>) const {
        if (i < 0) i += dims[D];
        assert(i >= 0);
        uint_t ui(i);
        assert(ui < dims[D]);
        return ui;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_<D>(ix, cte_t<std::is_signed<T>::value>());
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    Type& operator [] (T i) {
        return const_cast<Type&>(const_cast<const vec_t&>(*this)[i]);
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    const Type& operator [] (T i) const {
        return *data[to_idx(i)];
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,Type*> operator [] (const vec_t<1,T>& i) {
        vec_t<1,Type*> v(parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i[j])];
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,Type*> operator [] (const vec_t<1,T*>& i) {
        vec_t<1,Type*> v(parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i[j])];
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,const Type*> operator [] (const vec_t<1,T>& i) const {
        vec_t<1,const Type*> v(parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i[j])];
        }
        return v;
    }

    template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    vec_t<1,const Type*> operator [] (const vec_t<1,T*>& i) const {
        vec_t<1,const Type*> v(parent);
        v.data.resize(i.data.size());
        v.dims[0] = i.data.size();
        for (uint_t j = 0; j < i.data.size(); ++j) {
            v.data[j] = data[to_idx(i[j])];
        }
        return v;
    }

    vec_t<1,Type*> operator [] (placeholder_t) {
        vec_t<1,Type*> v(parent);
        v.data.resize(data.size());
        v.dims[0] = data.size();
        for (uint_t i = 0; i < data.size(); ++i) {
            v.data[i] = data[i];
        }
        return v;
    }

    vec_t<1,const Type*> operator [] (placeholder_t) const {
        vec_t<1,const Type*> v(parent);
        v.data.resize(data.size());
        v.dims[0] = data.size();
        for (uint_t i = 0; i < data.size(); ++i) {
            v.data[i] = data[i];
        }
        return v;
    }

    template<typename ... Args>
    auto operator [] (const std::tuple<Args...>& i) -> typename make_vrtype<Dim, Type*, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, Type*, Args...>::make(*this, i);
    }

    template<typename ... Args>
    auto operator [] (const std::tuple<Args...>& i) const -> typename make_vrtype<Dim, const Type*, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, const Type*, Args...>::make(*this, i);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) -> typename make_vrtype<Dim, Type*, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, Type*, Args...>::make(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const -> typename make_vrtype<Dim, const Type*, Args...>::type {
        static_assert(sizeof...(Args) == Dim, "wrong number of indices for this vector");
        return make_vrtype<Dim, const Type*, Args...>::make(*this, i...);
    }

    effective_type operator + () const {
        return *this;
    }

    effective_type operator - () const {
        effective_type v = *this;
        for (auto& t : v) {
            t = -t;
        }

        return *this;
    }

    #define OPERATOR(op) \
        template<typename U> \
        vec_t& operator op (const vec_t<Dim,U>& u) { \
            assert(data.size() == u.data.size()); \
            if (u.is_same(*this)) { \
                std::vector<dtype> t; t.resize(data.size()); \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    t[i] = dref(u.data[i]); \
                } \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    *data[i] op t[i]; \
                } \
            } else { \
                for (uint_t i = 0; i < data.size(); ++i) { \
                    *data[i] op dref(u.data[i]); \
                } \
            } \
            return *this; \
        } \
        \
        template<typename U> \
        vec_t& operator op (U u) { \
            for (auto& v : data) { \
                *v op u; \
            } \
            return *this; \
        }

    OPERATOR(*=)
    OPERATOR(/=)
    OPERATOR(+=)
    OPERATOR(-=)

    #undef OPERATOR

    using base_iterator = typename vtype::iterator;
    using base_const_iterator = typename vtype::const_iterator;

    using iterator = ptr_iterator_base<base_iterator, vec_t>;
    using const_iterator = const_ptr_iterator_base<base_const_iterator, vec_t>;

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

// Shortcuts for the most used types: vec1f, vec2d, ... up to vec6.
#define MAKE_TYPEDEFS(N) \
    using vec##N##f = vec_t<N, float>; \
    using vec##N##d = vec_t<N, double>; \
    using vec##N##i = vec_t<N, int_t>; \
    using vec##N##u = vec_t<N, std::size_t>; \
    using vec##N##s = vec_t<N, std::string>; \
    using vec##N##c = vec_t<N, char>; \
    using vec##N##b = vec_t<N, bool>; \
    using vec##N##rgb = vec_t<N, rgb>;

MAKE_TYPEDEFS(1)
MAKE_TYPEDEFS(2)
MAKE_TYPEDEFS(3)
MAKE_TYPEDEFS(4)
MAKE_TYPEDEFS(5)
MAKE_TYPEDEFS(6)

#undef MAKE_TYPEDEFS

// A few traits
template<typename T>
struct vec_dim;

template<std::size_t Dim, typename Type>
struct vec_dim<vec_t<Dim,Type>> {
    static const std::size_t value = Dim;
};

// Create vectors à la IDL.
template<typename T, typename ... Dims>
vec_t<sizeof...(Dims), T> arr(Dims ... ds) {
    vec_t<sizeof...(Dims), T> v;
    set_array(v.dims, ds...);
    v.resize();
    return v;
}

template<typename T, typename U, std::size_t N>
vec_t<N,T> arr(const std::array<U,N>& ds) {
    vec_t<N,T> v;
    v.dims = ds;
    v.resize();
    return v;
}

template<typename ... Dims>
auto fltarr(Dims ... ds) -> decltype(arr<float>(ds...)) {
    return arr<float>(ds...);
}

template<typename ... Dims>
auto dblarr(Dims ... ds) -> decltype(arr<double>(ds...)) {
    return arr<double>(ds...);
}

template<typename ... Dims>
auto intarr(Dims ... ds) -> decltype(arr<int_t>(ds...)) {
    return arr<int_t>(ds...);
}

template<typename ... Dims>
auto uintarr(Dims ... ds) -> decltype(arr<uint_t>(ds...)) {
    return arr<uint_t>(ds...);
}

template<typename ... Dims>
auto strarr(Dims ... ds) -> decltype(arr<std::string>(ds...)) {
    return arr<std::string>(ds...);
}

template<typename ... Dims>
auto bytarr(Dims ... ds) -> decltype(arr<char>(ds...)) {
    return arr<char>(ds...);
}

template<typename ... Dims>
auto boolarr(Dims ... ds) -> decltype(arr<bool>(ds...)) {
    return arr<bool>(ds...);
}

// Mathematical operators
struct op_mul_t;
struct op_div_t;
struct op_add_t;
struct op_sub_t;

struct op_node_t {
    op_mul_t operator * (op_node_t);
    op_div_t operator / (op_node_t);
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
using res_add_t = typename op_res_t<op_add_t, T, U>::type;
template<typename T, typename U>
using res_sub_t = typename op_res_t<op_sub_t, T, U>::type;

template<typename T>
const T& get_element_(const T& t, uint_t i) {
    return t;
}

template<std::size_t Dim, typename T>
auto get_element_(const vec_t<Dim,T>& t, uint_t i) {
    return dref(t.data[i]);
}

#define VECTORIZE(op, sop) \
    template<std::size_t Dim, typename T, typename U> \
    vec_t<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec_t<Dim,T>& v, const vec_t<Dim,U>& u) { \
        vec_t<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(get_element_(v, i) op get_element_(u, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U> \
    vec_t<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec_t<Dim,T>& v, const U& u) { \
        vec_t<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(get_element_(v, i) op u); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
    vec_t<Dim,T> operator op (vec_t<Dim,T>&& v, const vec_t<Dim,U>& u) { \
        for (uint_t i = 0; i < v.size(); ++i) { \
            v.data[i] sop get_element_(u, i); \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
    vec_t<Dim,T> operator op (vec_t<Dim,T>&& v, const U& u) { \
        for (auto& t : v) { \
            t sop u; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        !is_vec<U>::value>::type> \
    vec_t<Dim,typename op_res_t<OP_TYPE(op),U,T>::type> operator op (const U& u, const vec_t<Dim,T>& v) { \
        vec_t<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(u op get_element_(v, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
    vec_t<Dim,T> operator op (const vec_t<Dim,U>& u, vec_t<Dim,T>&& v) { \
        for (uint_t i = 0; i < v.size(); ++i) { \
            v.data[i] = get_element_(u, i) op v.data[i]; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
    vec_t<Dim,T> operator op (const U& u, vec_t<Dim,T>&& v) { \
        for (auto& t : v) { \
            t = u op t; \
        } \
        return std::move(v); \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value && \
        !std::is_pointer<U>::value>::type> \
    vec_t<Dim,T> operator op (vec_t<Dim,T>&& v, vec_t<Dim,U>&& u) { \
        for (uint_t i = 0; i < v.size(); ++i) { \
            v.data[i] sop u.data[i]; \
        } \
        return std::move(v); \
    } \

VECTORIZE(*, *=)
VECTORIZE(+, +=)
VECTORIZE(/, /=)
VECTORIZE(-, -=)

#undef VECTORIZE

// Logical operators
#define VECTORIZE(op) \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!is_vec<U>::value>::type> \
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v, const U& u) { \
        vec_t<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.data.size(); ++i) { \
            tv.data[i] = (dref(v.data[i]) op u); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!is_vec<U>::value>::type> \
    vec_t<Dim,bool> operator op (const U& u, const vec_t<Dim,T>& v) { \
        vec_t<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data[i] = (u op dref(v.data[i])); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U> \
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v, const vec_t<Dim,U>& u) { \
        assert(v.size() == u.size()); \
        vec_t<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data[i] = (dref(v.data[i]) op dref(u.data[i])); \
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
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v1, const vec_t<Dim,U>& v2) { \
        assert(v1.size() == v2.size()); \
        vec_t<Dim,bool> tv = v1; \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            tv.data[i] = tv.data[i] op dref(v2.data[i]); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename U, typename enable = typename std::enable_if< \
        is_bool_t<U>::value>::type> \
    vec_t<Dim,bool> operator op (vec_t<Dim,bool>&& v1, const vec_t<Dim,U>& v2) { \
        assert(v1.size() == v2.size()); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.data[i] = v1.data[i] op dref(v2.data[i]); \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v1, vec_t<Dim,bool>&& v2) { \
        assert(v1.size() == v2.size()); \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            v2.data[i] = dref(v1.data[i]) op v2.data[i]; \
        } \
        return std::move(v2); \
    } \
    template<std::size_t Dim> \
    vec_t<Dim,bool> operator op (vec_t<Dim,bool>&& v1, vec_t<Dim,bool>&& v2) { \
        assert(v1.size() == v2.size()); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.data[i] = v1.data[i] op dref(v2.data[i]); \
        } \
        return std::move(v1); \
    }

VECTORIZE(&&)
VECTORIZE(||)

#undef VECTORIZE

template<std::size_t Dim, typename T, typename enable = typename std::enable_if<
    is_bool_t<T>::value>::type>
vec_t<Dim,bool> operator ! (const vec_t<Dim,T>& v) {
    vec_t<Dim,bool> tv = v;
    for (uint_t i = 0; i < v.data.size(); ++i) {
        tv.data[i] = !tv.data[i];
    }

    return tv;
}

template<std::size_t Dim>
vec_t<Dim,bool> operator ! (vec_t<Dim,bool>&& v) {
    for (uint_t i = 0; i < v.data.size(); ++i) {
        v.data[i] = !v.data[i];
    }

    return std::move(v);
}

// Generate linearly increasing values.
template<typename T, typename A>
void indgen_vec_(std::vector<T,A>& v) {
    for (uint_t i = 0; i < v.size(); ++i) {
        v[i] = i;
    }
}

template<typename ... Dims>
auto findgen(Dims ... ds) -> vec_t<sizeof...(Dims), float> {
    auto v = fltarr(ds...);
    indgen_vec_(v.data);
    return v;
}

template<typename ... Dims>
auto dindgen(Dims ... ds) -> vec_t<sizeof...(Dims), double> {
    auto v = dblarr(ds...);
    indgen_vec_(v.data);
    return v;
}

template<typename ... Dims>
auto indgen(Dims ... ds) -> vec_t<sizeof...(Dims), int_t> {
    auto v = intarr(ds...);
    indgen_vec_(v.data);
    return v;
}

template<typename ... Dims>
auto uindgen(Dims ... ds) -> vec_t<sizeof...(Dims), uint_t> {
    auto v = uintarr(ds...);
    indgen_vec_(v.data);
    return v;
}

// Count the total number of elements in a vector.
template<std::size_t Dim, typename T>
uint_t n_elements(const vec_t<Dim,T>& v) {
    return v.data.size();
}

// Get the dimensions of a vector
template<std::size_t Dim, typename T>
vec1u dim(const vec_t<Dim,T>& v) {
    vec1u d = uintarr(Dim);
    for (uint_t i = 0; i < Dim; ++i) {
        d[i] = v.dims[i];
    }
    return d;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
uint_t n_elements(const T& v) {
    return 1;
}

template<typename T, typename U, typename ... Args>
bool same_size(const T& v1, const U& v2) {
    return n_elements(v1) && n_elements(v2);
}

template<typename T, typename U, typename ... Args>
bool same_size(const T& v1, const U& v2, const Args& ... args) {
    return n_elements(v1) && n_elements(v2) && same_size(v1, args...);
}

template<std::size_t Dim, typename T>
auto element(const vec_t<Dim,T>& v) -> decltype(v[0]) {
    return dref(v.data[0]);
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T& element(T& v) {
    return v;
}

template<std::size_t N, typename T>
std::array<T,N> rep(T t) {
    std::array<T,N> a;
    a.fill(t);
    return a;
}

// Create an single variable index from multiple indices to allow usage of operator[].
template<typename TP, typename T, typename ... Args>
auto ix__(TP&& tp, const T& t, cte_t<0>, Args&& ... args) {
    return ix_(std::move(tp), std::forward<Args>(args)...);
}

template<typename TP, typename T, std::size_t N, typename ... Args>
auto ix__(TP&& tp, const T& t, cte_t<N>, Args&& ... args) {
    return ix__(std::tuple_cat(std::move(tp), std::make_tuple(t)), t, cte_t<N-1>(), std::forward<Args>(args)...);
}

template<typename TP>
auto ix_(TP&& tp) {
    return std::move(tp);
}

template<typename TP, typename T, std::size_t N, typename ... Args>
auto ix_(TP&& tp, std::array<T,N>&& t, Args&& ... args) {
    return ix__(std::move(tp), t[0], cte_t<N>(), std::forward<Args>(args)...);
}

template<typename TP, typename T, typename ... Args>
auto ix_(TP&& tp, std::array<T,0>&& t, Args&& ... args) {
    return ix_(std::move(tp), std::forward<Args>(args)...);
}

template<typename TP, typename T, typename ... Args>
auto ix_(TP&& tp, T&& t, Args&& ... args) {
    return ix_(std::tuple_cat(std::move(tp), std::make_tuple(std::forward<T>(t))), std::forward<Args>(args)...);
}

template<typename ... Args>
auto ix(Args&& ... args) {
    return ix_(std::tuple<>(), std::forward<Args>(args)...);
}

// Create a range from a pair of indices.
template<typename T, typename U>
vec1u rx(T i, U j) {
    phypp_check(i >= 0 && j >= 0, "'rx' needs a positive or null values as arguments (got ", i, " and ", j, ")");

    if (i < T(j)) {
        uint_t n = j-i+1;
        vec1u v = uintarr(n);
        for (uint_t k = 0; k < n; ++k) {
            v[k] = i+k;
        }
        return v;
    } else {
        uint_t n = i-j+1;
        vec1u v = uintarr(n);
        for (uint_t k = 0; k < n; ++k) {
            v[k] = i-k;
        }
        return v;
    }
}

// Return the indices of the vector where the value is 1.
template<std::size_t Dim, typename Type>
vec1u where(const vec_t<Dim,Type>& v) {
    vec1u ids;
    ids.data.reserve(n_elements(v));
    uint_t i = 0;
    for (bool b : v) {
        if (b) {
            ids.data.push_back(i);
        }

        ++i;
    }

    ids.dims[0] = ids.data.size();
    ids.data.shrink_to_fit();
    return ids;
}

// In a sorted vector, return the first indices of each non unique sequence, effectively returning
// indices to all values that are different in the vector.
// By construction, the returned indices point to sorted values in the original vector.
template<std::size_t Dim, typename Type>
vec1u uniq(const vec_t<Dim,Type>& v) {
    vec1u r;
    if (v.empty()) return r;
    r.reserve(v.size()/4);

    rtype_t<Type> last = v[0];
    r.push_back(0);
    for (uint_t i = 1; i < v.size(); ++i) {
        if (v[i] != last) {
            r.push_back(i);
            last = v[i];
        }
    }

    r.data.shrink_to_fit();
    return r;
}

// In a vector, return indices to all values that are different. This version takes a second
// argument with indices that sort the input vector.
// The returned indices point to sorted values in the original vector.
template<std::size_t Dim, typename Type>
vec1u uniq(const vec_t<Dim,Type>& v, const vec1u& sid) {
    vec1u r;
    if (sid.empty()) return r;
    r.reserve(v.size()/4);

    rtype_t<Type> last = v[sid[0]];
    r.push_back(sid[0]);
    for (uint_t ti = 1; ti < sid.size(); ++ti) {
        uint_t i = sid[ti];
        if (v[i] != last) {
            r.push_back(i);
            last = v[i];
        }
    }

    r.data.shrink_to_fit();
    return r;
}

// For each value of the first vector, return 'true' if it is equal to any of the values of the
// second vector, and 'false' otherwise.
template<std::size_t Dim1, typename Type1, std::size_t Dim2 = Dim1, typename Type2 = Type1>
vec1b equal(const vec_t<Dim1,Type1>& v1, const vec_t<Dim2,Type2>& v2) {
    vec1b r(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r[i] = false;
        for (uint_t j = 0; j < v2.size(); ++j) {
            if (v1[i] == v2[j]) {
                r[i] = true;
            }
        }
    }

    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
void match(const vec_t<Dim,Type1>& v1, const vec_t<Dim,Type2>& v2, vec1u& id1, vec1u& id2) {
    uint_t n1 = v1.size();
    uint_t n2 = v2.size();
    if (n2 < n1) {
        match(v2, v1, id2, id1);
        return;
    }

    id1.data.reserve(n1);
    id2.data.reserve(n1);

    for (uint_t i = 0; i < n1; ++i) {
        vec1u r = where(v2 == v1[i]);
        if (r.empty()) continue;

        id1.data.push_back(i);
        id2.data.push_back(r[0]);
    }

    id1.dims[0] = id1.data.size();
    id1.data.shrink_to_fit();
    id2.dims[0] = id2.data.size();
    id2.data.shrink_to_fit();
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T flatten(T&& t) {
    return t;
}

template<std::size_t Dim, typename Type>
vec_t<1,Type> flatten(const vec_t<Dim,Type>& v) {
    vec_t<1,Type> r;
    r.dims[0] = v.data.size();
    r.data = v.data;
    return r;
}

template<std::size_t Dim, typename Type>
vec_t<1,Type> flatten(vec_t<Dim,Type>&& v) {
    vec_t<1,Type> r;
    r.dims[0] = v.data.size();
    r.data = std::move(v.data);
    return r;
}

template<std::size_t Dim, typename Type>
vec_t<1,Type*> flatten(const vec_t<Dim,Type*>& v) {
    vec_t<1,Type*> r(v.parent);
    r.dims[0] = v.data.size();
    r.data = v.data;
    return r;
}

template<std::size_t Dim, typename Type>
vec_t<1,Type*> flatten(vec_t<Dim,Type*>&& v) {
    vec_t<1,Type*> r(v.parent);
    r.dims[0] = v.data.size();
    r.data = std::move(v.data);
    return r;
}

template<typename Type>
vec_t<1,rtype_t<Type>> reverse(const vec_t<1,Type>& v) {
    vec_t<1,rtype_t<Type>> r = v;
    std::reverse(r.data.begin(), r.data.end());
    return r;
}

template<std::size_t Dim>
struct replicate_make_tuple_ {
    template<typename ... Args>
    static std::tuple<Args...> push_(const std::tuple<Args...>& t, cte_t<0>) {
        return t;
    }

    template<std::size_t N, typename ... Args>
    static auto push_(const std::tuple<Args...>& t, cte_t<N>) -> decltype(push_(std::tuple_cat(t, std::tuple<placeholder_t>()), cte_t<N-1>())) {
        return push_(std::tuple_cat(t, std::tuple<placeholder_t>()), cte_t<N-1>());
    }

    static auto make(const uint_t& i) -> decltype(push_(std::tuple<uint_t>(), cte_t<Dim>())) {
        std::tuple<uint_t> t;
        std::get<0>(t) = i;
        return push_(t, cte_t<Dim>());
    }
};

template<std::size_t Dim, typename Type, typename T>
vec_t<Dim+1,rtype_t<Type>> replicate(const vec_t<Dim,Type>& v, const T& n) {
    vec_t<Dim+1,rtype_t<Type>> r;
    r.dims[0] = n;
    for (uint_t i = 0; i < Dim; ++i) {
        r.dims[i+1] = v.dims[i];
    }
    r.resize();

    for (uint_t i = 0; i < std::size_t(n); ++i) {
        r[replicate_make_tuple_<Dim>::make(i)] = v;
    }

    return r;
}

template<std::size_t Dim, typename Type, typename T, typename ... Args>
vec_t<Dim+sizeof...(Args)+1,rtype_t<Type>> replicate(const vec_t<Dim,Type>& v, const Args& ... args, const T& n) {
    return replicate(replicate(v, args ...), n);
}

template<typename Type, typename ... Args, typename enable = typename std::enable_if<!is_vec<Type>::value>::type>
vec_t<sizeof...(Args),Type> replicate(const Type& v, const Args& ... args) {
    auto r = arr<Type>(args...);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<typename ... Args>
vec_t<sizeof...(Args),std::string> replicate(const char* v, const Args& ... args) {
    auto r = strarr(args...);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<std::size_t N, typename ... Args>
vec_t<sizeof...(Args),std::string> replicate(const char v[N], const Args& ... args) {
    auto r = strarr(args...);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<typename Type, std::size_t Dim>
vec_t<Dim,Type> replicate(const Type& v, const std::array<uint_t,Dim>& dims) {
    auto r = arr<Type>(dims);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<std::size_t Dim>
vec_t<Dim,std::string> replicate(const char* v, const std::array<uint_t,Dim>& dims) {
    auto r = strarr(dims);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<std::size_t N, std::size_t Dim>
vec_t<Dim,std::string> replicate(const char v[N], const std::array<uint_t,Dim>& dims) {
    auto r = strarr(dims);
    for (auto& t : r) {
        t = v;
    }
    return r;
}

template<std::size_t Dim, typename Type>
vec1u sort(const vec_t<Dim,Type>& v) {
    vec1u r = uindgen(v.size());
    std::sort(r.data.begin(), r.data.end(), [&v](uint_t i, uint_t j) {
        return typename vec_t<Dim,Type>::comparator()(v.data[i], v.data[j]);
    });

    return r;
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2, typename enable = typename std::enable_if<(N < Dim)>::type>
void append(vec_t<Dim,Type1>& t1, const vec_t<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
    } else {
        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1[ix(rep<N>(_),uindgen(n1),rep<Dim-N-1>(_))] = tmp;
        t1[ix(rep<N>(_),n1+uindgen(n2),rep<Dim-N-1>(_))] = t2;
    }
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2, typename enable = typename std::enable_if<(N < Dim)>::type>
void prepend(vec_t<Dim,Type1>& t1, const vec_t<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
    } else {
        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1[ix(rep<N>(_),uindgen(n2),rep<Dim-N-1>(_))] = t2;
        t1[ix(rep<N>(_),n2+uindgen(n1),rep<Dim-N-1>(_))] = tmp;
    }
}

template<typename Type1, typename Type2>
void append(vec_t<1,Type1>& t1, const vec_t<1,Type2>& t2) {
    t1.data.insert(t1.data.end(), t2.begin(), t2.end());
    t1.dims[0] += t2.dims[0];
}

template<typename Type1, typename Type2>
void prepend(vec_t<1,Type1>& t1, const vec_t<1,Type2>& t2) {
    t1.data.insert(t1.data.begin(), t2.begin(), t2.end());
    t1.dims[0] += t2.dims[0];
}

template<typename Type1, typename Type2>
vec_t<1,rtype_t<Type1>> merge(const vec_t<1,Type1>& t1, const vec_t<1,Type2>& t2) {
    using rtype = rtype_t<Type1>;
    vec_t<1,rtype> tv = t1;
    tv.data.insert(tv.data.end(), t2.begin(), t2.end());
    tv.dims[0] += t2.dims[0];
    return tv;
}

template<typename T, typename Type, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
vec_t<1,rtype_t<Type>> merge(T&& t, const vec_t<1,Type>& v) {
    using rtype = rtype_t<Type>;
    vec_t<1,rtype> tv = v;
    tv.data.insert(tv.data.begin(), std::forward<T>(t));
    ++tv.dims[0];
    return tv;
}

template<typename T, typename Type, typename enable = typename std::enable_if<!is_vec<T>::value && std::is_same<rtype_t<Type>, Type>::value>::type>
vec_t<1,Type> merge(T&& t, vec_t<1,Type>&& v) {
    v.data.insert(v.data.begin(), std::forward<T>(t));
    ++v.dims[0];
    return std::move(v);
}

template<typename T, typename Type, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
vec_t<1,rtype_t<Type>> merge(const vec_t<1,Type>& v, T&& t) {
    using rtype = rtype_t<Type>;
    vec_t<1,rtype> tv = v;
    tv.data.push_back(std::forward<T>(t));
    ++tv.dims[0];
    return tv;
}

template<typename T, typename Type, typename enable = typename std::enable_if<!is_vec<T>::value && std::is_same<rtype_t<Type>, Type>::value>::type>
vec_t<1,Type> merge(vec_t<1,Type>&& v, T&& t) {
    v.data.push_back(std::forward<T>(t));
    ++v.dims[0];
    return std::move(v);
}

// Print a vector into a stream.
template<typename O, std::size_t Dim, typename Type, typename enable = typename std::enable_if<!std::is_same<Type, bool>::value>::type>
O& operator << (O& o, const vec_t<Dim,Type>& v) {
    o << '[';
    for (uint_t i = 0; i < v.data.size(); ++i) {
        if (i != 0) o << ", ";
        o << v[i];
    }
    o << ']';

    return o;
}

template<typename O, std::size_t Dim>
O& operator << (O& o, const vec_t<Dim,bool>& v) {
    o << '[';
    for (uint_t i = 0; i < v.data.size(); ++i) {
        if (i != 0) o << ", ";
        o << bool(v[i]);
    }
    o << ']';

    return o;
}

#endif

