#ifndef VEC_HPP
#define VEC_HPP

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include "phypp/variadic.hpp"
#include "phypp/print.hpp"
#include "phypp/iterator.hpp"

using int_t = std::ptrdiff_t;
using uint_t = std::size_t;
static const uint_t npos = uint_t(-1);

template<typename T>
std::string strn(const T&);

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
struct is_vec_ : public std::false_type {};

template<std::size_t Dim, typename Type>
struct is_vec_<vec_t<Dim,Type>> : public std::true_type {};

template<typename T>
using is_vec = is_vec_<typename std::decay<T>::type>;

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

// Helper to match a type to the corresponding vector internal type.
template<typename T>
struct vec_type_ {
    using type = T;
};

template<>
struct vec_type_<char*> {
    using type = std::string;
};

template<>
struct vec_type_<const char*> {
    using type = std::string;
};

template<std::size_t N>
struct vec_type_<char[N]> {
    using type = std::string;
};

template<typename T>
using vec_type = typename vec_type_<typename std::decay<T>::type>::type;

namespace vec_access {
    template<std::size_t Dim, typename Type>
    const vec_t<Dim,Type>& get_parent(const vec_t<Dim,Type>& v) {
        return v;
    }

    template<std::size_t Dim, typename Type>
    void* get_parent(const vec_t<Dim,Type*>& v) {
        return v.parent;
    }

    // Generated vector dimension given index type
    template<typename T>
    struct output_dim_ :
        std::integral_constant<std::size_t, 0> {};

    template<std::size_t N, typename T>
    struct output_dim_<vec_t<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct output_dim_<placeholder_t> :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct output_dim_<repeated_value<N,T>> :
        std::integral_constant<std::size_t, N*output_dim_<T>::value> {};

    template<typename T>
    using output_dim = output_dim_<typename std::decay<T>::type>;

    template<typename T, typename ... Args>
    struct result_dim : std::integral_constant<std::size_t,
        output_dim<T>::value + result_dim<Args...>::value> {};

    template<typename T>
    struct result_dim<T> : std::integral_constant<std::size_t,
        output_dim<T>::value> {};

    // Needed vector dimension given index type
    template<typename T>
    struct input_dim_ :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct input_dim_<vec_t<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct input_dim_<placeholder_t> :
        std::integral_constant<std::size_t, 1> {};

    template<std::size_t N, typename T>
    struct input_dim_<repeated_value<N,T>> :
        std::integral_constant<std::size_t, N*input_dim_<T>::value> {};

    template<typename T>
    using input_dim = input_dim_<typename std::decay<T>::type>;

    template<typename T, typename ... Args>
    struct accessed_dim : std::integral_constant<std::size_t,
        input_dim<T>::value + accessed_dim<Args...>::value> {};

    template<typename T>
    struct accessed_dim<T> : std::integral_constant<std::size_t,
        input_dim<T>::value> {};

    // Helper to build the result of v(_, rgen(1,2), 5), i.e. when at least one index is not scalar.
    // The result is another array.
    template<bool IsConst, std::size_t Dim, std::size_t ODim, typename Type,
        typename ... Args>
    struct helper_ {
        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec_t<Dim,Type>, vec_t<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec_t<ODim, const rptype*>, vec_t<ODim, rptype*>>::type;

        // Functions to build the dimension of the resulting vector
        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type&, itype&, cte_t<IT>, cte_t<IV>, const T&) {}

        template<std::size_t IT, std::size_t IV>
        static void do_resize_(type& t, itype& v, cte_t<IT>, cte_t<IV>, placeholder_t) {
            t.dims[IT] = v.dims[IV];
        }

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type& t, itype&, cte_t<IT>, cte_t<IV>, const vec_t<1,T>& r) {
            t.dims[IT] = r.data.size();
        }

        template<std::size_t IT, std::size_t IV, typename T, typename ... Args2>
        static void resize_(type& t, itype& v, cte_t<IT> it, cte_t<IV> iv, const T& i,
            const Args2& ... args) {
            do_resize_(t, v, it, iv, i);
            resize_(t, v, cte_t<IT+output_dim<T>::value>(), cte_t<IV+input_dim<T>::value>(), args...);
        }

        static void resize_(type& t, itype& v, cte_t<ODim>, cte_t<Dim>) {
            t.resize();
        }

        // Functions to populate the resulting vector
        template<typename T>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, const T& ix) {

            t.data[itx] = ref(v.data[ivx+v.template to_idx<Dim-1>(ix)]);
            ++itx;
        }

        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, placeholder_t) {

            for (uint_t j = 0; j < v.dims[Dim-1]; ++j) {
                t.data[itx] = ref(v.data[ivx+j]);
                ++itx;
            }
        }

        template<typename T>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, const vec_t<1,T>& r) {

            for (uint_t j = 0; j < r.data.size(); ++j) {
                t.data[itx] = ref(v.data[ivx+v.template to_idx<Dim-1>(dref(r.data[j]))]);
                ++itx;
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, const T& ix, const Args2& ... i) {

            make_indices_(t, itx, v, ivx +
                v.template to_idx<IV>(ix)*pitch[IV], cte_t<IV+1>(), pitch, i...
            );
        }

        template<std::size_t IV, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, const placeholder_t&, const Args2& ... i) {

            for (uint_t j = 0; j < v.dims[IV]; ++j) {
                make_indices_(t, itx, v, ivx + j*pitch[IV], cte_t<IV+1>(), pitch, i...);
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, const vec_t<1,T>& r, const Args2& ... i) {

            for (uint_t j = 0; j < r.data.size(); ++j) {
                make_indices_(t, itx, v, ivx +
                    v.template to_idx<IV>(dref(r.data[j]))*pitch[IV], cte_t<IV+1>(), pitch, i...
                );
            }
        }

        template<typename ... UArgs>
        static type access_(itype& v, const UArgs& ... i) {
            type t(get_parent(v));
            resize_(t, v, cte_t<0>(), cte_t<0>(), i...);

            // TODO: cache this on construction
            std::array<uint_t, Dim> pitch;
            for (uint_t j = 0; j < Dim; ++j) {
                pitch[j] = 1;
                for (uint_t k = j+1; k < Dim; ++k) {
                    pitch[j] *= v.dims[k];
                }
            }

            uint_t itx = 0;
            make_indices_(t, itx, v, 0, cte_t<0>(), pitch, i...);
            return t;
        }

        template<typename ... UArgs>
        type operator() (itype& v, const UArgs& ... i) const {
            return access_(v, i...);
        }

        static type access(itype& v, const Args& ... i) {
            return unfold(helper_(), v, i...);
        }
    };

    template<bool IsConst, typename Type, typename Arg>
    struct helper_<IsConst, 1, 1, Type, Arg> {
        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec_t<1,Type>, vec_t<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec_t<1, const rptype*>, vec_t<1, rptype*>>::type;

        static type access(itype& v, placeholder_t) {
            type t(get_parent(v));
            t.dims[0] = v.dims[0];
            t.resize();

            for (uint_t i = 0; i < v.dims[0]; ++i) {
                t.data[i] = ref(v.data[i]);
            }

            return t;
        }

        template<typename T>
        static type access(itype& v, const vec_t<1,T>& r) {
            type t(get_parent(v));
            t.dims[0] = r.data.size();
            t.resize();

            for (uint_t i = 0; i < r.data.size(); ++i) {
                t.data[i] = ref(v.data[v.to_idx(dref(r.data[i]))]);
            }

            return t;
        }
    };

    // Helper to build the result of v(1,2,3), i.e. when all indices are scalars.
    // The result is a scalar.
    template<bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    struct helper_<IsConst, Dim, 0, Type, Args...>  {
        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec_t<Dim,Type>, vec_t<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            const dtype_t<rptype>&, dtype_t<rptype>&>::type;

        template<typename V>
        static uint_t get_index_(V& vi, uint_t idx) {
            return idx;
        }

        template<typename V, typename T, typename ... Args2>
        static uint_t get_index_(V& v, uint_t idx, const T& ix, const Args2& ... i) {
            uint_t pitch = 1;
            // TODO: cache this on construction
            for (uint_t j = Dim-sizeof...(Args2); j < Dim; ++j) {
                pitch *= v.dims[j];
            }

            idx += v.template to_idx<Dim-1-sizeof...(Args2)>(ix)*pitch;

            return get_index_(v, idx, i...);
        }

        static type access(itype& v, const Args& ... i) {
            return reinterpret_cast<type>(dref(v.data[get_index_(v, 0, i...)]));
        }
    };

    template<bool IsConst, typename Type, typename T>
    struct helper_<IsConst, 1, 0, Type, T> {
        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec_t<1,Type>, vec_t<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            const dtype_t<rptype>&, dtype_t<rptype>&>::type;

        static type access(itype& v, const T& i) {
            return dref(v.data[v.to_idx(i)]);
        }
    };

    template<bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    using helper = helper_<IsConst, Dim, result_dim<Args...>::value, Type, Args...>;
}

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

template<typename T>
struct is_dim_elem : std::is_arithmetic<T> {};


template<typename T, std::size_t N>
struct is_dim_elem<std::array<T,N>> : std::true_type {};

template<typename ... Args>
struct is_dim_list;

template<>
struct is_dim_list<> : std::true_type {};

template<typename T, typename ... Args>
struct is_dim_list<T, Args...> : std::integral_constant<bool,
    is_dim_elem<typename std::decay<T>::type>::value && is_dim_list<Args...>::value> {};

// The generic vector itself.
template<std::size_t Dim, typename Type>
struct vec_t;

template<typename Type>
struct vec_t<0,Type> {};

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

    template<typename T, typename ... Args, typename enable =
        typename std::enable_if<is_dim_list<T,Args...>::value>::type>
    explicit vec_t(T&& t, Args&& ... d) {
        set_array(dims, std::forward<T>(t), std::forward<Args>(d)...);
        resize();
    }

    template<std::size_t N>
    vec_t(const char t[N]) {
        static_assert(std::is_same<Type, std::string>::value,
            "character array constructor is only available for std::string vectors");
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 1;
        }
        resize();
        data[0] = t;
    }

    vec_t(const char* t) {
        static_assert(std::is_same<Type, std::string>::value,
            "character array constructor is only available for std::string vectors");
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

    template<typename ... Args>
    void resize(Args&& ... d) {
        set_array(dims, d...);
        resize();
    }

    void resize(const std::array<std::size_t, Dim> d) {
        dims = d;
        resize();
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

    template<typename T = Type, typename enable = typename std::enable_if<std::is_convertible<T,Type>::value>::type>
    void push_back(const vec_t<Dim-1,T>& t) {
        static_assert(Dim > 1, "cannot call push_back(vec_t<D-1>) on monodimensional vectors");
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

    template<typename T = Type, typename enable = typename std::enable_if<std::is_convertible<T,Type>::value>::type>
    void push_back(const vec_t<Dim-1,T*>& t) {
        static_assert(Dim > 1, "cannot call push_back(vec_t<D-1>) on monodimensional vectors");
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

    uint_t pitch(uint_t i) const {
        uint_t p = 1;
        for (uint_t j = i+1; j < Dim; ++j) {
            p *= dims[j];
        }
        return p;
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
    auto operator () (const Args& ... i) ->
        typename vec_access::helper<false, Dim, Type, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<false, Dim, Type, Args...>::access(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const ->
        typename vec_access::helper<true, Dim, Type, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, Dim, Type, Args...>::access(*this, i...);
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
            phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                "' ("+strn(dims)+" vs "+strn(u.dims)+")"); \
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

    uint_t pitch(uint_t i) const {
        uint_t p = 1;
        for (uint_t j = i+1; j < Dim; ++j) {
            p *= dims[j];
        }
        return p;
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
    auto operator () (const Args& ... i) ->
        typename vec_access::helper<false, Dim, Type*, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<false, Dim, Type*, Args...>::access(*this, i...);
    }

    template<typename ... Args>
    auto operator () (const Args& ... i) const ->
        typename vec_access::helper<true, Dim, Type*, Args...>::type {
        static_assert(vec_access::accessed_dim<Args...>::value == Dim,
            "wrong number of indices for this vector");
        return vec_access::helper<true, Dim, Type*, Args...>::access(*this, i...);
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
            phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                "' ("+strn(dims)+" vs "+strn(u.dims)+")"); \
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
    using vec##N##b = vec_t<N, bool>;

MAKE_TYPEDEFS(1)
MAKE_TYPEDEFS(2)
MAKE_TYPEDEFS(3)
MAKE_TYPEDEFS(4)
MAKE_TYPEDEFS(5)
MAKE_TYPEDEFS(6)

#undef MAKE_TYPEDEFS

// A few traits
template<typename T>
struct vec_dim_ : std::integral_constant<std::size_t, 0> {};

template<std::size_t Dim, typename Type>
struct vec_dim_<vec_t<Dim,Type>> : std::integral_constant<std::size_t, Dim> {};

template<typename T>
using vec_dim = vec_dim_<typename std::decay<T>::type>;

// Create vectors a la IDL.
template<typename T, typename ... Dims>
vec_t<dim_total<Dims...>::value, T> arr(Dims&& ... ds) {
    return vec_t<dim_total<Dims...>::value, T>(std::forward<Dims>(ds)...);
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
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
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
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
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
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
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
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
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
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
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
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        vec_t<Dim,bool> tv = v1; \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            tv.data[i] = tv.data[i] op dref(v2.data[i]); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename U, typename enable = typename std::enable_if< \
        is_bool_t<U>::value>::type> \
    vec_t<Dim,bool> operator op (vec_t<Dim,bool>&& v1, const vec_t<Dim,U>& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.data[i] = v1.data[i] op dref(v2.data[i]); \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v1, vec_t<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            v2.data[i] = dref(v1.data[i]) op v2.data[i]; \
        } \
        return std::move(v2); \
    } \
    template<std::size_t Dim> \
    vec_t<Dim,bool> operator op (vec_t<Dim,bool>&& v1, vec_t<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.data[i] = v1.data[i] op dref(v2.data[i]); \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec_t<Dim,bool> operator op (const vec_t<Dim,T>& v1, bool b) { \
        vec_t<Dim,bool> tv = v1; \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            tv.data[i] = tv.data[i] op b; \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec_t<Dim,bool> operator op (bool b, const vec_t<Dim,T>& v2) { \
        vec_t<Dim,bool> tv = v2; \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            tv.data[i] = b op tv.data[i]; \
        } \
        return tv; \
    } \
    template<std::size_t Dim> \
    vec_t<Dim,bool> operator op (vec_t<Dim,bool>&& v1, bool b) { \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.data[i] = v1.data[i] op b; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim> \
    vec_t<Dim,bool> operator op (bool b, vec_t<Dim,bool>&& v2) { \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            v2.data[i] = b op v2.data[i]; \
        } \
        return std::move(v2); \
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
template<typename T, typename ... Dims>
auto indgen_(Dims&& ... ds) -> decltype(arr<T>(std::forward<Dims>(ds)...)) {
    auto v = arr<T>(std::forward<Dims>(ds)...);
    for (uint_t i = 0; i < v.size(); ++i) {
        v[i] = i;
    }
    return v;
}

template<typename ... Dims>
auto findgen(Dims&& ... ds) -> decltype(indgen_<float>(std::forward<Dims>(ds)...)) {
    return indgen_<float>(std::forward<Dims>(ds)...);
}

template<typename ... Dims>
auto dindgen(Dims&& ... ds) -> decltype(indgen_<double>(std::forward<Dims>(ds)...)) {
    return indgen_<double>(std::forward<Dims>(ds)...);
}

template<typename ... Dims>
auto indgen(Dims&& ... ds) -> decltype(indgen_<int_t>(std::forward<Dims>(ds)...)) {
    return indgen_<int_t>(std::forward<Dims>(ds)...);
}

template<typename ... Dims>
auto uindgen(Dims&& ... ds) -> decltype(indgen_<uint_t>(std::forward<Dims>(ds)...)) {
    return indgen_<uint_t>(std::forward<Dims>(ds)...);
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

template<std::size_t Dim, typename Type, typename ... Args>
vec_t<sizeof...(Args), Type> reform(const vec_t<Dim,Type>& v, Args&& ... args) {
    auto r = arr<rtype_t<Type>>(std::forward<Args>(args)...);
    phypp_check(r.size() == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    for (uint_t i : range(v)) {
        r[i] = v[i];
    }

    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec_t<sizeof...(Args), Type> reform(vec_t<Dim,Type>&& v, Args&& ... args) {
    vec_t<sizeof...(Args), Type> r;
    set_array(r.dims, std::forward<Args>(args)...);
    std::size_t size = 1;
    for (uint_t i = 0; i < sizeof...(Args); ++i) {
        size *= r.dims[i];
    }

    phypp_check(size == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    r.data = std::move(v.data);

    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec_t<sizeof...(Args), Type*> reform(const vec_t<Dim,Type*>& v, Args&& ... args) {
    vec_t<sizeof...(Args), Type*> r(v.parent);
    r.resize(std::forward<Args>(args)...);
    phypp_check(r.size() == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    for (uint_t i : range(v)) {
        r[i] = v[i];
    }

    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec_t<sizeof...(Args), Type*> reform(vec_t<Dim,Type*>&& v, Args&& ... args) {
    vec_t<sizeof...(Args), Type*> r(v.parent);
    set_array(r.dims, std::forward<Args>(args)...);
    std::size_t size = 1;
    for (uint_t i = 0; i < sizeof...(Args); ++i) {
        size *= r.dims[i];
    }

    phypp_check(size == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    r.data = std::move(v.data);

    return r;
}

template<typename Type>
vec_t<1,Type> reverse(vec_t<1,Type> v) {
    std::reverse(v.data.begin(), v.data.end());
    return v;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec_t<Dim+dim_total<Args...>::value, rtype_t<Type>>
    replicate(const vec_t<Dim,Type>& t, Args&& ... args) {
    static const std::size_t FDim = Dim+dim_total<Args...>::value;
    vec_t<FDim, rtype_t<Type>> vec(std::forward<Args>(args)..., t.dims);

    std::size_t pitch = t.size();
    std::size_t n = vec.size()/pitch;
    for (uint_t i = 0; i < n; ++i) {
        for (uint_t j = 0; j < pitch; ++j) {
            vec.data[i*pitch + j] = t[j];
        }
    }

    return vec;
}

template<typename Type, typename ... Args>
vec_t<dim_total<Args...>::value, vec_type<Type>> replicate(const Type& t, Args&& ... args) {
    static const std::size_t FDim = dim_total<Args...>::value;
    vec_t<FDim, vec_type<Type>> vec(std::forward<Args>(args)...);

    for (auto& e : vec) {
        e = t;
    }

    return vec;
}

template<std::size_t Dim, typename Type>
vec1u sort(const vec_t<Dim,Type>& v) {
    vec1u r = uindgen(v.size());
    std::stable_sort(r.data.begin(), r.data.end(), [&v](uint_t i, uint_t j) {
        return typename vec_t<Dim,Type>::comparator()(v.data[i], v.data[j]);
    });

    return r;
}

template<std::size_t Dim, typename Type>
void inplace_sort(vec_t<Dim,Type>& v) {
    std::stable_sort(v.data.begin(), v.data.end(), typename vec_t<Dim,Type>::comparator());
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2, typename enable = typename std::enable_if<(N < Dim)>::type>
void append(vec_t<Dim,Type1>& t1, const vec_t<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
    } else {
        if (t2.empty()) return;

        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        phypp_check(t1.size()/n1 == t2.size()/n2, "cannot append dimension ", N, " in (", t1.dims,
            ") and (", t2.dims, ")");

        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1(repeat<N>(_), uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
        t1(repeat<N>(_), n1+uindgen(n2), repeat<Dim-N-1>(_)) = t2;
    }
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2, typename enable = typename std::enable_if<(N < Dim)>::type>
void prepend(vec_t<Dim,Type1>& t1, const vec_t<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
    } else {
        if (t2.empty()) return;

        std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
        phypp_check(t1.size()/n1 == t2.size()/n2, "cannot prepend dimension ", N, " in (", t1.dims,
            ") and (", t2.dims, ")");

        auto tmp = t1;
        t1.dims[N] += n2;
        t1.resize();

        t1(repeat<N>(_), uindgen(n2), repeat<Dim-N-1>(_)) = t2;
        t1(repeat<N>(_), n2+uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
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

template<typename Type>
vec_t<1,rtype_t<Type>> shift(const vec_t<1,Type>& v, int_t n, rtype_t<Type> def = 0) {
    vec_t<1,rtype_t<Type>> tmp;

    if (uint_t(abs(n)) > v.dims[0]) {
        tmp = replicate(def, v.dims[0]);
        return tmp;
    }

    tmp = v.concretise();

    if (n < 0) {
        tmp.data.erase(tmp.data.begin(), tmp.data.begin()-n);
        tmp.data.insert(tmp.data.end(), -n, def);
    } else if (n > 0) {
        tmp.data.erase(tmp.data.end()-n, tmp.data.end());
        tmp.data.insert(tmp.data.begin(), n, def);
    }

    return tmp;
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

