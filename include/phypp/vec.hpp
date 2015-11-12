#ifndef VEC_HPP
#define VEC_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <initializer_list>
#include "phypp/range.hpp"
#include "phypp/variadic.hpp"
#include "phypp/print.hpp"
#include "phypp/iterator.hpp"

template<typename T>
std::string strn(const T&);

template<typename T>
std::string pretty_type_(type_list<T>) {
    return typeid(T).name();
};

template<typename T>
std::string pretty_type_(type_list<T&>) {
    return pretty_type_(type_list<T>{})+"&";
};

template<typename T>
std::string pretty_type_(type_list<T*>) {
    return pretty_type_(type_list<T>{})+"*";
};

template<typename T>
std::string pretty_type_(type_list<const T>) {
    return "const "+pretty_type_(type_list<T>{});
};

std::string pretty_type_(type_list<char>) {
    return "char";
};

std::string pretty_type_(type_list<short>) {
    return "short";
};

std::string pretty_type_(type_list<int>) {
    return "int";
};

std::string pretty_type_(type_list<long>) {
    return "long";
};

std::string pretty_type_(type_list<long long>) {
    return "llong";
};

std::string pretty_type_(type_list<unsigned char>) {
    return "uchar";
};

std::string pretty_type_(type_list<unsigned short>) {
    return "ushort";
};

std::string pretty_type_(type_list<unsigned int>) {
    return "uint";
};

std::string pretty_type_(type_list<unsigned long>) {
    return "ulong";
};

std::string pretty_type_(type_list<unsigned long long>) {
    return "ullong";
};

std::string pretty_type_(type_list<bool>) {
    return "bool";
};

std::string pretty_type_(type_list<float>) {
    return "float";
};

std::string pretty_type_(type_list<double>) {
    return "double";
};

std::string pretty_type_(type_list<std::string>) {
    return "string";
};

template<std::size_t Dim, typename T>
std::string pretty_type_(type_list<vec<Dim,T>>) {
    return "vec<"+strn(Dim)+","+pretty_type_(type_list<T>{})+">";
};

#define pretty_type(x) pretty_type_(type_list<decltype(x)>{})
#define pretty_type_t(x) pretty_type_(type_list<x>{})

// Generic vector type
template<std::size_t Dim, typename Type>
struct vec;

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
using rtype_t = typename std::remove_cv<
    typename std::remove_pointer<
        typename std::remove_cv<T>::type
    >::type
>::type;

// Helpers to reference/dereference variables only when necessary.
// 'Type' must be the second template argument of the vector from which this
// value comes, i.e. vec<Dim,Type>.
template<typename Type, typename T>
rtype_t<Type>& dref(T& t) {
    return reinterpret_cast<rtype_t<Type>&>(t);
}
template<typename Type, typename T>
const rtype_t<Type>& dref(const T& t) {
    return reinterpret_cast<const rtype_t<Type>&>(t);
}

template<typename Type, typename T>
rtype_t<Type>& dref(T* t) {
    return reinterpret_cast<rtype_t<Type>&>(*t);
}
template<typename Type, typename T>
const rtype_t<Type>& dref(const T* t) {
    return reinterpret_cast<const rtype_t<Type>&>(*t);
}

template<typename Type, typename T>
rtype_t<Type>* ptr(T* t) {
    return reinterpret_cast<rtype_t<Type>*>(t);
}
template<typename Type, typename T>
const rtype_t<Type>* ptr(const T* t) {
    return reinterpret_cast<const rtype_t<Type>*>(t);
}

template<typename Type, typename T>
rtype_t<Type>* ptr(T& t) {
    return reinterpret_cast<rtype_t<Type>*>(&t);
}
template<typename Type, typename T>
const rtype_t<Type>* ptr(const T& t) {
    return reinterpret_cast<const rtype_t<Type>*>(&t);
}

// Iterator for bool vectors (necessary to avoid std::vector<bool>)
template<typename I>
struct bool_iterator_policy {
    static auto get_obj(I& i) -> decltype(dref<bool>(*i)) {
        return dref<bool>(*i);
    }
    static auto get_ptr(I& i) -> decltype(ptr<bool>(*i)) {
        return ptr<bool>(*i);
    }
};

template<typename T, typename C>
using bool_iterator_base = iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using const_bool_iterator_base = const_iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using reverse_bool_iterator_base = reverse_iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using const_reverse_bool_iterator_base = const_reverse_iterator_base<T,C,bool_iterator_policy<T>>;

namespace std {
    template<typename T, typename C>
    struct iterator_traits<bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = const bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<reverse_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_reverse_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = const bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };
}

// Helper to define the vector iterator types.
// Necessary because of tricks to avoid std::vector<bool>.
template<typename T>
struct vec_iterator_type {
    using vtype = typename T::vtype;
    using iterator = typename vtype::iterator;
    using const_iterator = typename vtype::const_iterator;
    using reverse_iterator = typename vtype::reverse_iterator;
    using const_reverse_iterator = typename vtype::const_reverse_iterator;
};

template<std::size_t Dim>
struct vec_iterator_type<vec<Dim,bool>> {
    using vtype = typename vec<Dim,bool>::vtype;
    using iterator = bool_iterator_base<typename vtype::iterator,vec<Dim,bool>>;
    using const_iterator = const_bool_iterator_base<typename vtype::const_iterator,vec<Dim,bool>>;
    using reverse_iterator = reverse_bool_iterator_base<
        typename vtype::reverse_iterator,vec<Dim,bool>>;
    using const_reverse_iterator = const_reverse_bool_iterator_base<
        typename vtype::const_reverse_iterator,vec<Dim,bool>>;
};

// Tag type to mark initialization of a reference vector.
struct vec_ref_tag_t {} vec_ref_tag;
// Tag type to permit copy initialization without data copy.
struct vec_nocopy_tag_t {} vec_nocopy_tag;

// Helper to check if a given type is a generic vector.
template<typename T>
struct is_vec_ : public std::false_type {};

template<std::size_t Dim, typename Type>
struct is_vec_<vec<Dim,Type>> : public std::true_type {};

template<typename T>
using is_vec = is_vec_<typename std::decay<T>::type>;

// Return the data type of the provided type
template<typename T>
struct data_type_t_ {
    using type = T;
};

template<std::size_t Dim, typename T>
struct data_type_t_<vec<Dim,T>> {
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
    const vec<Dim,Type>& get_parent(const vec<Dim,Type>& v) {
        return v;
    }

    template<std::size_t Dim, typename Type>
    void* get_parent(const vec<Dim,Type*>& v) {
        return v.parent;
    }

    // Generated vector dimension given index type
    template<typename T>
    struct output_dim_ :
        std::integral_constant<std::size_t, 0> {};

    template<std::size_t N, typename T>
    struct output_dim_<vec<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct output_dim_<full_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<left_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<right_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct output_dim_<left_right_range_t> :
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
    struct input_dim_<vec<N,T>> :
        std::integral_constant<std::size_t, N> {};

    template<>
    struct input_dim_<full_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<left_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<right_range_t> :
        std::integral_constant<std::size_t, 1> {};
    template<>
    struct input_dim_<left_right_range_t> :
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

    template<typename T>
    struct is_index_vector : std::false_type {};

    template<std::size_t Dim, typename T>
    struct is_index_vector<vec<Dim,T>> : std::integral_constant<bool,
        std::is_integral<typename std::decay<typename std::remove_pointer<T>::type>::type>::value> {};

    template<typename T>
    struct is_index_base : std::integral_constant<bool,
        std::is_integral<T>::value || is_range<T>::value || is_index_vector<T>::value> {};

    template<typename T>
    struct is_repeated_index : std::false_type {};

    template<std::size_t N, typename T>
    struct is_repeated_index<repeated_value<N,T>> : is_index_base<T> {};

    template<typename T>
    struct is_index : std::integral_constant<bool,
        is_index_base<T>::value || is_repeated_index<T>::value> {};

    template<typename ... Args>
    struct are_indices;

    template<>
    struct are_indices<> : std::true_type {};

    template<typename T, typename ... Args>
    struct are_indices<T, Args...> : std::integral_constant<bool,
        is_index<T>::value && are_indices<Args...>::value> {};

    // Helper to build the result of v(_, rgen(1,2), 5), i.e. when at least one index is not scalar.
    // The result is another array.
    template<bool IsSafe, bool IsConst, std::size_t Dim, std::size_t ODim, typename Type,
        typename ... Args>
    struct helper_ {
        static_assert(are_indices<Args...>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<Dim,Type>, vec<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec<ODim, const rptype*>, vec<ODim, rptype*>>::type;

        // Functions to build the dimension of the resulting vector
        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_impl_(type&, itype&, cte_t<IT>, cte_t<IV>, const T&, std::false_type) {}

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_impl_(type& t, itype& v, cte_t<IT>, cte_t<IV>, const T& rng, std::true_type) {
            t.dims[IT] = range_impl::range_size(rng, v.dims[IV]);
        }

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type& t, itype& v, cte_t<IT> d1, cte_t<IV> d2, const T& i) {
            do_resize_impl_(t, v, d1, d2, i, is_range<T>{});
        }

        template<std::size_t IT, std::size_t IV, typename T>
        static void do_resize_(type& t, itype&, cte_t<IT>, cte_t<IV>, const vec<1,T>& ids) {
            t.dims[IT] = ids.size();
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

        // Adapter to switch between safe/unsafe array indexing
        template<std::size_t D, typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.template to_idx<D>(t);
        }

        template<std::size_t D>
        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<std::size_t D, typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_<D>(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Functions to populate the resulting vector
        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, std::false_type, const T& ix) {

            t.data[itx] = ptr<Type>(v.data[ivx+to_idx<Dim-1>(v,ix)]);
            ++itx;
        }

        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, std::true_type, const T& rng) {

            for (uint_t j : range(rng, v.dims[Dim-1])) {
                t.data[itx] = ptr<Type>(v.data[ivx+j]);
                ++itx;
            }
        }

        template<typename T>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1> d,
            const std::array<uint_t, Dim>& pitch, const T& ix) {

            make_indices_(t, itx, v, ivx, d, pitch, is_range<T>{}, ix);
        }

        template<typename T>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<Dim-1>,
            const std::array<uint_t, Dim>& pitch, const vec<1,T>& ids) {

            for (uint_t j : ids) {
                t.data[itx] = ptr<Type>(v.data[ivx+to_idx<Dim-1>(v,j)]);
                ++itx;
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, std::false_type, const T& ix, const Args2& ... i) {

            make_indices_(t, itx, v, ivx +
                to_idx<IV>(v,ix)*pitch[IV], cte_t<IV+1>(), pitch, i...
            );
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_impl_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, std::true_type, const T& rng, const Args2& ... i) {

            for (uint_t j : range(rng, v.dims[IV])) {
                make_indices_(t, itx, v, ivx + j*pitch[IV], cte_t<IV+1>(), pitch, i...);
            }
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV> d,
            const std::array<uint_t, Dim>& pitch, const T& ix, const Args2& ... i) {

            make_indices_impl_(t, itx, v, ivx, d, pitch, is_range<T>{}, ix, i...);
        }

        template<std::size_t IV, typename T, typename ... Args2>
        static void make_indices_(type& t, uint_t& itx, itype& v, uint_t ivx, cte_t<IV>,
            const std::array<uint_t, Dim>& pitch, const vec<1,T>& ids, const Args2& ... i) {

            for (uint_t j : ids) {
                make_indices_(t, itx, v, ivx + to_idx<IV>(v,j)*pitch[IV], cte_t<IV+1>(),
                    pitch, i...
                );
            }
        }

        template<typename ... UArgs>
        static type access_(itype& v, const UArgs& ... i) {
            type t(vec_ref_tag, get_parent(v));
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

    template<bool IsSafe, bool IsConst, typename Type, typename Arg>
    struct helper_<IsSafe, IsConst, 1, 1, Type, Arg> {
        static_assert(is_index<Arg>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<1,Type>, vec<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            vec<1, const rptype*>, vec<1, rptype*>>::type;

        // Adapter to switch between safe/unsafe array indexing
        template<typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.to_idx(t);
        }

        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access functions
        template<typename T, typename enable = typename std::enable_if<is_range<T>::value>::type>
        static type access(itype& v, const T& rng) {
            type t(vec_ref_tag, get_parent(v));
            t.dims[0] = range_impl::range_size(rng, v.dims[0]);
            t.resize();

            uint_t itx = 0;
            for (uint_t i : range(rng, v.dims[0])) {
                t.data[itx] = ptr<Type>(v.data[i]);
                ++itx;
            }

            return t;
        }

        template<typename T>
        static type access(itype& v, const vec<1,T>& ids) {
            type t(vec_ref_tag, get_parent(v));
            t.dims[0] = ids.size();
            t.resize();

            uint_t itx = 0;
            for (uint_t i : ids) {
                t.data[itx] = ptr<Type>(v.data[to_idx(v,i)]);
                ++itx;
            }

            return t;
        }
    };

    // Helper to build the result of v(1,2,3), i.e. when all indices are scalars.
    // The result is a scalar.
    template<bool IsSafe, bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    struct helper_<IsSafe, IsConst, Dim, 0, Type, Args...>  {
        static_assert(are_indices<Args...>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<Dim,Type>, vec<Dim,Type>>::type;
        using type = typename std::conditional<IsConst,
            const dtype_t<rptype>&, dtype_t<rptype>&>::type;

        // Adapter to switch between safe/unsafe array indexing
        template<std::size_t D, typename T>
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.template to_idx<D>(t);
        }

        template<std::size_t D>
        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        template<std::size_t D, typename T>
        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_<D>(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access functions
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

            idx += to_idx<Dim-1-sizeof...(Args2)>(v, ix)*pitch;

            return get_index_(v, idx, i...);
        }

        static type access(itype& v, const Args& ... i) {
            return reinterpret_cast<type>(v.safe[get_index_(v, 0, i...)]);
        }
    };

    // Helper to build the result of v(1), i.e. when all indices are scalars.
    // The result is a scalar.
    template<bool IsSafe, bool IsConst, typename Type, typename T>
    struct helper_<IsSafe, IsConst, 1, 0, Type, T> {
        static_assert(is_index<T>::value, "vector access can only be made with "
            "integers or '_'");

        using rptype = typename std::remove_pointer<Type>::type;
        using itype = typename std::conditional<IsConst,
            const vec<1,Type>, vec<1,Type>>::type;
        using type = typename std::conditional<IsConst,
            const dtype_t<rptype>&, dtype_t<rptype>&>::type;

        // Adapter to switch between safe/unsafe array indexing
        static uint_t to_idx_(itype& v, const T& t, std::true_type) {
            return v.to_idx(t);
        }

        static uint_t to_idx_(itype&, uint_t t, std::false_type) {
            return t;
        }

        static uint_t to_idx(itype& v, const T& t) {
            return to_idx_(v, t, std::integral_constant<bool,IsSafe>{});
        }

        // Access function
        static type access(itype& v, const T& i) {
            return v.safe[to_idx(v, i)];
        }
    };

    template<bool IsSafe, bool IsConst, std::size_t Dim, typename Type, typename ... Args>
    using helper = helper_<IsSafe, IsConst, Dim, result_dim<Args...>::value, Type, Args...>;

    template<typename V, typename T>
    auto bracket_access(V& parent, const T& rng) ->
        vec<1,constify<typename V::rtype, V>*> {
        vec<1,constify<typename V::rtype, V>*> v(vec_ref_tag, parent);
        v.dims[0] = range_impl::range_size(rng, parent.size());
        v.data.resize(v.dims[0]);

        uint_t itx = 0;
        for (uint_t i : range(rng, parent.size())) {
            v.data[itx] = ptr<typename V::rtype>(parent.data[i]);
            ++itx;
        }

        return v;
    }
}

// Helper to intialize a generic vector from an initializer_list.
template<std::size_t D, typename T>
struct make_nested_initializer_list {
    using type = typename make_nested_initializer_list<D-1,std::initializer_list<T>>::type;
};

template<typename T>
struct make_nested_initializer_list<0,T> {
    using type = T;
};

template<std::size_t D, typename T>
using nested_initializer_list = typename make_nested_initializer_list<D,dtype_t<T>>::type;

template<std::size_t Dim, typename Type>
struct ilist_t {
    static void set_dim_(vec<Dim,Type>& v, nested_initializer_list<1,Type> il, cte_t<Dim-1>) {
        if (il.size() == 0) {
            for (uint_t i = 0; i < Dim; ++i) {
                v.dims[i] = 0;
            }
        } else {
            v.dims[Dim-1] = il.size();
        }
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void set_dim_(vec<Dim,Type>& v, nested_initializer_list<Dim-N,Type> il, cte_t<N>) {
        if (il.size() == 0) {
            for (uint_t i = 0; i < Dim; ++i) {
                v.dims[i] = 0;
            }
        } else {
            v.dims[N] = il.size();
            set_dim_(v, *il.begin(), cte_t<N+1>());
        }
    }

    static void fill_(vec<Dim,Type>& v, nested_initializer_list<1,Type> il, uint_t& idx, cte_t<Dim-1>) {
        phypp_check(il.size() == v.dims[Dim-1], "heterogeneous intializer lists are not allowed");
        for (auto& t : il) {
            v.data[idx] = t;
            ++idx;
        }
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void fill_(vec<Dim,Type>& v, nested_initializer_list<Dim-N,Type> il, uint_t& idx, cte_t<N>) {
        phypp_check(il.size() == v.dims[N], "heterogeneous intializer lists are not allowed");
        for (auto& t : il) {
            fill_(v, t, idx, cte_t<N+1>());
        }
    }

    static void fill(vec<Dim,Type>& v, nested_initializer_list<Dim,Type> il) {
        set_dim_(v, il, cte_t<0>());
        v.resize();
        if (!v.empty()) {
            uint_t idx = 0;
            fill_(v, il, idx, cte_t<0>());
        }
    }
};

template<std::size_t Dim, typename Type>
struct ilist_t<Dim, Type*> {
    static void fill_(vec<Dim,Type*>& v, nested_initializer_list<1,Type> il, uint_t& idx, cte_t<Dim-1>) {
        phypp_check(il.size() == v.dims[Dim-1], "heterogeneous intializer lists are not allowed");
        for (auto& t : il) {
            *v.data[idx] = t;
            ++idx;
        }
    }

    template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
    static void fill_(vec<Dim,Type*>& v, nested_initializer_list<Dim-N,Type> il, uint_t& idx, cte_t<N>) {
        phypp_check(il.size() == v.dims[N], "heterogeneous intializer lists are not allowed");
        for (auto& t : il) {
            fill_(v, t, idx, cte_t<N+1>());
        }
    }

    static void fill(vec<Dim,Type*>& v, nested_initializer_list<Dim,Type> il) {
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

template<typename TFrom, typename TTo>
struct vec_convertible_ : std::is_convertible<TFrom,TTo> {};

// Implicit conversion to/from bool is disabled
template<typename TFrom>
struct vec_convertible_<TFrom,bool> : std::is_same<TFrom,bool> {};

template<typename TFrom, typename TTo>
using vec_convertible = vec_convertible_<
    typename std::decay<typename std::remove_pointer<TFrom>::type>::type,
    typename std::decay<typename std::remove_pointer<TTo>::type>::type
>;

// The generic vector itself.
template<std::size_t Dim, typename Type>
struct vec;

template<typename Type>
struct vec<0,Type> {};

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

    vec() : safe(*this) {}
    vec(const vec& v) : data(v.data), dims(v.dims), safe(*this) {}
    vec(vec_nocopy_tag_t, const vec& v) : dims(v.dims), safe(*this) {}

    vec(vec&& v) : data(std::move(v.data)), dims(v.dims), safe(*this) {
        for (uint_t i = 0; i < Dim; ++i) {
            v.dims[i] = 0;
        }
    }

    template<typename ... Args, typename enable =
        typename std::enable_if<is_dim_list<Args...>::value>::type>
    explicit vec(Args&& ... d) : safe(*this) {
        set_array(dims, std::forward<Args>(d)...);
        resize();
    }

    vec(nested_initializer_list<Dim,Type> il) : safe(*this) {
        ilist_t<Dim, Type>::fill(*this, il);
    }

    template<typename T, typename enable = typename std::enable_if<!std::is_same<rtype_t<T>,bool>::value>::type>
    vec(const vec<Dim,T>& v) : dims(v.dims), safe(*this) {
        static_assert(vec_convertible<T,Type>::value, "could not assign vectors of "
            "non-convertible types");

        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = v.safe[i];
        }
    }

    template<typename T = Type, typename enable = typename std::enable_if<!std::is_same<T,bool>::value>::type>
    explicit vec(const vec<Dim,bool>& v) : dims(v.dims), safe(*this) {
        static_assert(vec_convertible<bool,Type>::value, "could not assign vectors of "
            "non-convertible types");
        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = v.safe[i];
        }
    }

    template<typename T = Type, typename enable = typename std::enable_if<!std::is_same<T,bool>::value>::type>
    explicit vec(const vec<Dim,bool*>& v) : dims(v.dims), safe(*this) {
        static_assert(vec_convertible<bool,Type>::value, "could not assign vectors of "
            "non-convertible types");
        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = v.safe[i];
        }
    }

    vec(const vec<Dim,Type*>& v) : dims(v.dims), safe(*this) {
        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = *v.data[i];
        }
    }

    vec& operator = (const Type& t) {
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = 1;
        }
        resize();
        data[0] = t;
        return *this;
    }

    vec& operator = (nested_initializer_list<Dim,Type> il) {
        ilist_t<Dim, Type>::fill(*this, il);
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

    template<typename T, typename enable = typename std::enable_if<!std::is_same<rtype_t<T>,bool>::value>::type>
    vec& operator = (const vec<Dim,T>& v) {
        static_assert(vec_convertible<T,Type>::value, "could not assign vectors of "
            "non-convertible types");

        dims = v.dims;
        data.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            data[i] = v.safe[i];
        }

        return *this;
    }

    vec& operator = (const vec<Dim,Type*>& v) {
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
    bool is_same(const vec<D,T>&) const {
        return false;
    }

    template<std::size_t D>
    bool is_same(const vec<D,Type*>& v) const {
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

    void clear() {
        data.clear();
        for (uint_t i = 0; i < Dim; ++i) {
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            if (u.is_same(*this)) { \
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

// Iterators for vector view
template<typename I>
struct ptr_iterator_policy {
    static auto get_obj(I& i) -> decltype(**i) {
        return **i;
    }
    static auto get_ptr(I& i) -> decltype(*i) {
        return *i;
    }
};

template<typename T, typename C>
using ptr_iterator_base = iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using const_ptr_iterator_base = const_iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using reverse_ptr_iterator_base = reverse_iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using const_reverse_ptr_iterator_base = const_reverse_iterator_base<T,C,ptr_iterator_policy<T>>;

namespace std {
    template<typename T, typename C>
    struct iterator_traits<ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<reverse_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_reverse_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };
}

// "View" into a generic vector, allowing modification of indexed arrays: v[rgen(1,3)] = indgen(3)
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

    vec& operator = (nested_initializer_list<Dim,Type> il) {
        ilist_t<Dim, Type*>::fill(*this, il);
        return *this;
    }

    vec& operator = (const vec<Dim,Type*>& v) {
        phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
            v.data.size(), " to ", data.size(), ")");

        // Make a copy to prevent aliasing
        // TODO: can this be optimized out if v.parent != this ?
        std::vector<dtype> t; t.resize(v.data.size());
        for (uint_t i = 0; i < v.data.size(); ++i) {
            t[i] = *v.data[i];
        }

        // Actual assignment
        for (uint_t i = 0; i < v.data.size(); ++i) {
            *data[i] = t[i];
        }

        return *this;
    }

    template<typename T, typename enable =
        typename std::enable_if<!std::is_same<rtype_t<T>,bool>::value || std::is_same<Type,bool>::value>::type>
    vec& operator = (const vec<Dim,T>& v) {
        static_assert(vec_convertible<T,Type>::value, "could not assign vectors of "
            "non-convertible types");
        phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
            v.data.size(), " to ", data.size(), ")");
        for (uint_t i = 0; i < v.data.size(); ++i) {
            *data[i] = v.safe[i];
        }
        return *this;
    }

    template<std::size_t D, typename T>
    bool is_same(const vec<D,T>& v) const {
        return static_cast<void*>(const_cast<vec<D,T>*>(&v)) == parent;
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
    typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
        return to_idx_(ix, cte_t<std::is_signed<T>::value>());
    }

    template<std::size_t D, typename T, typename enable =
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
    Type& operator [] (T i) {
        return const_cast<Type&>(const_cast<const vec&>(*this)[i]);
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
    const Type& operator [] (T i) const {
        return *data[to_idx(i)];
    }

    template<typename T, typename enable =
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
        typename std::enable_if<std::is_arithmetic<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            typename std::enable_if<std::is_integral<T>::value &&
                                    std::is_unsigned<T>::value>::type>
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
            if (u.is_same(*this)) { \
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

// A few traits
template<typename T>
struct vec_dim_ : std::integral_constant<std::size_t, 0> {};

template<std::size_t Dim, typename Type>
struct vec_dim_<vec<Dim,Type>> : std::integral_constant<std::size_t, Dim> {};

template<typename T>
using vec_dim = vec_dim_<typename std::decay<T>::type>;

// Create vectors a la IDL.
template<typename T, typename ... Dims>
vec<dim_total<Dims...>::value, T> arr(Dims&& ... ds) {
    return vec<dim_total<Dims...>::value, T>(std::forward<Dims>(ds)...);
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
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(get_element_(v, i) op get_element_(u, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U> \
    vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> operator op (const vec<Dim,T>& v, const U& u) { \
        vec<Dim,typename op_res_t<OP_TYPE(op),T,U>::type> tv; tv.dims = v.dims; tv.reserve(v.size()); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(get_element_(v, i) op u); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),T,U>::type, T>::value>::type> \
    vec<Dim,T> operator op (vec<Dim,T>&& v, const vec<Dim,U>& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        for (uint_t i = 0; i < v.size(); ++i) { \
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
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.data.push_back(u op get_element_(v, i)); \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if< \
        std::is_same<typename op_res_t<OP_TYPE(op),U,T>::type, T>::value>::type> \
    vec<Dim,T> operator op (const vec<Dim,U>& u, vec<Dim,T>&& v) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        for (uint_t i = 0; i < v.size(); ++i) { \
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
        for (uint_t i = 0; i < v.size(); ++i) { \
            v.data[i] sop u.data[i]; \
        } \
        return std::move(v); \
    } \

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
        vec<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.data.size(); ++i) { \
            tv.safe[i] = (v.safe[i] op u); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U, typename enable = typename std::enable_if<!is_vec<U>::value>::type> \
    vec<Dim,bool> operator op (const U& u, const vec<Dim,T>& v) { \
        vec<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            tv.safe[i] = (u op v.safe[i]); \
        } \
        return tv; \
    } \
    \
    template<std::size_t Dim, typename T, typename U> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v, const vec<Dim,U>& u) { \
        phypp_check(v.dims == u.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v.dims)+" vs "+strn(u.dims)+")"); \
        vec<Dim,bool> tv = boolarr(v.dims); \
        for (uint_t i = 0; i < v.size(); ++i) { \
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
        for (uint_t i = 0; i < v1.size(); ++i) { \
            tv.safe[i] = tv.safe[i] op v2.safe[i]; \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename U, typename enable = typename std::enable_if< \
        is_bool_t<U>::value>::type> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, const vec<Dim,U>& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v1, vec<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            v2.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v2); \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, vec<Dim,bool>&& v2) { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions in operator '" #op \
            "' ("+strn(v1.dims)+" vs "+strn(v2.dims)+")"); \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.safe[i] = v1.safe[i] op v2.safe[i]; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (const vec<Dim,T>& v1, bool b) { \
        vec<Dim,bool> tv = v1; \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            tv.safe[i] = tv.safe[i] op b; \
        } \
        return tv; \
    } \
    template<std::size_t Dim, typename T, typename enable = typename std::enable_if< \
        is_bool_t<T>::value>::type> \
    vec<Dim,bool> operator op (bool b, const vec<Dim,T>& v2) { \
        vec<Dim,bool> tv = v2; \
        for (uint_t i = 0; i < v2.size(); ++i) { \
            tv.safe[i] = b op tv.safe[i]; \
        } \
        return tv; \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (vec<Dim,bool>&& v1, bool b) { \
        for (uint_t i = 0; i < v1.size(); ++i) { \
            v1.safe[i] = v1.safe[i] op b; \
        } \
        return std::move(v1); \
    } \
    template<std::size_t Dim> \
    vec<Dim,bool> operator op (bool b, vec<Dim,bool>&& v2) { \
        for (uint_t i = 0; i < v2.size(); ++i) { \
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
    for (uint_t i = 0; i < v.data.size(); ++i) {
        tv.safe[i] = !tv.safe[i];
    }

    return tv;
}

template<std::size_t Dim>
vec<Dim,bool> operator ! (vec<Dim,bool>&& v) {
    for (uint_t i = 0; i < v.data.size(); ++i) {
        v.safe[i] = !v.safe[i];
    }

    return std::move(v);
}

// Generate linearly increasing values.
template<typename T, typename ... Dims>
auto indgen_(Dims&& ... ds) -> decltype(arr<T>(std::forward<Dims>(ds)...)) {
    auto v = arr<T>(std::forward<Dims>(ds)...);
    for (uint_t i = 0; i < v.size(); ++i) {
        v.safe[i] = i;
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
uint_t n_elements(const vec<Dim,T>& v) {
    return v.data.size();
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
uint_t n_elements(const T& v) {
    return 1;
}

// Get the dimensions of a vector
template<std::size_t Dim, typename T>
vec1u dim(const vec<Dim,T>& v) {
    vec1u d = uintarr(Dim);
    for (uint_t i = 0; i < Dim; ++i) {
        d.safe[i] = v.dims[i];
    }
    return d;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
vec1u dim(const T& t) {
    return {1u};
}

// Get multi-dim IDs from a flat ID
template<std::size_t D, typename T>
vec1u mult_ids(const vec<D,T>& v, uint_t i) {
    vec1u r(D);

    for (uint_t j : range(D)) {
        r.safe[D-1-j] = i % v.dims[D-1-j];
        i /= v.dims[D-1-j];
    }

    return r;
}

template<std::size_t D>
vec1u mult_ids(const std::array<uint_t,D>& dims, uint_t i) {
    vec1u r(D);

    for (uint_t j : range(D)) {
        r.safe[D-1-j] = i % dims[D-1-j];
        i /= dims[D-1-j];
    }

    return r;
}

// Get flat ID from multi-dim IDs
template<std::size_t D, typename T>
uint_t flat_id_(const vec<D,T>& v, uint_t ret, cte_t<D>) {
    return ret;
}

template<std::size_t D, typename T, std::size_t I, typename U, typename ... Args>
uint_t flat_id_(const vec<D,T>& v, uint_t ret, cte_t<I>, U i, Args&& ... args) {
    return flat_id_(v, ret + v.pitch(I)*v.template to_idx<I>(i), cte_t<I+1>{},
        std::forward<Args>(args)...);
}

template<std::size_t D, typename T, typename ... Args>
uint_t flat_id(const vec<D,T>& v, Args&& ... args) {
    static_assert(sizeof...(Args) == D, "wrong number of IDs provided");
    return flat_id_(v, 0, cte_t<0>{}, std::forward<Args>(args)...);
}


template<typename T, typename U, typename ... Args>
bool same_size(const T& v1, const U& v2) {
    return n_elements(v1) && n_elements(v2);
}

template<typename T, typename U, typename ... Args>
bool same_size(const T& v1, const U& v2, const Args& ... args) {
    return n_elements(v1) && n_elements(v2) && same_size(v1, args...);
}

bool same_dims_or_scalar_(uint_t size) {
    return true;
}

template<std::size_t Dim, typename T, typename ... Args>
bool same_dims_or_scalar_(uint_t size, const vec<Dim,T>& v1, const Args& ... args);

template<typename T, typename ... Args,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool same_dims_or_scalar_(uint_t size, const T& v1, const Args& ... args) {
    return same_dims_or_scalar_(size, args...);
}

template<std::size_t Dim, typename T, typename ... Args>
bool same_dims_or_scalar_(uint_t size, const vec<Dim,T>& v1, const Args& ... args) {
    return v1.size() == size && same_dims_or_scalar_(size, args...);
}

uint_t same_dims_or_scalar_get_size_() {
    return 0u;
}

template<std::size_t Dim, typename T, typename ... Args>
uint_t same_dims_or_scalar_get_size_(const vec<Dim,T>& v1, const Args& ... args) {
    return v1.size();
}

template<typename T, typename ... Args,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
uint_t same_dims_or_scalar_get_size_(const T& v1, const Args& ... args) {
    return same_dims_or_scalar_get_size_(args...);
}

template<typename ... Args>
bool same_dims_or_scalar(const Args& ... args) {
    uint_t size = same_dims_or_scalar_get_size_(args...);
    return size == 0 || same_dims_or_scalar_(size, args...);
}

template<std::size_t Dim, typename T>
auto element(const vec<Dim,T>& v) -> decltype(v.safe[0]) {
    phypp_check(!v.empty(), "cannot get element of empty array");
    return v.safe[0];
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T& element(T& v) {
    return v;
}

template<typename T>
auto first(const vec<1,T>& v) -> decltype(v.safe[0]) {
    phypp_check(!v.empty(), "cannot get first element of empty array");
    return v.safe[0];
}

template<typename T>
auto last(const vec<1,T>& v) -> decltype(v.safe[0]) {
    phypp_check(!v.empty(), "cannot get last element of empty array");
    return v.safe[v.data.size()-1];
}

// Return the indices of the vector where the value is 'true'.
template<std::size_t Dim, typename Type, typename enable =
    typename std::enable_if<std::is_same<rtype_t<Type>,bool>::value>::type>
vec1u where(const vec<Dim,Type>& v) {
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
    return ids;
}

// Build the complement of a set of indices, i.e. return the indices that are not part of the
// provided set.
template<std::size_t Dim, typename Type>
vec1u complement(const vec<Dim,Type>& v, const vec1u& ids) {
    if (ids.size() == v.size()) return vec1u();
    phypp_check(ids.size() < v.size(), "incompatible size of ids (", ids.size(),
        " vs. ", v.size(), ")");

    vec1b sel(v.size());
    for (uint_t i : ids) {
        sel[i] = true;
    }

    vec1u res; res.reserve(v.size() - ids.size());
    for (uint_t i = 0; i < v.size(); ++i) {
        if (!sel.safe[i]) {
            res.push_back(i);
        }
    }

    return res;
}

// In a sorted vector, return the first indices of each non unique sequence, effectively returning
// indices to all values that are different in the vector.
// By construction, the returned indices point to sorted values in the original vector.
template<std::size_t Dim, typename Type>
vec1u uniq(const vec<Dim,Type>& v) {
    vec1u r;
    if (v.empty()) return r;
    r.reserve(v.size()/4);

    rtype_t<Type> last = v.safe[0];
    r.push_back(0);
    for (uint_t i = 1; i < v.size(); ++i) {
        if (v.safe[i] != last) {
            r.push_back(i);
            last = v.safe[i];
        }
    }

    r.data.shrink_to_fit();
    return r;
}

// In a vector, return indices to all values that are different. This version takes a second
// argument with indices that sort the input vector.
// The returned indices point to sorted values in the original vector.
template<std::size_t Dim, typename Type>
vec1u uniq(const vec<Dim,Type>& v, const vec1u& sid) {
    vec1u r;
    if (sid.empty()) return r;
    r.reserve(v.size()/4);

    rtype_t<Type> last = v[sid[0]];
    r.push_back(sid[0]);
    for (uint_t ti = 1; ti < sid.size(); ++ti) {
        uint_t i = sid.safe[ti];
        if (v[i] != last) {
            r.push_back(i);
            last = v.safe[i];
        }
    }

    r.data.shrink_to_fit();
    return r;
}

// For each value of the first vector, return 'true' if it is equal to any of the values of the
// second vector, and 'false' otherwise.
template<typename Type1, std::size_t Dim2, typename Type2 = Type1>
bool is_any_of(const Type1& v1, const vec<Dim2,Type2>& v2) {
    for (uint_t j = 0; j < v2.size(); ++j) {
        if (v1 == v2.safe[j]) {
            return true;
        }
    }

    return false;
}

// For each value of the first vector, return 'true' if it is equal to any of the values of the
// second vector, and 'false' otherwise.
template<std::size_t Dim1, typename Type1, std::size_t Dim2 = Dim1, typename Type2 = Type1>
vec<Dim1,bool> is_any_of(const vec<Dim1,Type1>& v1, const vec<Dim2,Type2>& v2) {
    vec<Dim1,bool> r(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i)
    for (uint_t j = 0; j < v2.size(); ++j) {
        if (v1.safe[i] == v2.safe[j]) {
            r.safe[i] = true;
            break;
        }
    }

    return r;
}

// Compare the two provided vectors and push indices where the two match into 'id1' and 'id2'.
// Each value is matched once, and cannot be used again. For example, if the content of the
// two vectors is: [12,12,-1,-1] and [0,12,-1,5], then the function will return [0,2] and [1,2].
template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
void match(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2, vec1u& id1, vec1u& id2) {
    uint_t n1 = v1.size();
    uint_t n2 = v2.size();
    if (n2 < n1) {
        match(v2, v1, id2, id1);
        return;
    }

    id1.data.reserve(n1 + id1.size());
    id2.data.reserve(n1 + id2.size());

    for (uint_t i = 0; i < n1; ++i) {
        vec1u r = where(v2 == v1.safe[i]);
        if (r.empty()) continue;

        id1.data.push_back(i);
        id2.data.push_back(r.safe[0]);
    }

    id1.dims[0] = id1.data.size();
    id1.data.shrink_to_fit();
    id2.dims[0] = id2.data.size();
    id2.data.shrink_to_fit();
}

// Compare the two provided vectors and push indices where the two match into 'id1' and 'id2'.
// Indices are returned for every matching pair (assuming no repetitions in 'v2'). Using
// the sample example as 'match', with the two vectors containing [12,12,-1,-1] and
// [0,12,-1,5], the function will return [0,1,2,3] and [1,2,1,2]. Contrary to 'match',
// this function is not symetric, in the sense that the result will be different if the two
// input vectors are swapped. The order of the returned indices is unspecified.
template<std::size_t D1, std::size_t D2, typename Type1, typename Type2>
void match_dictionary(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2, vec1u& id1, vec1u& id2) {
    vec1u ids = sort(v1);
    auto sv1 = v1[ids].concretise();
    vec1u idu = uniq(sv1);

    vec1u tid1, tid2;
    match(sv1.safe[idu], v2, tid1, tid2);

    const uint_t nm = tid1.size();
    const uint_t nplus = nm*(float(idu.size())/v1.size());
    id1.reserve(id1.size() + nplus);
    id2.reserve(id2.size() + nplus);

    for (uint_t tu = 0; tu < nm; ++tu) {
        uint_t iu = tid1.safe[tu];
        uint_t u = idu.safe[iu];
        uint_t n = iu == idu.size()-1 ? sv1.size() - u : idu.safe[iu+1] - u;
        for (uint_t i = 0; i < n; ++i) {
            id1.push_back(ids.safe[u+i]);
            id2.push_back(tid2.safe[tu]);
        }
    }
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T flatten(T&& t) {
    return t;
}

template<std::size_t Dim, typename Type>
vec<1,Type> flatten(const vec<Dim,Type>& v) {
    vec<1,Type> r;
    r.dims[0] = v.data.size();
    r.data = v.data;
    return r;
}

template<std::size_t Dim, typename Type>
vec<1,Type> flatten(vec<Dim,Type>&& v) {
    vec<1,Type> r;
    r.dims[0] = v.data.size();
    r.data = std::move(v.data);
    return r;
}

template<std::size_t Dim, typename Type>
vec<1,Type*> flatten(const vec<Dim,Type*>& v) {
    vec<1,Type*> r(vec_ref_tag, v.parent);
    r.dims[0] = v.data.size();
    r.data = v.data;
    return r;
}

template<std::size_t Dim, typename Type>
vec<1,Type*> flatten(vec<Dim,Type*>&& v) {
    vec<1,Type*> r(vec_ref_tag, v.parent);
    r.dims[0] = v.data.size();
    r.data = std::move(v.data);
    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec<sizeof...(Args), Type> reform(const vec<Dim,Type>& v, Args&& ... args) {
    auto r = arr<rtype_t<Type>>(std::forward<Args>(args)...);
    phypp_check(r.size() == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    for (uint_t i : range(v)) {
        r.safe[i] = v.safe[i];
    }

    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec<sizeof...(Args), Type> reform(vec<Dim,Type>&& v, Args&& ... args) {
    vec<sizeof...(Args), Type> r;
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
vec<sizeof...(Args), Type*> reform(const vec<Dim,Type*>& v, Args&& ... args) {
    vec<sizeof...(Args), Type*> r(vec_ref_tag, v.parent);
    r.resize(std::forward<Args>(args)...);
    phypp_check(r.size() == v.size(),
        "incompatible dimensions ("+strn(v.dims)+" vs "+strn(r.dims)+")");

    for (uint_t i : range(v)) {
        r.safe[i] = v.safe[i];
    }

    return r;
}

template<std::size_t Dim, typename Type, typename ... Args>
vec<sizeof...(Args), Type*> reform(vec<Dim,Type*>&& v, Args&& ... args) {
    vec<sizeof...(Args), Type*> r(vec_ref_tag, v.parent);
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
vec<1,Type> reverse(vec<1,Type> v) {
    std::reverse(v.data.begin(), v.data.end());
    return v;
}

template<typename Type>
vec<2,Type> transpose(const vec<2,Type>& v) {
    vec<2,Type> r(vec_nocopy_tag, v);
    std::swap(r.dims[0], r.dims[1]);

    for (uint_t i : range(v.size())) {
        r.data.push_back(v.data[(i%v.dims[0])*v.dims[1] + i/v.dims[0]]);
    }

    // TODO: see who's faster
    // r.resize();
    // for (uint_t i : range(r.dims[0]))
    // for (uint_t j : range(r.dims[1])) {
    //     r.data[j+i*r.dims[1]] = v.data[i+j*v.dims[1]];
    // }

    return r;
}

template<std::size_t Dim, typename Type = double, typename ... Args>
vec<Dim+dim_total<Args...>::value, rtype_t<Type>>
    replicate(const vec<Dim,Type>& t, Args&& ... args) {
    static const std::size_t FDim = Dim+dim_total<Args...>::value;
    vec<FDim, rtype_t<Type>> v(std::forward<Args>(args)..., t.dims);

    std::size_t pitch = t.size();
    std::size_t n = v.size()/pitch;
    for (uint_t i = 0; i < n; ++i) {
        for (uint_t j = 0; j < pitch; ++j) {
            v.safe[i*pitch + j] = t.safe[j];
        }
    }

    return v;
}

template<typename Type, typename ... Args>
vec<dim_total<Args...>::value, vec_type<Type>> replicate(const Type& t, Args&& ... args) {
    static const std::size_t FDim = dim_total<Args...>::value;
    vec<FDim, vec_type<Type>> v(std::forward<Args>(args)...);

    for (auto& e : v) {
        e = t;
    }

    return v;
}

template<std::size_t Dim, typename Type>
vec1u sort(const vec<Dim,Type>& v) {
    vec1u r = uindgen(v.size());
    std::stable_sort(r.data.begin(), r.data.end(), [&v](uint_t i, uint_t j) {
        return typename vec<Dim,Type>::comparator()(v.data[i], v.data[j]);
    });

    return r;
}

template<std::size_t Dim, typename Type>
void inplace_sort(vec<Dim,Type>& v) {
    std::stable_sort(v.data.begin(), v.data.end(), typename vec<Dim,Type>::comparator());
}

// Check if a given array is sorted or not
template<std::size_t Dim, typename Type>
bool is_sorted(const vec<Dim,Type>& v) {
    for (uint_t i = 0; i < v.size()-1; ++i) {
        if (v.safe[i] >= v.safe[i+1]) return false;
    }

    return true;
}

// Returns the position of the last value in the array that is less than or equal to 'x'.
// Returns 'npos' if no value satisfy this criterium.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, std::size_t Dim, typename Type>
uint_t lower_bound(T x, const vec<Dim,Type>& v) {
    if (v.empty()) return npos;

    auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
        typename vec<Dim,Type>::comparator());

    if (iter == v.data.begin()) {
        return npos;
    } else {
        return iter - v.data.begin() - 1;
    }
}

// Returns the position of the first value in the array that is greater than 'x'.
// Returns 'npos' if no value satisfy this criterium.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, std::size_t Dim, typename Type>
uint_t upper_bound(T x, const vec<Dim,Type>& v) {
    if (v.empty()) return npos;

    auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
        typename vec<Dim,Type>::comparator());

    if (iter == v.data.end()) {
        return npos;
    } else {
        return iter - v.data.begin();
    }
}

// Return the position of the last value in 'v' that is less than or equal to 'x' and
// the position of the first value in 'v' that is greater than 'x'.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, std::size_t Dim, typename Type>
std::array<uint_t,2> bounds(T x, const vec<Dim,Type>& v) {
    std::array<uint_t,2> res;

    if (v.empty()) {
        res[0] = npos;
        res[1] = npos;
    } else {
        auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
            typename vec<Dim,Type>::comparator());

        if (iter == v.data.end()) {
            res[0] = v.size() - 1;
            res[1] = npos;
        } else if (iter == v.data.begin()) {
            res[0] = npos;
            res[1] = 0;
        } else {
            res[1] = iter - v.data.begin();
            res[0] = res[1] - 1;
        }
    }

    return res;
}

// Return the position of the last value in 'v' that is less than or equal to 'x1' and
// the position of the first value in 'v' that is greater than 'x2'.
// Note: assumes that:
//  1) 'v' is sorted and does not contain NaN values,
//  2) 'x2' is greater than or equal to 'x1'.
template<typename T, typename U, std::size_t Dim, typename Type>
std::array<uint_t,2> bounds(T x1, U x2, const vec<Dim,Type>& v) {
    std::array<uint_t,2> res;

    if (v.empty()) {
        res[0] = npos;
        res[1] = npos;
    } else {
        auto iter = std::upper_bound(v.data.begin(), v.data.end(), x1,
            typename vec<Dim,Type>::comparator());

        if (iter == v.data.begin()) {
            res[0] = npos;
        } else {
            res[0] = iter - v.data.begin() - 1;
        }

        iter = std::upper_bound(iter, v.data.end(), x2,
            typename vec<1,Type>::comparator());

        if (iter == v.data.end()) {
            res[1] = npos;
        } else {
            res[1] = iter - v.data.begin();
        }
    }

    return res;
}

// Return the indices of all the values in the array that are equal to 'x'.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, std::size_t Dim, typename Type>
vec1u equal_range(T x, const vec<Dim,Type>& v) {
    auto res = std::equal_range(v.data.begin(), v.data.end(), x,
        typename vec<Dim,Type>::comparator());

    return uindgen(1 + (res.second - res.first)) + (res.first - v.data.begin());
}

template<std::size_t Dim, typename Type>
void inplace_remove(vec<Dim,Type>& v, vec1u ids) {
    inplace_sort(ids);
    uint_t i = 0;
    uint_t pitch = v.pitch(0);
    while (i < ids.size()) {
        uint_t i1 = ids.safe[ids.size()-1-i];
        uint_t i0 = i1;

        ++i;
        while (i < ids.size() && i0 - ids.safe[ids.size()-1-i] == 1) {
            i0 = ids.safe[ids.size()-1-i];
            ++i;
        }

        v.data.erase(v.data.begin()+i0*pitch, v.data.begin()+(i1+1)*pitch);
    }

    v.dims[0] -= ids.size();
}

template<std::size_t Dim, typename Type>
vec<Dim,Type> remove(vec<Dim,Type> v, const vec1u& ids) {
    inplace_remove(v, ids);
    return v;
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2 = Type1,
    typename enable = typename std::enable_if<(N < Dim)>::type>
void append(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
        return;
    }

    if (t2.empty()) return;

    std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
    phypp_check(t1.size()/n1 == t2.size()/n2, "cannot append dimension ", N, " in (", t1.dims,
        ") and (", t2.dims, ")");

    // TODO: optimize this copy
    auto tmp = t1;
    t1.dims[N] += n2;
    t1.resize();

    t1(repeat<N>(_), uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
    t1(repeat<N>(_), n1+uindgen(n2), repeat<Dim-N-1>(_)) = t2;
}

template<std::size_t N, std::size_t Dim, typename Type1, typename Type2 = Type1,
    typename enable = typename std::enable_if<(N < Dim)>::type>
void prepend(vec<Dim,Type1>& t1, const vec<Dim,Type2>& t2) {
    if (t1.empty()) {
        t1 = t2;
        return;
    }

    if (t2.empty()) return;

    std::size_t n1 = t1.dims[N], n2 = t2.dims[N];
    phypp_check(t1.size()/n1 == t2.size()/n2, "cannot prepend dimension ", N, " in (", t1.dims,
        ") and (", t2.dims, ")");

    // TODO: optimize this copy
    auto tmp = t1;
    t1.dims[N] += n2;
    t1.resize();

    t1(repeat<N>(_), uindgen(n2), repeat<Dim-N-1>(_)) = t2;
    t1(repeat<N>(_), n2+uindgen(n1), repeat<Dim-N-1>(_)) = tmp;
}

template<typename Type1, typename Type2 = Type1>
void append(vec<1,Type1>& t1, const vec<1,Type2>& t2) {
    t1.data.insert(t1.data.end(), t2.begin(), t2.end());
    t1.dims[0] += t2.dims[0];
}

template<typename Type1, typename Type2 = Type1>
void prepend(vec<1,Type1>& t1, const vec<1,Type2>& t2) {
    t1.data.insert(t1.data.begin(), t2.begin(), t2.end());
    t1.dims[0] += t2.dims[0];
}

template<typename Type>
vec<1,rtype_t<Type>> shift(const vec<1,Type>& v, int_t n, rtype_t<Type> def = 0) {
    vec<1,rtype_t<Type>> tmp;

    if (uint_t(std::abs(n)) > v.dims[0]) {
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

// Find the closest point in a 2D array that satisfies a given criterium
bool astar_find(const vec2b& map, uint_t& x, uint_t& y) {
    phypp_check(!map.empty(), "this algorithm requires a non empty 2D vector");

    if (x >= map.dims[0]) x = map.dims[0]-1;
    if (y >= map.dims[1]) y = map.dims[1]-1;

    if (map.safe(x,y)) return true;

    using vec_pair = vec<1,std::pair<uint_t,uint_t>>;
    vec_pair open;
    open.push_back(std::make_pair(x,y));

    vec2b visit(map.dims);
    visit.safe(x,y) = true;

    while (!open.empty()) {
        vec_pair old_open = std::move(open);

        for (auto p : old_open) {
            int_t ox = p.first, oy = p.second;

            for (uint_t d : range(4)) {
                int_t tnx, tny;
                if (d == 0) {
                    tnx = ox;   tny = oy+1;
                } else if (d == 1) {
                    tnx = ox+1; tny = oy;
                } else if (d == 2) {
                    tnx = ox;   tny = oy-1;
                } else {
                    tnx = ox-1; tny = oy;
                }

                if (tnx < 0 || tny < 0) continue;

                x = tnx, y = tny;
                if (x >= map.dims[0] || y >= map.dims[1] || visit.safe(x,y)) continue;

                if (!map.safe(x,y)) {
                    open.push_back(std::make_pair(x,y));
                    visit.safe(x,y) = true;
                } else {
                    return true;
                }
            }
        }
    }

    return false;
}

// Transform a scalar lambda into an overloaded lambda that supports both scalar and
// vector calls.
template<typename L>
struct vectorized_lambda_t {
    L lambda;

    vectorized_lambda_t(L tlam) : lambda(tlam) {}

    template<typename T, typename ... Args, typename enable =
        typename std::enable_if<!is_vec<typename std::decay<T>::type>::value>::type>
    auto operator()(T&& t, Args&& ... args) ->
        decltype(lambda(std::forward<T>(t), std::forward<Args>(args)...)) {
        return lambda(std::forward<T>(t), std::forward<Args>(args)...);
    }

    template<typename T, std::size_t D, typename ... Args>
    auto operator()(const vec<D,T>& t, Args&& ... args) ->
        vec<D,decltype(lambda(std::declval<const rtype_t<T>&>(), std::forward<Args>(args)...))> {
        vec<D,decltype(lambda(std::declval<const rtype_t<T>&>(), std::forward<Args>(args)...))> ret(t.dims);
        for (uint_t i : range(t)) {
            ret.safe[i] = lambda(t.safe[i], std::forward<Args>(args)...);
        }
        return ret;
    }

    template<typename T, std::size_t D, typename ... Args>
    auto operator()(vec<D,T>&& t, Args&& ... args) ->
        vec<D,decltype(lambda(std::declval<rtype_t<T>>(), std::forward<Args>(args)...))> {
        vec<D,decltype(lambda(std::declval<rtype_t<T>>(), std::forward<Args>(args)...))> ret(t.dims);
        for (uint_t i : range(t)) {
            ret.safe[i] = lambda(std::move(t.safe[i]), std::forward<Args>(args)...);
        }
        return ret;
    }
};

template<typename T>
vectorized_lambda_t<typename std::decay<T>::type> vectorize_lambda(T&& t) {
    vectorized_lambda_t<typename std::decay<T>::type> func(std::move(t));
    return func;
}

// Increment a list of indices, step by step
// This is an alternative to a recursive loop.
void increment_index_list(vec1u& ids, const vec1u& n) {
    uint_t i = ids.size();
    do {
        --i;
        ++(ids.safe[i]);
        if (ids.safe[i] == n.safe[i]) {
            ids.safe[i] = 0;
        } else {
            break;
        }
    } while (i != 0);
}

void increment_index_list(vec1u& ids, const uint_t& n) {
    uint_t i = ids.size();
    do {
        --i;
        ++(ids.safe[i]);
        if (ids.safe[i] == n) {
            ids.safe[i] = 0;
        } else {
            break;
        }
    } while (i != 0);
}

// Print a vector into a stream.
template<typename O, std::size_t Dim, typename Type, typename enable = typename std::enable_if<!std::is_same<Type, bool>::value>::type>
O& operator << (O& o, const vec<Dim,Type>& v) {
    o << '{';
    for (uint_t i = 0; i < v.data.size(); ++i) {
        if (i != 0) o << ", ";
        o << v.safe[i];
    }
    o << '}';

    return o;
}

template<typename O, std::size_t Dim>
O& operator << (O& o, const vec<Dim,bool>& v) {
    o << '{';
    for (uint_t i = 0; i < v.data.size(); ++i) {
        if (i != 0) o << ", ";
        o << bool(v.safe[i]);
    }
    o << '}';

    return o;
}

#endif

