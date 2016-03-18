template<typename T>
std::string strn(const T&);

namespace impl {
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

    inline std::string pretty_type_(type_list<char>) {
        return "char";
    };

    inline std::string pretty_type_(type_list<short>) {
        return "short";
    };

    inline std::string pretty_type_(type_list<int>) {
        return "int";
    };

    inline std::string pretty_type_(type_list<long>) {
        return "long";
    };

    inline std::string pretty_type_(type_list<long long>) {
        return "llong";
    };

    inline std::string pretty_type_(type_list<unsigned char>) {
        return "uchar";
    };

    inline std::string pretty_type_(type_list<unsigned short>) {
        return "ushort";
    };

    inline std::string pretty_type_(type_list<unsigned int>) {
        return "uint";
    };

    inline std::string pretty_type_(type_list<unsigned long>) {
        return "ulong";
    };

    inline std::string pretty_type_(type_list<unsigned long long>) {
        return "ullong";
    };

    inline std::string pretty_type_(type_list<bool>) {
        return "bool";
    };

    inline std::string pretty_type_(type_list<float>) {
        return "float";
    };

    inline std::string pretty_type_(type_list<double>) {
        return "double";
    };

    inline std::string pretty_type_(type_list<std::string>) {
        return "string";
    };

    template<std::size_t Dim, typename T>
    std::string pretty_type_(type_list<vec<Dim,T>>) {
        return "vec<"+strn(Dim)+","+pretty_type_(type_list<T>{})+">";
    };
}

#define pretty_type(x) ::impl::pretty_type_(type_list<decltype(x)>{})
#define pretty_type_t(x) ::impl::pretty_type_(type_list<x>{})

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

// Helper to get the decayed data type of a vector
// vec<D,T>       = T
// vec<D,T*>      = T
// vec<D,const T> = T
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

// Compute the total number of dimensions generated by a set of arguments (for indexing)
template<typename T>
struct dim_index : std::integral_constant<std::size_t, 1> {};
template<std::size_t N, typename T>
struct dim_index<std::array<T,N>> : std::integral_constant<std::size_t, N> {};

template<typename T, typename ... Args>
struct dim_total : std::integral_constant<std::size_t,
    dim_index<typename std::decay<T>::type>::value + dim_total<Args...>::value> {};

template<typename T>
struct dim_total<T> : std::integral_constant<std::size_t,
    dim_index<typename std::decay<T>::type>::value> {};

// Assign an arbitrary list of values to an array.
// Supports scalars and other arrays.
template<std::size_t N, typename T, std::size_t I>
void set_array_(std::array<T,N>& v, cte_t<I>) {}

template<std::size_t N, typename T, std::size_t I, typename U, typename ... Args>
void set_array_(std::array<T,N>& v, cte_t<I>, const U& t, Args&& ... args) {
    v[I] = static_cast<T>(t);
    set_array_(v, cte_t<I+1>{}, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, std::size_t I, typename U, std::size_t M, typename ... Args>
void set_array_(std::array<T,N>& v, cte_t<I>, const std::array<U,M>& t, Args&& ... args) {
    for (uint_t i : range(M)) {
        v[I+i] = static_cast<T>(t[i]);
    }

    set_array_(v, cte_t<I+M>{}, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, typename ... Args>
void set_array(std::array<T,N>& v, Args&& ... args) {
    static_assert(N == dim_total<Args...>::value, "wrong number of elements for this array");
    set_array_(v, cte_t<0>{}, std::forward<Args>(args)...);
}

// Trait to figure out if type list matches an dimension list
template<typename T>
struct is_dim_elem : std::is_arithmetic<T> {};

template<typename T, std::size_t N>
struct is_dim_elem<std::array<T,N>> : std::true_type {};

template<typename ... Args>
struct is_dim_list : std::integral_constant<bool,
    are_all_true<bool_list<is_dim_elem<typename std::decay<Args>::type>::value...>>::value> {};

template<>
struct is_dim_list<> : std::true_type {};

// Trait to define if a vector type can be converted into another
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
