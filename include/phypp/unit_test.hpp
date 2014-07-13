#ifndef UNIT_TEST_HPP
#define UNIT_TEST_HPP

uint_t tested = 0u;
uint_t failed = 0u;

bool reduce_and(bool b) {
    return b;
}

template<std::size_t Dim>
bool reduce_and(const vec_t<Dim,bool>& b) {
    bool res = true;
    for (bool v : b) {
        res = res && v;
    }
    return res;
}

template<typename T1>
struct is_float_ : std::false_type {
};
template<>
struct is_float_<float> : std::true_type {
};
template<>
struct is_float_<double> : std::true_type {
};
template<std::size_t D, typename T>
struct is_float_<vec_t<D,T>> : is_float_<typename std::remove_pointer<T>::type> {
};

template<typename T1, typename T2>
struct is_float : std::integral_constant<bool, is_float_<T1>::value || is_float_<T2>::value> {
};

template<typename T1, typename T2>
auto is_same_(const T1& v1, const T2& v2, std::false_type) -> decltype(v1 == v2) {
    return v1 == v2;
}

template<typename T1, typename T2>
auto is_same_(const T1& v1, const T2& v2, std::true_type) -> decltype(v1 == v2) {
    return fabs(v1 - v2) <= std::numeric_limits<float>::epsilon() ||
        (nan(v1) && nan(v2)) || nan(v1 - v2);
}

template<typename T1, typename T2>
auto is_same(const T1& v1, const T2& v2) -> decltype(v1 == v2) {
    return is_same_(v1, v2, is_float<T1,T2>{});
}

#define check_base_(cond, msg, line) do { \
    print("vec.cpp:"+strn(line)); \
    if (!cond) { \
        print(msg); \
        ++failed; \
    } \
    ++tested; \
} while (false)

#define check_base(cond, msg) check_base_(reduce_and(cond), msg, __LINE__)

#define check(t, s) \
    check_base(is_same(t, s), "  failed: "+std::string(#t)+" = "+strn(t)+" != "+strn(s))

#endif
