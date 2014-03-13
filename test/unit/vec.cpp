#include <phypp.hpp>

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

template<typename T>
void assert_is_vec() {
    static_assert(is_vec<T>::value, "failed is_vec");
}

template<typename T>
void assert_not_is_vec() {
    static_assert(!is_vec<T>::value, "failed !is_vec");
}

template<typename T>
void test_is_vec() {
    // Metaprogramming
    vec1u r;

    assert_is_vec<vec_t<1,T>>();
    assert_is_vec<vec_t<2,T>>();
    assert_is_vec<vec_t<3,T>>();
    assert_is_vec<vec_t<4,T>>();

    assert_is_vec<decltype(vec_t<1,T>()(_))>();
    assert_is_vec<decltype(vec_t<1,T>()(r))>();
    assert_is_vec<decltype(vec_t<1,T>()[_])>();
    assert_is_vec<decltype(vec_t<1,T>()[r])>();
    assert_is_vec<decltype(vec_t<2,T>()[_])>();
    assert_is_vec<decltype(vec_t<2,T>()[r])>();

    assert_is_vec<decltype(vec_t<2,T>()(_,0))>();
    assert_is_vec<decltype(vec_t<2,T>()(r,0))>();
    assert_is_vec<decltype(vec_t<2,T>()(0,_))>();
    assert_is_vec<decltype(vec_t<2,T>()(0,r))>();

    assert_is_vec<decltype(vec_t<3,T>()(_,_,0))>();
    assert_is_vec<decltype(vec_t<3,T>()(r,_,0))>();
    assert_is_vec<decltype(vec_t<3,T>()(_,r,0))>();
    assert_is_vec<decltype(vec_t<3,T>()(r,r,0))>();
    assert_is_vec<decltype(vec_t<3,T>()(_,0,_))>();
    assert_is_vec<decltype(vec_t<3,T>()(r,0,_))>();
    assert_is_vec<decltype(vec_t<3,T>()(_,0,r))>();
    assert_is_vec<decltype(vec_t<3,T>()(r,0,r))>();
    assert_is_vec<decltype(vec_t<3,T>()(0,_,_))>();
    assert_is_vec<decltype(vec_t<3,T>()(0,r,_))>();
    assert_is_vec<decltype(vec_t<3,T>()(0,_,r))>();
    assert_is_vec<decltype(vec_t<3,T>()(0,r,r))>();

    assert_not_is_vec<T>();
    assert_not_is_vec<std::vector<T>>();
    assert_not_is_vec<decltype(vec_t<1,T>()(0))>();
    assert_not_is_vec<decltype(vec_t<2,T>()[0])>();
    assert_not_is_vec<decltype(vec_t<2,T>()(0,0))>();
    assert_not_is_vec<decltype(vec_t<3,T>()(0,0,0))>();
}

template<typename T, std::size_t D>
void assert_vec_dim() {
    static_assert(vec_dim<T>::value == D, "failed vec_dim");
}

template<typename T>
void test_vec_dim() {
    // Metaprogramming
    vec1u r;

    assert_vec_dim<vec_t<1,T>,1>();
    assert_vec_dim<vec_t<2,T>,2>();
    assert_vec_dim<vec_t<3,T>,3>();

    assert_vec_dim<decltype(vec_t<1,T>()(_)),1>();
    assert_vec_dim<decltype(vec_t<1,T>()(r)),1>();
    assert_vec_dim<decltype(vec_t<1,T>()[_]),1>();
    assert_vec_dim<decltype(vec_t<1,T>()[r]),1>();
    assert_vec_dim<decltype(vec_t<2,T>()[_]),1>();
    assert_vec_dim<decltype(vec_t<2,T>()[r]),1>();
    assert_vec_dim<decltype(vec_t<3,T>()[_]),1>();
    assert_vec_dim<decltype(vec_t<3,T>()[r]),1>();
    assert_vec_dim<decltype(vec_t<2,T>()(_,0)),1>();
    assert_vec_dim<decltype(vec_t<2,T>()(r,0)),1>();
    assert_vec_dim<decltype(vec_t<2,T>()(0,_)),1>();
    assert_vec_dim<decltype(vec_t<2,T>()(0,r)),1>();

    assert_vec_dim<decltype(vec_t<3,T>()(_,_,0)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(r,_,0)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(_,r,0)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(r,r,0)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(_,0,_)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(r,0,_)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(_,0,r)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(r,0,r)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(0,_,_)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(0,r,_)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(0,_,r)),2>();
    assert_vec_dim<decltype(vec_t<3,T>()(0,r,r)),2>();

    assert_vec_dim<T,0>();
    assert_vec_dim<std::vector<T>,0>();
    assert_vec_dim<std::array<T,5>,0>();
}

template<typename T>
struct generator {
    generator() {}

    static const std::array<T,20> list;

    T operator[] (uint_t i) const {
        return list[i];
    }
};

template<>
const std::array<uint_t,20> generator<uint_t>::list = {{
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
    10, 11, 12, 13, 14, 15, 16, 17, 18, 19
}};

template<>
const std::array<int_t,20> generator<int_t>::list = {{
   -9, -8, -7, -6, -5, -4, -3, -2, -1,  0,
    1,  2,  3,  4,  5,  6,  7,  8,  9, 10
}};

template<>
const std::array<float,20> generator<float>::list = {{
    (float)cos(1.0f), (float)sin(2.0f), (float)exp(1.0f), (float)log(2.0f),
    fpi,   0.0f, -1.0f,  fnan, +finf, -finf,
    1.0f, -2.0f,  3.0f, -4.0f,  5.0f, -6.0f, 7.0f, -8.0f, 9.0f, 1e6f
}};

template<>
const std::array<double,20> generator<double>::list = {{
    cos(1.0), sin(2.0), exp(1.0), log(2.0), dpi,  0.0, -1.0,  dnan, +dinf, -dinf,
    1.0,      -2.0,     3.0,      -4.0,     5.0, -6.0,  7.0, -8.0,   9.0,   1e6
}};

template<>
const std::array<std::string,20> generator<std::string>::list = {{
    "h",   "e",   "l",   "l",     "o",     "w",     "o", "r",  "l",  "d",
    "foo", "bar", "bee", "!blop", "flop#", "m lop", "",  "\t", "\n", "10"
}};

template<>
const std::array<bool,20> generator<bool>::list = {{
    true,  false, true, false, false, false, true,  true,  false, true,
    false, false, true, true,  true,  true,  false, false, false, true
}};

template<typename T>
void test_vec_constructor() {
    print("test_vec_constructor...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    const T def = std::vector<T>(1)[0];
    const generator<T> gen;

    {
        // Default constructor
        vec_t<1,T> v01;
        vec_t<2,T> v02;
        vec_t<3,T> v03;

        check(v01.dims[0], 0u);
        check(v02.dims[0], 0u);
        check(v02.dims[1], 0u);
        check(v03.dims[0], 0u);
        check(v03.dims[1], 0u);
        check(v03.dims[2], 0u);

        check(v01.size(), 0u);
        check(v02.size(), 0u);
        check(v03.size(), 0u);

        check(v01.empty(), true);
        check(v02.empty(), true);
        check(v03.empty(), true);
    }

    {
        // Size constructor
        const uint_t d1  = 10u;
        const uint_t d21 = 10u;
        const uint_t d22 = 20u;
        const uint_t d31 = 5u;
        const uint_t d32 = 7u;
        const uint_t d33 = 15u;
        const uint_t dr  = 3u;

        vec1u r(dr);
        vec_t<1,T> v1(d1);
        vec_t<2,T> v2(d21,d22);
        vec_t<3,T> v3(d31,d32,d33);

        check(v1.dims[0], d1);
        check(v2.dims[0], d21);
        check(v2.dims[1], d22);
        check(v3.dims[0], d31);
        check(v3.dims[1], d32);
        check(v3.dims[2], d33);
        check(v1[_].dims[0], d1);
        check(v1[r].dims[0], dr);
        check(v2[_].dims[0], d21*d22);
        check(v2[r].dims[0], dr);
        check(v3[_].dims[0], d31*d32*d33);
        check(v3[r].dims[0], dr);
        check(v2(_,0).dims[0], d21);
        check(v2(r,0).dims[0], dr);
        check(v2(0,_).dims[0], d22);
        check(v2(0,r).dims[0], dr);
        check(v3(_,_,0).dims[0], d31);
        check(v3(r,_,0).dims[0], dr);
        check(v3(_,r,0).dims[0], d31);
        check(v3(r,r,0).dims[0], dr);
        check(v3(_,_,0).dims[1], d32);
        check(v3(r,_,0).dims[1], d32);
        check(v3(_,r,0).dims[1], dr);
        check(v3(r,r,0).dims[1], dr);
        check(v3(_,0,_).dims[0], d31);
        check(v3(r,0,_).dims[0], dr);
        check(v3(_,0,r).dims[0], d31);
        check(v3(r,0,r).dims[0], dr);
        check(v3(_,0,_).dims[1], d33);
        check(v3(r,0,_).dims[1], d33);
        check(v3(_,0,r).dims[1], dr);
        check(v3(r,0,r).dims[1], dr);
        check(v3(0,_,_).dims[0], d32);
        check(v3(0,r,_).dims[0], dr);
        check(v3(0,_,r).dims[0], d32);
        check(v3(0,r,r).dims[0], dr);
        check(v3(0,_,_).dims[1], d33);
        check(v3(0,r,_).dims[1], d33);
        check(v3(0,_,r).dims[1], dr);
        check(v3(0,r,r).dims[1], dr);

        check(v1.size(), d1);
        check(v2.size(), d21*d22);
        check(v3.size(), d31*d32*d33);
        check(v1[_].size(), d1);
        check(v1[r].size(), dr);
        check(v2[_].size(), d21*d22);
        check(v2[r].size(), dr);
        check(v3[_].size(), d31*d32*d33);
        check(v3[r].size(), dr);
        check(v2(_,0).size(), d21);
        check(v2(r,0).size(), dr);
        check(v2(0,_).size(), d22);
        check(v2(0,r).size(), dr);
        check(v3(_,_,0).size(), d31*d32);
        check(v3(r,_,0).size(), dr*d32);
        check(v3(_,r,0).size(), d31*dr);
        check(v3(r,r,0).size(), dr*dr);
        check(v3(_,0,_).size(), d31*d33);
        check(v3(r,0,_).size(), dr*d33);
        check(v3(_,0,r).size(), d31*dr);
        check(v3(r,0,r).size(), dr*dr);
        check(v3(0,_,_).size(), d32*d33);
        check(v3(0,r,_).size(), dr*d33);
        check(v3(0,_,r).size(), d32*dr);
        check(v3(0,r,r).size(), dr*dr);

        check(v1.empty(), false);
        check(v2.empty(), false);
        check(v3.empty(), false);
        check(v1[_].empty(), false);
        check(v1[r].empty(), false);
        check(v2[_].empty(), false);
        check(v2[r].empty(), false);
        check(v3[_].empty(), false);
        check(v3[r].empty(), false);
        check(v2(_,0).empty(), false);
        check(v2(r,0).empty(), false);
        check(v2(0,_).empty(), false);
        check(v2(0,r).empty(), false);
        check(v3(_,_,0).empty(), false);
        check(v3(r,_,0).empty(), false);
        check(v3(_,r,0).empty(), false);
        check(v3(r,r,0).empty(), false);
        check(v3(_,0,_).empty(), false);
        check(v3(r,0,_).empty(), false);
        check(v3(_,0,r).empty(), false);
        check(v3(r,0,r).empty(), false);
        check(v3(0,_,_).empty(), false);
        check(v3(0,r,_).empty(), false);
        check(v3(0,_,r).empty(), false);
        check(v3(0,r,r).empty(), false);

        bool res = true;
        for (uint_t& v : r) {
            res = res && v == 0u;
        }
        check(res, true);

        res = true;
        for (T& v : v1) {
            res = res && v == def;
        }
        check(res, true);

        res = true;
        for (T& v : v2) {
            res = res && v == def;
        }
        check(res, true);

        res = true;
        for (T& v : v3) {
            res = res && v == def;
        }
        check(res, true);
    }

    {
        // Initializer list constructor
        vec_t<1,T> v1({gen[0], gen[1], gen[2], gen[3], gen[4]});
        check(v1.dims[0], 5u);
        check(v1.size(),  5u);
        check(v1.empty(), false);
        for (uint_t i = 0; i < 5u; ++i) {
            check(v1[i], gen[i]);
        }

        vec_t<1,T> v1e({});
        check(v1e.dims[0], 0u);
        check(v1e.size(),  0u);
        check(v1e.empty(), true);

        vec_t<2,T> v2({
            {gen[0], gen[1]},
            {gen[2], gen[3]},
            {gen[4], gen[5]}
        });
        check(v2.dims[0], 3u);
        check(v2.dims[1], 2u);
        check(v2.size(),  6u);
        check(v2.empty(), false);
        for (uint_t i = 0; i < 6u; ++i) {
            check(v2[i], gen[i]);
        }
        for (uint_t i = 0; i < 3u; ++i)
        for (uint_t j = 0; j < 2u; ++j) {
            check(v2(i,j), gen[i*2+j]);
        }

        vec_t<2,T> v2e({});
        check(v2e.dims[0], 0u);
        check(v2e.dims[1], 0u);
        check(v2e.size(),  0u);
        check(v2e.empty(), true);

        vec_t<3,T> v3({{
            {gen[0],  gen[1]},
            {gen[2],  gen[3]},
            {gen[4],  gen[5]}
        }, {
            {gen[6],  gen[7]},
            {gen[8],  gen[9]},
            {gen[10], gen[11]}
        }});
        check(v3.dims[0], 2u);
        check(v3.dims[1], 3u);
        check(v3.dims[2], 2u);
        check(v3.size(),  12u);
        check(v3.empty(), false);
        for (uint_t i = 0; i < 12u; ++i) {
            check(v3[i], gen[i]);
        }
        for (uint_t i = 0; i < 2u; ++i)
        for (uint_t j = 0; j < 3u; ++j)
        for (uint_t k = 0; k < 2u; ++k) {
            check(v3(i,j,k), gen[i*3*2+j*2+k]);
        }

        vec_t<3,T> v3e({});
        check(v3e.dims[0], 0u);
        check(v3e.dims[1], 0u);
        check(v3e.dims[2], 0u);
        check(v3e.size(),  0u);
        check(v3e.empty(), true);
    }

    {
        // Move constructor
        const uint_t d1  = 10u;
        const uint_t d21 = 10u;
        const uint_t d22 = 20u;
        const uint_t d31 = 5u;
        const uint_t d32 = 7u;
        const uint_t d33 = 15u;

        vec_t<1,T> v1(d1);
        vec_t<2,T> v2(d21,d22);
        vec_t<3,T> v3(d31,d32,d33);

        vec_t<1,T> nv1(std::move(v1));
        vec_t<2,T> nv2(std::move(v2));
        vec_t<3,T> nv3(std::move(v3));

        check(nv1.dims[0], d1);
        check(v1.dims[0],  0u);
        check(nv2.dims[0], d21);
        check(nv2.dims[1], d22);
        check(v2.dims[0],  0u);
        check(v2.dims[1],  0u);
        check(nv3.dims[0], d31);
        check(nv3.dims[1], d32);
        check(nv3.dims[2], d33);
        check(v3.dims[0],  0u);
        check(v3.dims[1],  0u);
        check(v3.dims[2],  0u);

        check(nv1.size(), d1);
        check(v1.size(),  0u);
        check(nv2.size(), d21*d22);
        check(v2.size(),  0u);
        check(nv3.size(), d31*d32*d33);
        check(v3.size(),  0u);

        bool res = true;
        for (T& v : nv1) {
            res = res && v == def;
        }
        check(res, true);

        res = true;
        for (T& v : nv2) {
            res = res && v == def;
        }
        check(res, true);

        res = true;
        for (T& v : nv3) {
            res = res && v == def;
        }
        check(res, true);
    }

    {
        // Move constructor (with actual values)
        vec_t<1,T> w1({gen[0], gen[1], gen[2], gen[3]});
        vec_t<2,T> w2({
            {gen[0], gen[1], gen[2]},
            {gen[3], gen[4], gen[5]}
        });
        vec_t<3,T> w3({{
            {gen[0],  gen[1]},
            {gen[2],  gen[3]},
            {gen[4],  gen[5]}
        }, {
            {gen[6],  gen[7]},
            {gen[8],  gen[9]},
            {gen[10], gen[11]}
        }});

        vec_t<1,T> nw1(std::move(w1));
        vec_t<2,T> nw2(std::move(w2));
        vec_t<3,T> nw3(std::move(w3));

        for (uint_t i = 0; i < 4; ++i) {
            check(nw1[i], gen[i]);
        }
        for (uint_t i = 0; i < 6; ++i) {
            check(nw2[i], gen[i]);
        }
        for (uint_t i = 0; i < 12; ++i) {
            check(nw3[i], gen[i]);
        }
    }

    {
        // Copy constructor
        const uint_t d1 = 10u;
        const uint_t d21 = 10u;
        const uint_t d22 = 20u;
        const uint_t d31 = 5u;
        const uint_t d32 = 7u;
        const uint_t d33 = 15u;

        vec_t<1,T> v1(d1);
        vec_t<2,T> v2(d21,d22);
        vec_t<3,T> v3(d31,d32,d33);

        vec_t<1,T> nv1(v1);
        vec_t<2,T> nv2(v2);
        vec_t<3,T> nv3(v3);

        check(nv1.dims[0], d1);
        check(v1.dims[0],  d1);
        check(nv2.dims[0], d21);
        check(nv2.dims[1], d22);
        check(v2.dims[0],  d21);
        check(v2.dims[1],  d22);
        check(nv3.dims[0], d31);
        check(nv3.dims[1], d32);
        check(nv3.dims[2], d33);
        check(v3.dims[0],  d31);
        check(v3.dims[1],  d32);
        check(v3.dims[2],  d33);

        check(nv1.size(), d1);
        check(v1.size(),  d1);
        check(nv2.size(), d21*d22);
        check(v2.size(),  d21*d22);
        check(nv3.size(), d31*d32*d33);
        check(v3.size(),  d31*d32*d33);

        bool res = true;
        for (T& v : nv1) {
            res = res && v == def;
        }
        check(res, true);
        res = true;
        for (T& v : nv2) {
            res = res && v == def;
        }
        check(res, true);
        res = true;
        for (T& v : nv3) {
            res = res && v == def;
        }
        check(res, true);
    }

    {
        // Copy constructor (with actual values)
        vec_t<1,T> w1({gen[0], gen[1], gen[2], gen[3]});
        vec_t<2,T> w2({
            {gen[0], gen[1], gen[2]},
            {gen[3], gen[4], gen[5]}
        });
        vec_t<3,T> w3({{
            {gen[0],  gen[1]},
            {gen[2],  gen[3]},
            {gen[4],  gen[5]}
        }, {
            {gen[6],  gen[7]},
            {gen[8],  gen[9]},
            {gen[10], gen[11]}
        }});

        vec_t<1,T> nw1(w1);
        vec_t<2,T> nw2(w2);
        vec_t<3,T> nw3(w3);

        for (uint_t i : range(w1)) {
            check(nw1[i], w1[i]);
        }
        for (uint_t i : range(w2)) {
            check(nw2[i], w2[i]);
        }
        for (uint_t i : range(w3)) {
            check(nw3[i], w3[i]);
        }
    }

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T>
void test_vec_index() {
    print("test_vec_index...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    const generator<T> gen;

    vec_t<2,T> v = {
        {gen[0], gen[1], gen[2], gen[3]},
        {gen[4], gen[5], gen[6], gen[7]}
    };

    {
        vec_t<2,T> v1 = v;

        // Pure integer index
        for (uint_t i = 0; i < 8; ++i) {
            check(v1[i], gen[i]);
        }
        for (uint_t i = 0; i < 2; ++i)
        for (uint_t j = 0; j < 4; ++j) {
            check(v1(i,j), gen[i*4+j]);
        }

        v1[0] = gen[8];
        check(v1[0],   gen[8]);
        check(v1(0,0), gen[8]);

        v1(0,0) = gen[9];
        check(v1[0],   gen[9]);
        check(v1(0,0), gen[9]);

        v1(1,2) = gen[10];
        check(v1[6],   gen[10]);
        check(v1(1,2), gen[10]);

        v1(1,3) = gen[11];
        check(v1[7],   gen[11]);
        check(v1(1,3), gen[11]);

        vec_t<2,T> tv = {
            {gen[9], gen[1], gen[2],  gen[3]},
            {gen[4], gen[5], gen[10], gen[11]}
        };
        check(v1, tv);
    }

    {
        vec_t<2,T> v1 = v;

        // Placeholder index
        vec_t<1,T> sv0 = {gen[0], gen[1], gen[2], gen[3], gen[4], gen[5], gen[6], gen[7]};
        vec_t<1,T> sv1 = {gen[0], gen[1], gen[2], gen[3]};
        vec_t<1,T> sv2 = {gen[4], gen[5], gen[6], gen[7]};
        vec_t<1,T> sv3 = {gen[0], gen[4]};
        vec_t<1,T> sv4 = {gen[1], gen[5]};
        vec_t<1,T> sv5 = {gen[2], gen[6]};
        vec_t<1,T> sv6 = {gen[3], gen[7]};

        check(v1[_],   sv0);
        check(v1(0,_), sv1);
        check(v1(1,_), sv2);
        check(v1(_,0), sv3);
        check(v1(_,1), sv4);
        check(v1(_,2), sv5);
        check(v1(_,3), sv6);
        check(v1(_,_), v1);

        v1(0,_) = sv2;
        check(v1(0,_), sv2);

        v1(1,_) = sv1;
        check(v1(1,_), sv1);

        vec_t<2,T> tv = {
            {gen[4], gen[5], gen[6], gen[7]},
            {gen[0], gen[1], gen[2], gen[3]}
        };
        check(v1, tv);

        v1(0,_) = v1(1,_);
        tv = {
            {gen[0], gen[1], gen[2], gen[3]},
            {gen[0], gen[1], gen[2], gen[3]}
        };
        check(v1, tv);

        v1(1,_) = v(1,_);
        check(v1, v);

        v1(_,0) = sv4;
        check(v1(_,0), sv4);

        v1(_,1) = sv5;
        check(v1(_,1), sv5);

        v1(_,2) = sv6;
        check(v1(_,2), sv6);

        v1(_,3) = sv3;
        check(v1(_,3), sv3);

        tv = {
            {gen[1], gen[2], gen[3], gen[0]},
            {gen[5], gen[6], gen[7], gen[4]}
        };
        check(v1, tv);
    }

    {
        vec_t<2,T> v1 = v;

        // Vector index
        vec1u r1 = {0,1,2,3};
        vec1u r2 = {1,0};

        vec_t<1,T> sv0 = {gen[1], gen[0]};
        vec_t<1,T> sv1 = {gen[0], gen[1], gen[2], gen[3]};
        vec_t<1,T> sv2 = {gen[5], gen[4]};
        vec_t<1,T> sv3 = {gen[4], gen[5], gen[6], gen[7]};
        vec_t<1,T> sv4 = {gen[4], gen[0]};
        vec_t<1,T> sv5 = {gen[5], gen[1]};
        vec_t<1,T> sv6 = {gen[6], gen[2]};
        vec_t<1,T> sv7 = {gen[7], gen[3]};
        vec_t<2,T> sv8 = {
            {gen[5], gen[4]},
            {gen[1], gen[0]}
        };
        vec_t<2,T> sv9 = {
            {gen[4], gen[5], gen[6], gen[7]},
            {gen[0], gen[1], gen[2], gen[3]}
        };

        check(v1(0,r2),  sv0);
        check(v1(0,r1),  sv1);
        check(v1(1,r2),  sv2);
        check(v1(1,r1),  sv3);
        check(v1(r2,0),  sv4);
        check(v1(r2,1),  sv5);
        check(v1(r2,2),  sv6);
        check(v1(r2,3),  sv7);
        check(v1(r2,r2), sv8);
        check(v1(r2,r1), sv9);

        v1(0,r1) = sv3;
        check(v1(0,r1), sv3);

        v1(1,r1) = sv1;
        check(v1(1,r1), sv1);
        check(v1, sv9);

        v1(r2,0) = sv4;
        check(v1(r2,0), sv4);

        v1(r2,1) = sv5;
        check(v1(r2,1), sv5);

        v1(r2,2) = sv6;
        check(v1(r2,2), sv6);

        v1(r2,3) = sv7;
        check(v1(r2,3), sv7);

        check(v1, v);

        v1(vec1u{1,0},_) = v1(vec1u{0,1},_);
        check(v1, sv9);
    }

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T>
void test_vec_assign() {
    print("test_vec_assign...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    const generator<T> gen;

    // TODO

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T, typename U>
void test_vec_convert_to() {
    // TODO
}

template<typename T>
struct convert_to {
    using list = type_list<>;
};

template<>
struct convert_to<uint_t> {
    using list = type_list<int_t,float,double>;
};

template<>
struct convert_to<int_t> {
    using list = type_list<uint_t,float,double>;
};

template<>
struct convert_to<float> {
    using list = type_list<uint_t,int_t,double>;
};

template<>
struct convert_to<double> {
    using list = type_list<uint_t,int_t,float>;
};

template<typename T>
void test_vec_convert_iter(type_list<> list) {}

template<typename T, typename U, typename ... Args>
void test_vec_convert_iter(type_list<U,Args...> list) {
    test_vec_convert_to<T,U>();
    test_vec_convert_iter<T>(pop_front(list));
}

template<typename T>
void test_vec_convert() {
    print("test_vec_convert...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    test_vec_convert_iter<T>(typename convert_to<T>::list{});

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T>
void test_vec_iterator() {
    print("test_vec_iterator...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    const generator<T> gen;

    // TODO

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T>
void test_vec_operator() {
    print("test_vec_operator...");

    uint_t old_tested = tested;
    uint_t old_failed = failed;

    const generator<T> gen;

    // TODO

    print("> ", tested - failed - (old_tested - old_failed), "/", tested - old_tested," passed");
}

template<typename T>
void test() {
    test_is_vec<T>();
    test_vec_dim<T>();
    test_vec_constructor<T>();
    test_vec_index<T>();
    test_vec_assign<T>();
    test_vec_convert<T>();
    test_vec_iterator<T>();
    test_vec_operator<T>();
}

int main(int argc, char* argv[]) {
    test<uint_t>();
    test<float>();
    test<double>();
    test<std::string>();
    test<bool>();

    print("total:");
    print("> ", tested - failed, "/", tested," passed");

    return failed == 0u ? 0 : 1;
}
