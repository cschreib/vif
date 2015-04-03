#include <phypp.hpp>

template<typename T>
void prof(const std::string& msg, T&& t) {
    double d = profile(std::forward<T>(t), 1);
    print(msg, ": ", seconds_str(d));
}

template<typename T>
struct range_t;

template<typename T>
struct range_iterator_t;

template<typename T>
struct range_iterator_t<range_t<T>> {
    const range_t<T>& range;
    std::size_t i;

    range_iterator_t& operator++ (int) {
        ++i;
        return *this;
    }

    range_iterator_t operator++ () {
        return range_iterator_t{range, i++};
    }

    T operator * () const {
        return range.b + range.d*i;
    }

    bool operator == (const range_iterator_t& iter) const {
        return iter.i == i && iter.range == range;
    }

    bool operator != (const range_iterator_t& iter) const {
        return iter.i != i || &iter.range != &range;
    }
};

template<typename T>
struct range_t {
    T b, e, d;
    std::size_t n;

    range_t(T b_, T e_, std::size_t n_) : b(b_), e(e_), d((e_-b_)/n_), n(n_) {}

    using iterator = range_iterator_t<range_t>;

    iterator begin() const { return iterator{*this, 0}; }
    iterator end() const { return iterator{*this, n}; }
};

template<typename T>
range_t<T> range(T i, T e, std::size_t n) {
    return range_t<T>(i, e, n);
}

template<typename T>
range_t<T> range(T i, T e) {
    return range(i, e, e-i);
}

template<typename T>
range_t<T> range(T n) {
    return range(T(0), n);
}

template<typename T, typename enable = typename std::enable_if<is_vec<T>::value>::type>
range_t<T> range(T n) {
    return range(0u, n.size()-1);
}

int main(int argc, char* argv[]) {
    uint_t test = 0;

    read_args(argc, argv, arg_list(test));

    if (test == 0) {
        uint_t n = 1e7;
        vec1i u = indgen(n);

        int* v0;
        vec1i v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11;
        std::vector<int> v12, v12b, v13;
        prof("C ref", [&]() {
            v0 = (int*)malloc(n*sizeof(int));
            for (uint_t i = 0; i < n; ++i) {
                v0[i] = 2*u.data[i];
            }
        });
        prof("std::vector ref", [&]() {
            v12.resize(n);
            for (uint_t i = 0; i < n; ++i) {
                v12[i] = 2*u.data[i];
            }
        });
        prof("std::vector ref2", [&]() {
            v12b = std::vector<int>(n);
            for (uint_t i = 0; i < n; ++i) {
                v12b[i] = 2*u.data[i];
            }
        });
        prof("Simple loop .data", [&]() {
            v1 = intarr(n);
            for (uint_t i = 0; i < n; ++i) {
                v1.data[i] = 2*u.data[i];
            }
        });
        prof("Simple loop", [&]() {
            v2 = intarr(n);
            for (uint_t i = 0; i < n; ++i) {
                v2[i] = 2*u[i];
            }
        });
        prof("Vectorized loop", [&]() {
            v3 = 2*u;
        });
        prof("Copy and modify simple .data", [&]() {
            v4 = u;
            for (uint_t i = 0; i < n; ++i) {
                v4.data[i] *= 2;
            }
        });
        prof("Copy and modify simple", [&]() {
            v5 = u;
            for (uint_t i = 0; i < n; ++i) {
                v5[i] *= 2;
            }
        });
        prof("Copy and modify auto", [&]() {
            v6 = u;
            for (auto& t : v6) {
                t *= 2;
            }
        });
        prof("Copy and modify vectorized", [&]() {
            v7 = u;
            v7 *= 2;
        });
        prof("Push_back std::vector ref", [&]() {
            v13.reserve(n);
            for (uint_t i = 0; i < n; ++i) {
                v13.push_back(2*u.data[i]);
            }
        });
        prof("Push_back std::vector auto ref", [&]() {
            v13.reserve(n);
            for (auto& t : u) {
                v13.push_back(2*t);
            }
        });
        prof("Push_back", [&]() {
            v8.reserve(n);
            for (uint_t i = 0; i < n; ++i) {
                v8.push_back(2*u[i]);
            }
        });
        prof("Push_back auto", [&]() {
            v9.reserve(n);
            for (auto& t : u) {
                v9.push_back(2*t);
            }
        });
        prof("Push_back .data", [&]() {
            v10.reserve(n);
            for (uint_t i = 0; i < n; ++i) {
                v10.data.push_back(2*u.data[i]);
            }
            v10.dims[0] = v10.data.size();
        });
        prof("Push_back auto .data", [&]() {
            v11.reserve(n);
            for (auto& t : u) {
                v11.data.push_back(2*t);
            }
            v11.dims[0] = v11.data.size();
        });

        print(vec1i{v0[0], v0[1], v0[2], v0[3]});
        print(vec1i{v1[0], v1[1], v1[2], v1[3]});
        print(vec1i{v2[0], v2[1], v2[2], v2[3]});
        print(vec1i{v3[0], v3[1], v3[2], v3[3]});
        print(vec1i{v4[0], v4[1], v4[2], v4[3]});
        print(vec1i{v5[0], v5[1], v5[2], v5[3]});
        print(vec1i{v6[0], v6[1], v6[2], v6[3]});
        print(vec1i{v7[0], v7[1], v7[2], v7[3]});
        print(vec1i{v8[0], v8[1], v8[2], v8[3]});
        print(vec1i{v9[0], v9[1], v9[2], v9[3]});
        print(vec1i{v10[0], v10[1], v10[2], v10[3]});
        print(vec1i{v11[0], v11[1], v11[2], v11[3]});
        print(vec1i{v12[0], v12[1], v12[2], v12[3]});
        print(vec1i{v12b[0], v12b[1], v12b[2], v12b[3]});
        print(vec1i{v13[0], v13[1], v13[2], v13[3]});
    } else if (test == 1) {
        double d1 = 0.1, d2 = 0.1, d3 = 0.1;
        prof("rgen", [&]() {
            for (auto u : rgen(1e6)) {
                d1 += u*u + fabs(log10(d1));
            }
        });
        prof("range", [&]() {
            for (auto u : range(uint_t(1e6))) {
                d2 += u*u + fabs(log10(d2));
            }
        });
        prof("C", [&]() {
            for (uint_t u = 0; u < 1e6; ++u) {
                d3 += u*u + fabs(log10(d3));
            }
        });

        print(d1);
        print(d2);
        print(d3);
    }

    return 0;
}
