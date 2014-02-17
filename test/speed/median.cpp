#include <phypp.hpp>

namespace speed_test {
    template<std::size_t D, typename T>
    rtype_t<T> median1(const vec_t<D,T>& v) {
        vec1u ok = where(!nan(v));
        if (ok.empty()) return dnan;

        typename vec_t<1,T>::effective_type t = v[ok];
        std::ptrdiff_t offset = t.size()/2;
        std::nth_element(t.begin(), t.begin() + offset, t.end());
        return *(t.begin() + offset);
    }

    template<std::size_t D, typename T>
    rtype_t<T> median2(const vec_t<D,T>& v) {
        vec1u ok = where(!nan(v));
        if (ok.empty()) return dnan;

        std::ptrdiff_t offset = ok.size()/2;
        std::nth_element(ok.begin(), ok.begin() + offset, ok.end(), [&v](uint_t i, uint_t j) {
            return v[i] < v[j];
        });
        return v[*(ok.begin() + offset)];
    }

    template<std::size_t D, typename T>
    rtype_t<T> median3(vec_t<D,T> v) {
        uint_t nwrong = 0;
        for (uint_t i : range(v)) {
            nwrong += nan(v[i]);
        }

        if (nwrong == v.size()) return dnan;

        std::ptrdiff_t offset = (v.size()-nwrong)/2;
        std::nth_element(v.begin(), v.begin() + offset, v.end(), [&v](rtype_t<T> i, rtype_t<T> j) {
            if (nan(i)) return false;
            if (nan(j)) return true;
            return i < j;
        });
        return *(v.begin() + offset);
    }

    template<std::size_t D, typename T>
    rtype_t<T> inplace_median3(vec_t<D,T>& v) {
        uint_t nwrong = 0;
        for (uint_t i : range(v)) {
            nwrong += nan(v[i]);
        }

        if (nwrong == v.size()) return dnan;

        std::ptrdiff_t offset = (v.size()-nwrong)/2;
        std::nth_element(v.begin(), v.begin() + offset, v.end(), [&v](rtype_t<T> i, rtype_t<T> j) {
            if (nan(i)) return false;
            if (nan(j)) return true;
            return i < j;
        });
        return *(v.begin() + offset);
    }
}

int main(int argc, char* argv[]) {
    uint_t nsrc = 1000000;
    uint_t navg = 1;
    bool check = false;

    read_args(argc, argv, arg_list(nsrc, navg, check));

    auto seed = make_seed(42);
    if (!check) {
        // Time
        vec1d data = randomn(seed, nsrc);
        vec1u id = randomi(seed, 0, nsrc-1, nsrc/10);
        data[id] = dnan;

        double res = 0.0;
        double t = profile([&data,&res]() {
            res += speed_test::median1(data);
        }, navg);

        print(t);

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median2(data);
        }, navg);

        print(t);

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median3(data);
        }, navg);

        print(t);

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::inplace_median3(data);
        }, navg);

        print(t);
    } else {
        // Check
        uint_t bad2 = 0;
        uint_t bad3 = 0;
        uint_t bad4 = 0;

        for (uint_t i = 0; i < navg; ++i) {
            vec1d data = randomn(seed, nsrc);
            vec1u id = randomi(seed, 0, nsrc-1, nsrc/10);
            data[id] = dnan;

            double ref = speed_test::median1(data);
            bad2 += speed_test::median2(data) != ref;
            bad3 += speed_test::median3(data) != ref;
            bad4 += speed_test::inplace_median3(data) != ref;
        }

        print(bad2);
        print(bad3);
        print(bad4);
    }

    return 0;
}
