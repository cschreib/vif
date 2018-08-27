#include <vif.hpp>

using namespace vif;

namespace speed_test {
    template<std::size_t D, typename T>
    rtype_t<T> median1(const vec<D,T>& v) {
        vec1u ok = where(!is_nan(v));
        if (ok.empty()) return dnan;

        typename vec<1,T>::effective_type t = v[ok];
        std::ptrdiff_t offset = t.size()/2;
        std::nth_element(t.begin(), t.begin() + offset, t.end());
        return *(t.begin() + offset);
    }

    template<std::size_t D, typename T>
    rtype_t<T> median2(const vec<D,T>& v) {
        // vec1u ok = where(!is_nan(v));
        // if (ok.empty()) return dnan;

        vec1u ok; ok.data.reserve(v.size());
        for (uint_t i : range(v)) {
            if (!is_nan(v.safe[i])) ok.data.push_back(i);
        }

        ok.dims[0] = ok.data.size();

        std::ptrdiff_t offset = ok.size()/2;
        std::nth_element(ok.begin(), ok.begin() + offset, ok.end(),
            [&v](uint_t i, uint_t j) {
            return v.safe[i] < v.safe[j];
        });
        return v.safe[*(ok.begin() + offset)];
    }

    template<std::size_t D, typename T>
    rtype_t<T> median3(vec<D,T> v) {
        uint_t nwrong = 0;
        for (auto d : v) {
            nwrong += is_nan(d);
        }

        if (nwrong == v.size()) return dnan;

        std::ptrdiff_t offset = (v.size()-nwrong)/2;
        std::nth_element(v.begin(), v.begin() + offset, v.end(),
            [&v](rtype_t<T> i, rtype_t<T> j) {
            if (is_nan(i)) return false;
            if (is_nan(j)) return true;
            return i < j;
        });
        return *(v.begin() + offset);
    }

    template<std::size_t D, typename T>
    rtype_t<T> inplace_median3(vec<D,T>& v) {
        uint_t nwrong = 0;
        for (auto d : v) {
            nwrong += is_nan(d);
        }

        if (nwrong == v.size()) return dnan;

        std::ptrdiff_t offset = (v.size()-nwrong)/2;
        std::nth_element(v.begin(), v.begin() + offset, v.end(),
            [&v](rtype_t<T> i, rtype_t<T> j) {
            if (is_nan(i)) return false;
            if (is_nan(j)) return true;
            return i < j;
        });
        return *(v.begin() + offset);
    }
}

int vif_main(int argc, char* argv[]) {
    uint_t nsrc = 1000000;
    uint_t navg = 1;
    bool check = false;

    read_args(argc, argv, arg_list(nsrc, navg, check));

    auto seed = make_seed(42);
    if (!check) {
        // Time
        vec1d odata = randomn(seed, nsrc);
        vec1d data = odata;
        vec1u id = randomi(seed, 0, nsrc-1, nsrc/10);
        data[id] = dnan;

        double res = 0.0;
        double t = profile([&data,&res]() {
            res += speed_test::median1(data);
        }, navg);

        print(t);
        data = odata;

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median2(data);
        }, navg);

        print(t);
        data = odata;

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median3(data);
        }, navg);

        print(t);
        data = odata;

        res = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::inplace_median3(data);
        }, navg);

        print(t);
        data = odata;

        res = 0.0;
        t = profile([&data,&res]() {
            res += median(data);
        }, navg);

        print(t);
        data = odata;

        res = 0.0;
        t = profile([&data,&res]() {
            res += inplace_median(data);
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
