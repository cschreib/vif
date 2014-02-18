#include <phypp.hpp>

namespace speed_test {
    template<std::size_t Dim, typename Type>
    vec_t<Dim-1,rtype_t<Type>> median1(const vec_t<Dim,Type>& v, uint_t dim) {
        using fptr = rtype_t<Type> (*)(vec_t<1,rtype_t<Type>>);
        return run_index_<fptr, &median<1,rtype_t<Type>>>(v, dim);
    }

    template<typename I>
    struct stride_iterator {
        I iter_;
        uint_t stride_;

        stride_iterator(I i, uint_t stride) : iter_(i), stride_(stride) {}

        stride_iterator& operator ++ (int) {
            iter_ += stride_;
            return *this;
        }
        stride_iterator operator ++ () {
            iter_ += stride_;
            return *this;
        }
        stride_iterator& operator -- (int) {
            iter_ -= stride_;
            return *this;
        }
        stride_iterator operator -- () {
            iter_ -= stride_;
            return *this;
        }
        void operator += (int_t n) {
            iter_ += n*int_t(stride_);
        }
        void operator -= (int_t n) {
            iter_ -= n*int_t(stride_);
        }
        stride_iterator operator + (int_t n) const {
            stride_iterator i = *this;
            i += n;
            return i;
        }
        stride_iterator operator - (int_t n) const {
            stride_iterator i = *this;
            i -= n;
            return i;
        }
        std::ptrdiff_t operator - (const stride_iterator& iter) const {
            return (iter_ - iter.iter_)/stride_;
        }
        bool operator == (const stride_iterator& iter) const {
            return iter_ == iter.iter_;
        }
        bool operator != (const stride_iterator& iter) const {
            return iter_ != iter.iter_;
        }
        bool operator < (const stride_iterator& iter) const {
            return iter_ < iter.iter_;
        }
        bool operator <= (const stride_iterator& iter) const {
            return iter_ <= iter.iter_;
        }
        bool operator > (const stride_iterator& iter) const {
            return iter_ > iter.iter_;
        }
        bool operator >= (const stride_iterator& iter) const {
            return iter_ >= iter.iter_;
        }
        auto operator * () -> decltype(*iter_) {
            return *iter_;
        }
    };

    template<typename I>
    std::pair<stride_iterator<I>,stride_iterator<I>> stride(I i0, I i1, uint_t stride) {
        uint_t mod = (i1 - i0) % stride;
        return std::make_pair(
            stride_iterator<I>(i0, stride),
            stride_iterator<I>(i1 + (mod != 0 ? stride - mod : 0), stride)
        );
    }

    template<typename I>
    std::pair<stride_iterator<I>,stride_iterator<I>> stride(I i0, uint_t stride, uint_t niter) {
        return std::make_pair(
            stride_iterator<I>(i0, stride),
            stride_iterator<I>(i0 + niter*stride, stride)
        );
    }
}

namespace std {
    template<typename I>
    struct iterator_traits<speed_test::stride_iterator<I>> : iterator_traits<I> {};
}

namespace speed_test {
    template<std::size_t Dim, typename Type>
    vec_t<Dim-1,rtype_t<Type>> median2(vec_t<Dim,Type> v, uint_t dim) {
        vec_t<Dim-1,rtype_t<Type>> r;
        for (uint_t i = 0; i < dim; ++i) {
            r.dims[i] = v.dims[i];
        }
        for (uint_t i = dim+1; i < Dim; ++i) {
            r.dims[i-1] = v.dims[i];
        }
        r.resize();

        uint_t mpitch = 1;
        for (uint_t i = dim+1; i < Dim; ++i) {
            mpitch *= v.dims[i];
        }

        for (uint_t i = 0; i < r.size(); ++i) {
            uint_t base = (i%mpitch) + (i/mpitch)*v.dims[dim]*mpitch;

            auto p = stride(v.begin() + base, mpitch, v.dims[dim]);

            uint_t nwrong = 0;
            for (auto iter = p.first; iter != p.second; ++iter) {
                nwrong += nan(*iter);
            }

            if (nwrong == v.dims[dim]) {
                r[i] = dnan;
            } else {
                std::ptrdiff_t offset = (v.dims[dim]-nwrong)/2;
                std::nth_element(p.first, p.first + offset, p.second,
                    [](rtype_t<Type> x0, rtype_t<Type> x1) {
                        if (nan(x0)) return false;
                        if (nan(x1)) return true;
                        return x0 < x1;
                    }
                );

                r[i] = *(p.first + offset);
            }
        }

        return r;
    }
}

int main(int argc, char* argv[]) {
    uint_t nsrc = 10000;
    uint_t navg = 1;
    bool check = false;

    read_args(argc, argv, arg_list(nsrc, navg, check));

    auto seed = make_seed(42);

    if (check) {
        uint_t bad2 = 0;

        for (uint_t i : range(navg)) {
            vec3d data = randomn(seed, nsrc, 61, 61);
            vec2d ref = speed_test::median1(data, 0);
            bad2 += total(speed_test::median2(data, 0) != ref);
        }

        print(bad2);
    } else {
        vec3d data = randomn(seed, nsrc, 61, 61);
        vec1u id = randomi(seed, 0, nsrc*61*61-1, nsrc/10);
        data[id] = dnan;

        vec2d res(61,61);
        double t;

        res[_] = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median1(data, 0);
        }, navg);

        print(t);

        res[_] = 0.0;
        t = profile([&data,&res]() {
            res += speed_test::median2(data, 0);
        }, navg);

        print(t);
    }

    return 0;
}
