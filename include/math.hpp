#ifndef MATH_HPP
#define MATH_HPP

#include "vec.hpp"
#include "string.hpp"
#include <random>
#include <cmath>

#ifndef NO_GSL
#include <gsl/gsl_sf_bessel.h>
#endif

static const double dnan = std::numeric_limits<double>::quiet_NaN();
static const float  fnan = std::numeric_limits<float>::quiet_NaN();
static const double dinf = std::numeric_limits<double>::infinity();
static const float  finf = std::numeric_limits<float>::infinity();
static const double dpi = 3.14159265359;
static const float  fpi = 3.14159265359;

template<typename T>
auto e10(const T& t) {
    return pow(10.0, t);
}

// Create a range.
template<typename T>
vec1u rgen(T n) {
    phypp_check(n >= 0, "'rgen(n)' needs a positive or null value for 'n' (got ", n, ")");

    vec1u v = uintarr(n);
    for (uint_t k = 0; k < uint_t(n); ++k) {
        v[k] = k;
    }

    return v;
}

template<typename T, typename U>
vec_t<1,T> rgen(T i, U j) {
    if (i < T(j)) {
        uint_t n = j-i+1;
        vec_t<1,T> v = arr<T>(n);
        for (uint_t k = 0; k < n; ++k) {
            v[k] = i+k;
        }
        return v;
    } else {
        uint_t n = i-j+1;
        vec_t<1,T> v = arr<T>(n);
        for (uint_t k = 0; k < n; ++k) {
            v[k] = i-k;
        }
        return v;
    }
}

template<typename T, typename U, typename V>
vec1d rgen(T i, U j, V n) {
    phypp_check(n >= 0, "'rgen(a,b,n)' needs a positive or null value for 'n' (got ", n, ")");

    if (n == 1) {
        vec1d v = dindgen(1);
        v[0] = i;
        return v;
    } else {
        vec1d v = dindgen(n);

        double dx = (j-i)/double(n-1);
        for (uint_t k = 0; k < uint_t(n); ++k) {
            v[k] = i + k*dx;
        }

        return v;
    }
}

template<typename T, typename U, typename V>
vec1d rgen_log(T i, U j, V n) {
    phypp_check(n >= 0, "'rgen_log(a,b,n)' needs a positive or null value for 'n' (got ", n, ")");

    if (n == 1) {
        vec1d v = dindgen(1);
        v[0] = i;
        return v;
    } else {
        vec1d v = dindgen(n);

        double dx = log10(j/i)/double(n-1);
        for (uint_t k = 0; k < uint_t(n); ++k) {
            v[k] = i*e10(k*dx);
        }

        return v;
    }
}

template<typename T>
vec_t<2,T> make_bins(T mi, T ma) {
    return {{mi}, {ma}};
}

template<typename T>
vec_t<2,T> make_bins(T mi, T ma, uint_t n) {
    vec_t<2,T> b(2,n);
    T d = (ma - mi)/n;
    for (uint_t i = 0; i < n; ++i) {
        b(0,i) = mi + i*d;
        b(1,i) = mi + (i+1)*d;
    }

    return b;
}

template<typename T>
vec_t<2,T> make_bins(const vec_t<1,T>& v) {
    vec_t<2,T> b(2, v.size()-1);
    vec1u ids = uindgen(v.size()-1);
    b(0,_) = v[ids];
    b(1,_) = v[ids+1];

    return b;
}

template<std::size_t Dim, typename Type, typename B>
vec_t<Dim,bool> in_bin(const vec_t<Dim,Type>& v, const vec_t<1,B>& b) {
    return v >= b[0] && v < b[1];
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool finite(const T& t) {
    return std::isfinite(t);
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> finite(const vec_t<Dim,Type>& v) {
    vec_t<Dim,bool> r = boolarr(v.dims);
    for (uint_t i = 0; i < v.size(); ++i) {
        r[i] = std::isfinite(v[i]);
    }

    return r;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool nan(const T& t) {
    return std::isnan(t);
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> nan(const vec_t<Dim,Type>& v) {
    vec_t<Dim,bool> r = boolarr(v.dims);
    for (uint_t i = 0; i < v.size(); ++i) {
        r[i] = std::isnan(v[i]);
    }

    return r;
}

using seed_t = std::mt19937;

template<typename T>
seed_t make_seed(T seed) {
    return std::mt19937(seed);
}

template<typename T>
double randomn(T& seed) {
    std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(seed);
}

template<typename T, typename ... Args>
auto randomn(T& seed, Args&& ... args) {
    auto v = dblarr(std::forward<Args>(args)...);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (uint_t i = 0; i < v.size(); ++i) {
        v.data[i] = distribution(seed);
    }

    return v;
}

template<typename T>
double randomu(T& seed) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(seed);
}

template<typename T, typename ... Args>
auto randomu(T& seed, Args&& ... args) {
    auto v = dblarr(std::forward<Args>(args)...);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (uint_t i = 0; i < v.size(); ++i) {
        v.data[i] = distribution(seed);
    }

    return v;
}

template<typename T, typename TMi, typename TMa>
auto randomi(T& seed, TMi mi, TMa ma) -> decltype(mi + ma) {
    using rtype = decltype(mi + ma);
    std::uniform_int_distribution<rtype> distribution(mi, ma);
    return distribution(seed);
}

template<typename T, typename TMi, typename TMa, typename ... Args>
auto randomi(T& seed, TMi mi, TMa ma, Args&& ... args) {
    auto v = randomu(seed, std::forward<Args>(args)...);
    using rtype = decltype(mi + ma);
    return vec_t<vec_dim<decltype(v)>::value,rtype>(v*(ma + 1 - mi) + mi);
}

template<std::size_t Dim, typename Type, typename T>
auto shuffle(vec_t<Dim,Type> v, T& seed) {
    std::shuffle(v.begin(), v.end(), seed);
    return v;
}

template<typename F, F f, std::size_t Dim, typename Type>
auto run_index_(const vec_t<Dim,Type>& v, uint_t dim) -> vec_t<Dim-1,typename return_type<F>::type> {
    vec_t<Dim-1,typename return_type<F>::type> r;
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

//  Example demonstration of the index computation:
//
//  Assume we have a 4D array of dimensions d1, d2, d3 and d4.
//  We want to compute an operation on the array generated by:
//      data[i*d2*d3*d4 + j*d3*d4 + k*d4 + l]
//  and by running the 'dim'th index (i, j, k or l), assigning
//  the value to another 3D array using a continuous index 'u'.
//
//  - (i,j,k,_) : k+j*d3+i*d2*d3 = u : mean(((u/1)*d4 + l)*1 + (u%1))
//      mpitch = 1
//
//  - (i,j,_,l) : l+j*d4+i*d2*d4 = u : mean(((u/d4)*d3 + k)*d4 + (u%d4))
//      mpitch = d4
//
//  - (i,_,k,l) : l+k*d4+i*d3*d4 = u : mean(((u/(d3*d4))*d2 + j)*d3*d4 + (u%(d3*d4)))
//      mpitch = d4*d3
//
//  - (_,j,k,l) : l+k*d4+j*d3*d4 = u : mean(((u/(d4*d3*d2))*d1 + i)*d2*d3*d4 + (u%(d4*d3*d2)))
//      mpitch = d4*d3*d2
//
//  Final recipe:
//      ((u/mpitch)*dim[d] + i)*mpitch + (u%mpitch)

    vec_t<1,rtype_t<Type>> tmp(v.dims[dim]);
    for (uint_t i = 0; i < r.size(); ++i) {
        uint_t base = (i%mpitch) + (i/mpitch)*v.dims[dim]*mpitch;
        for (uint_t j = 0; j < v.dims[dim]; ++j) {
            tmp[j] = dref(v.data[base + j*mpitch]);
        }

        r[i] = (*f)(tmp);
    }

    return r;
}

template<typename F, std::size_t Dim, typename Type>
void run_dim_idx(const vec_t<Dim,Type>& v, uint_t dim, F&& func) {
    uint_t nint = v.dims[dim];
    uint_t np = v.size()/nint;
    uint_t mpitch = 1;
    for (uint_t i = dim+1; i < Dim; ++i) {
        mpitch *= v.dims[i];
    }

    vec_t<1,rtype_t<Type>> tmp(nint);
    for (uint_t i = 0; i < np; ++i) {
        uint_t base = (i%mpitch) + (i/mpitch)*nint*mpitch;
        for (uint_t j = 0; j < nint; ++j) {
            tmp[j] = dref(v.data[base + j*mpitch]);
        }

        func(i, tmp);
    }
}

template<std::size_t Dim, typename Type>
vec_t<1,rtype_t<Type>> run_dim_idx_apply_ids_(const vec1u& ids, const vec_t<Dim,Type>& v) {
    return v[ids].concretise();
}

template<std::size_t Dim, typename Type, typename ... Args>
std::array<uint_t,Dim> run_dim_idx_get_dim_(const vec_t<Dim,Type>& v, const Args& ... vs) {
    return v.dims;
}

template<typename F, typename ... Args>
void run_dim_idx(uint_t dim, F&& func, const Args& ... vs) {
    auto ds = run_dim_idx_get_dim_(vs...);
    const uint_t N = array_size<decltype(ds)>::size;

    uint_t nint = ds[dim];
    uint_t mpitch = 1;
    uint_t np = 1;
    for (uint_t i = 0; i < N; ++i) {
        if (i != dim) np *= ds[i];
        if (i > dim) mpitch *= ds[i];
    }

    vec1u ids(nint);
    for (uint_t i = 0; i < np; ++i) {
        uint_t base = (i%mpitch) + (i/mpitch)*nint*mpitch;
        for (uint_t j = 0; j < nint; ++j) {
            ids[j] = base + j*mpitch;
        }

        func(i, run_dim_idx_apply_ids_<N>(ids, vs)...);
    }
}

template<typename F, std::size_t Dim, typename Type>
auto run_dim(const vec_t<Dim,Type>& v, uint_t dim, F&& func) -> vec_t<Dim-1,typename return_type<F>::type> {
    vec_t<Dim-1,typename return_type<F>::type> r;
    for (uint_t i = 0; i < dim; ++i) {
        r.dims[i] = v.dims[i];
    }
    for (uint_t i = dim+1; i < Dim; ++i) {
        r.dims[i-1] = v.dims[i];
    }
    r.resize();

    run_dim_idx(v, dim, [&](uint_t i, const vec_t<1, rtype_t<Type>>& tv) {
        r[i] = func(tv);
    });

    return r;
}

// template<typename F, std::size_t Dim, typename Type>
// auto run_dim(const vec_t<Dim,Type>& v, uint_t dim, F&& f) -> vec_t<Dim-1,typename return_type<F>::type> {
//     vec_t<Dim-1,typename return_type<F>::type> r;
//     for (uint_t i = 0; i < dim; ++i) {
//         r.dims[i] = v.dims[i];
//     }
//     for (uint_t i = dim+1; i < Dim; ++i) {
//         r.dims[i-1] = v.dims[i];
//     }
//     r.resize();

//     uint_t mpitch = 1;
//     for (uint_t i = dim+1; i < Dim; ++i) {
//         mpitch *= v.dims[i];
//     }

//     typename vec_t<1,Type>::effective_type tmp = arr<rtype_t<Type>>(v.dims[dim]);
//     for (uint_t i = 0; i < r.size(); ++i) {
//         for (uint_t j = 0; j < v.dims[dim]; ++j) {
//             tmp[j] = dref(v.data[(i%mpitch) + ((i/mpitch)*v.dims[dim] + j)*mpitch]);
//         }

//         r[i] = f(tmp);
//     }

//     return r;
// }

template<std::size_t Dim, typename Type>
double total(const vec_t<Dim,Type>& v) {
    double total = 0;
    for (auto& t : v) {
        total += t;
    }

    return total;
}


template<std::size_t Dim, typename Type>
vec_t<Dim-1,double> total(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = double (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &total<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim>
double fraction(const vec_t<Dim,bool>& b) {
    return total(b)/b.size();
}

template<std::size_t Dim, typename Type>
double mean(const vec_t<Dim,Type>& v) {
    double total = 0.0;
    for (auto& t : v) {
        total += t;
    }

    return total/n_elements(v);
}

template<std::size_t Dim, typename Type>
vec_t<Dim-1,double> mean(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = double (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &mean<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> median(const vec_t<Dim,Type>& v) {
    vec1u ok = where(finite(v));
    if (ok.empty()) return 0;

    typename vec_t<1,Type>::effective_type t = v[ok];
    std::ptrdiff_t offset = t.size()/2;
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    return *(t.begin() + offset);
}

template<std::size_t Dim, typename Type>
vec_t<Dim-1,rtype_t<Type>> median(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = rtype_t<Type> (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &median<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim, typename Type, typename U>
rtype_t<Type> percentile(const vec_t<Dim,Type>& v, const U& u) {
    vec1u ok = where(finite(v));
    if (ok.empty()) return 0;

    typename vec_t<1,Type>::effective_type t = v[ok];
    std::ptrdiff_t offset = t.size()*u;
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    return *(t.begin() + offset);
}

template<std::size_t Dim, typename Type>
void percentiles_(vec_t<1,Type>& r, uint_t i, vec_t<Dim,Type>& t) {}

template<std::size_t Dim, typename Type, typename U, typename ... Args>
void percentiles_(vec_t<1,Type>& r, uint_t i, vec_t<Dim,Type>& t, const U& u, const Args& ... args) {
    std::ptrdiff_t offset = t.size()*u;
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    r.data[i] = *(t.begin() + offset);
    ++i;

    percentiles_(r, i, t, args...);
}

template<std::size_t Dim, typename Type, typename ... Args>
typename vec_t<1,Type>::effective_type percentiles(const vec_t<Dim,Type>& v, const Args& ... args) {
    vec1u ok = where(finite(v));
    typename vec_t<1,Type>::effective_type t;
    if (ok.empty()) return t;
    t = v[ok];

    typename vec_t<1,Type>::effective_type r = arr<rtype_t<Type>>(sizeof...(Args));
    percentiles_(r, 0, t, args...);

    return r;
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> sigma_clip(const vec_t<Dim,Type>& v, double percl = 0.15, double percu = -1.0) {
    if (percl > 0.5) percl = 1 - percl;
    if (percu < 0.0) percu = 1 - percl;
    else if (percu < percl) std::swap(percu, percl);

    auto p = percentiles(v, percl, percu);
    return v > p[0] && v < p[1];
}

template<std::size_t Dim, typename Type>
rtype_t<Type> min(const vec_t<Dim,Type>& v) {
    return dref(*std::min_element(v.data.begin(), v.data.end(), [](Type t1, Type t2){
        if (nan(dref(t1))) return false;
        if (nan(dref(t2))) return true;
        return dref(t1) < dref(t2);
    }));
}

template<std::size_t Dim, typename Type>
rtype_t<Type> max(const vec_t<Dim,Type>& v) {
    return dref(*std::max_element(v.data.begin(), v.data.end(), [](Type t1, Type t2){
        if (nan(dref(t1))) return true;
        if (nan(dref(t2))) return false;
        return dref(t1) < dref(t2);
    }));
}

template<std::size_t Dim, typename Type>
uint_t min_id(const vec_t<Dim,Type>& v) {
    vec1u tmp = uindgen(v.size());
    return *std::min_element(tmp.begin(), tmp.end(), [&](uint_t i1, uint_t i2){
        if (nan(v[i1])) return false;
        if (nan(v[i2])) return true;
        return v[i1] < v[i2];
    });
}

template<std::size_t Dim, typename Type>
uint_t max_id(const vec_t<Dim,Type>& v) {
    vec1u tmp = uindgen(v.size());
    return *std::max_element(tmp.begin(), tmp.end(), [&](uint_t i1, uint_t i2){
        if (nan(v[i1])) return true;
        if (nan(v[i2])) return false;
        return v[i1] < v[i2];
    });
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const vec_t<Dim,Type1>& v1, const vec_t<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "min: incompatible vector dimensions");
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r[i] = std::min(v1[i], v2[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const vec_t<Dim,Type1>& v1, const vec_t<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "max: incompatible vector dimensions");
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r[i] = std::max(v1[i], v2[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const vec_t<Dim,Type1>& v1, const Type2& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r[i] = std::min(v1[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const vec_t<Dim,Type1>& v1, const Type2& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r[i] = std::max(v1[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const Type1& v1, const vec_t<Dim,Type2>& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r[i] = std::min(v1, v2[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const Type1& v1, const vec_t<Dim,Type2>& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r[i] = std::max(v1, v2[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type>
double rms(const vec_t<Dim,Type>& v) {
    double sum = 0;
    for (auto& t : v) {
        sum += t*t;
    }

    return sqrt(sum/v.size());
}

template<std::size_t Dim, typename Type>
vec_t<Dim-1,double> rms(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = double (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &rms<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim, typename Type>
double stddev(const vec_t<Dim,Type>& v) {
    return rms(v - mean(v));
}

template<std::size_t Dim, typename Type>
vec_t<Dim-1,double> stddev(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = double (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &stddev<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim, typename Type>
double mad(const vec_t<Dim,Type>& v) {
    return median(fabs(v - median(v)));
}

template<std::size_t Dim, typename Type>
vec_t<Dim-1,double> mad(const vec_t<Dim,Type>& v, uint_t dim) {
    using fptr = double (*)(const vec_t<1,rtype_t<Type>>&);
    return run_index_<fptr, &mad<1,rtype_t<Type>>>(v, dim);
}

template<std::size_t Dim, typename Type, typename TypeB>
vec1u histogram(const vec_t<Dim,Type>& data, const vec_t<2,TypeB>& bins) {
    using rtype = dtype_t<Type>;
    uint_t nbin = bins.dims[1];
    vec1u counts(nbin);
    for (uint_t i = 0; i < nbin; ++i) {
        counts[i] = total(in_bin(data, bins(_,i)));
    }

    return counts;
}

template<std::size_t Dim, typename Type, typename TypeB, typename TypeW>
vec1d histogram(const vec_t<Dim,Type>& data, const vec_t<Dim,TypeW>& weight,
    const vec_t<2,TypeB>& bins) {

    using rtype = dtype_t<Type>;
    uint_t nbin = bins.dims[1];
    vec1d counts(nbin);
    for (uint_t i = 0; i < nbin; ++i) {
        vec1u ids = where(in_bin(data, bins(_,i)));
        counts[i] = total(weight[ids]);
    }

    return counts;
}

template<std::size_t Dim, typename Type>
void data_info_(const vec_t<Dim,Type>& v) {
    vec1u idok = where(finite(v));
    print(n_elements(idok), "/", n_elements(v), " valid values (dims: ", dim(v), ")");
    if (n_elements(idok) == 0) return;
    print(" min : ", min(v[idok]));
    print(" 15% : ", percentile(v[idok], 0.15));
    print(" 50% : ", median(v[idok]));
    print(" mean: ", mean(v[idok]));
    print(" 85% : ", percentile(v[idok], 0.85));
    print(" max : ", max(v[idok]));
    print(" rms : ", stddev(v[idok] - median(v[idok])));
}

#define data_info(x) \
    print("data info: ", #x); \
    data_info_(x);

template<typename T, typename U, typename V,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
auto clamp(const T& t, const U& mi, const V& ma) {
    return (t < mi ? mi : (t > ma ? ma : t));
}

template<std::size_t Dim, typename Type, typename U,
    typename enable = typename std::enable_if<!is_vec<U>::value>::type>
auto pow(const U& u, const vec_t<Dim,Type>& v) -> vec_t<Dim, decltype(pow(u,v(0)))> {
    vec_t<Dim, decltype(pow(u,v(0)))> r = v;
    for (auto& t : r) {
        t = pow(u, t);
    }
    return r;
}

template<std::size_t Dim, typename Type, typename U>
auto pow(const U& u, vec_t<Dim,Type>&& v) -> typename std::enable_if<!is_vec<U>::value &&
    !std::is_pointer<Type>::value && std::is_same<decltype(pow(u,v(0))), Type>::value,
    vec_t<Dim,Type>>::type {
    for (auto& t : v) {
        t = pow(u, t);
    }
    return std::move(v);
}

template<typename T>
auto sqr(T&& t) {
    return t*t;
}

template<std::size_t Dim, typename Type, typename enable =
    typename std::enable_if<!std::is_pointer<Type>::value>::type>
auto sqr(vec_t<Dim,Type>&& v) {
    for (auto& t : v) {
        t *= t;
    }
    return std::move(v);
}

template<typename T>
auto invsqr(T&& t) {
    return 1.0/(t*t);
}

template<std::size_t Dim, typename Type, typename enable =
    typename std::enable_if<!std::is_pointer<Type>::value>::type>
auto invsqr(vec_t<Dim,Type>&& v) {
    for (auto& t : v) {
        t = 1.0/(t*t);
    }
    return std::move(v);
}

#define VECTORIZE(name) \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(const vec_t<Dim,Type>& v, const Args& ... args) { \
        using ntype = decltype(name(v[0], args...)); \
        vec_t<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(name(dref(t), args...)); \
        } \
        return r; \
    } \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(vec_t<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<Type>::value && std::is_same<decltype(name(v[0], args...)), Type>::value, \
        vec_t<Dim,Type>>::type { \
        for (auto& t : v) { \
            t = name(t, args...); \
        } \
        return std::move(v); \
    }

#define VECTORIZE_REN(name, orig) \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(const vec_t<Dim,Type>& v, const Args& ... args) { \
        using ntype = decltype(orig(v[0], args...)); \
        vec_t<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(orig(dref(t), args...)); \
        } \
        return r; \
    } \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(vec_t<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<Type>::value && std::is_same<decltype(orig(v[0], args...)), Type>::value, \
        vec_t<Dim,Type>>::type { \
        for (auto& t : v) { \
            t = orig(t, args...); \
        } \
        return std::move(v); \
    } \
    template<typename ... Args> \
    auto name(Args&& ... args) { \
        return orig(std::forward<Args>(args)...); \
    }

VECTORIZE(sqrt);
VECTORIZE(sqr);
VECTORIZE(invsqr);
VECTORIZE(pow);
VECTORIZE(cos);
VECTORIZE(sin);
VECTORIZE(tan);
VECTORIZE(acos);
VECTORIZE(asin);
VECTORIZE(atan);
VECTORIZE(cosh);
VECTORIZE(sinh);
VECTORIZE(tanh);
VECTORIZE(acosh);
VECTORIZE(asinh);
VECTORIZE(atanh);
VECTORIZE(exp);
VECTORIZE(log);
VECTORIZE(log2);
VECTORIZE(log10);
VECTORIZE(erf);
VECTORIZE(erfc);
VECTORIZE(tgamma);
VECTORIZE(ceil);
VECTORIZE(floor);
VECTORIZE(round);
VECTORIZE(fabs);
VECTORIZE(clamp);
VECTORIZE_REN(bessel_j0, j0);
VECTORIZE_REN(bessel_j1, j1);
VECTORIZE_REN(bessel_y0, y0);
VECTORIZE_REN(bessel_y1, y1);

#ifndef NO_GSL
VECTORIZE_REN(bessel_i0, gsl_sf_bessel_I0);
VECTORIZE_REN(bessel_i1, gsl_sf_bessel_I1);
VECTORIZE_REN(bessel_k0, gsl_sf_bessel_K0);
VECTORIZE_REN(bessel_k1, gsl_sf_bessel_K1);
#endif

#undef VECTORIZE

template<typename F>
auto derivate1(F func, const double& x, const double ep) -> decltype(func(x)) {
    static const double a[5] = {1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0};

    double tmp = x - 2*ep;
    decltype(func(x)) res = a[0]*func(tmp)/ep;
    for (uint_t i = 1; i < 5; ++i) {
        tmp += ep;
        res += a[i]*func(tmp)/ep;
    }

    return res;
}

template<typename F>
auto derivate2(F func, const double& x, const double ep) -> decltype(func(x)) {
    static const double a[5] = {-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0};

    double tmp = x - 2*ep;
    decltype(func(x)) res = a[0]*func(tmp)/(ep*ep);
    for (uint_t i = 1; i < 5; ++i) {
        tmp += ep;
        res += a[i]*func(tmp)/(ep*ep);
    }

    return res;
}

template<typename F>
auto derivate1(F func, const vec1d& x, const double ep, uint_t ip) -> decltype(func(x)) {
    static const double a[5] = {1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0};

    vec1d tmp = x;
    tmp[ip] -= 2*ep;

    decltype(func(x)) res = a[0]*func(tmp)/ep;
    for (uint_t i = 1; i < 5; ++i) {
        tmp[ip] += ep;
        res += a[i]*func(tmp)/ep;
    }

    return res;
}

template<typename F>
auto derivate2(F func, const vec1d& x, const double ep, uint_t ip) -> decltype(func(x)) {
    static const double a[5] = {-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0};

    vec1d tmp = x;
    tmp[ip] -= 2*ep;

    decltype(func(x)) res = a[0]*func(tmp)/(ep*ep);
    for (uint_t i = 1; i < 5; ++i) {
        tmp[ip] += ep;
        res += a[i]*func(tmp)/(ep*ep);
    }

    return res;
}

template<typename F>
auto derivate2(F func, const vec1d& x, const double ep, uint_t ip1, uint_t ip2) -> decltype(func(x)) {
    static const double a[5] = {1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0};

    vec1d tmp = x;
    tmp[ip1] -= 2*ep;
    tmp[ip2] -= 2*ep;

    decltype(func(x)) res = a[0]*a[0]*func(tmp)/(ep*ep);
    for (uint_t j = 1; j < 5; ++j) {
        tmp[ip2] += ep;
        res += a[0]*a[j]*func(tmp)/(ep*ep);
    }

    for (uint_t i = 1; i < 5; ++i) {
        tmp[ip1] += ep;
        tmp[ip2] = x(ip2) - 2*ep;

        res += a[i]*a[0]*func(tmp)/(ep*ep);
        for (uint_t j = 1; j < 5; ++j) {
            tmp[ip2] += ep;
            res += a[i]*a[j]*func(tmp)/(ep*ep);
        }
    }

    return res;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec_t<2,TypeA>& a, const vec_t<2,TypeB>& b) -> vec_t<2,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying matrices "
        "("+strn(a.dims)+" x "+strn(b.dims)+")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];
    const uint_t m = b.dims[1];

    vec_t<2,ntype_t> r(n,m);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t j = 0; j < m; ++j)
    for (uint_t k = 0; k < o; ++k) {
        r(i,j) += a(i,k)*b(k,j);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec_t<2,TypeA>& a, const vec_t<1,TypeB>& b) -> vec_t<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying matrix by vector "
        "("+strn(a.dims)+" x "+strn(b.dims)+")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];

    vec_t<1,ntype_t> r(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r(i) += a(i,k)*b(k);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec_t<1,TypeB>& b, const vec_t<2,TypeA>& a) -> vec_t<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying vector by matrix "
        "("+strn(a.dims)+" x "+strn(b.dims)+")");

    const uint_t o = a.dims[0];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[1];

    vec_t<1,ntype_t> r = arr<ntype_t>(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r(i) += b(k)*a(k,i);
    }

    return r;
}

template<typename Type>
void mprint(const vec_t<2,Type>& m) {
    for (uint_t i = 0; i < m.dims[0]; ++i) {
        for (uint_t j = 0; j < m.dims[1]; ++j) {
            if (j != 0) std::cout << ", ";
            std::cout << m(i,j);
        }

        std::cout << "\n";
    }

    std::cout << std::flush;
}

template<std::size_t Dim, typename Type, typename ... Args>
void transpose_(vec_t<Dim,rtype_t<Type>>& r, const vec_t<Dim,Type>& v, cte_t<Dim>, Args ... args) {
    auto t = std::make_tuple(args...);
    r[tuple_reverse(t)] = v[t];
}

template<std::size_t Dim, typename Type, std::size_t I, typename ... Args>
void transpose_(vec_t<Dim,rtype_t<Type>>& r, const vec_t<Dim,Type>& v, cte_t<I>, Args ... args) {
    for (uint_t i = 0; i < v.dims[I]; ++i) {
        transpose_(r, v, cte_t<I+1>(), args..., i);
    }
}

template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<(Dim>1)>::type>
vec_t<Dim,rtype_t<Type>> transpose(const vec_t<Dim,Type>& v) {
    auto d = v.dims;
    std::reverse(d.begin(), d.end());
    vec_t<Dim,rtype_t<Type>> r(d);
    transpose_(r, v, cte_t<0>());
    return r;
}

template<typename Type>
auto diag(const vec_t<2,Type>& v) -> decltype(v(_,0)) {
    assert(v.dims[0] == v.dims[1]);
    decltype(v(_,0)) d(get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i = 0; i < v.dims[0]; ++i) {
        d.data[i] = &v(i,i);
    }

    return d;
}

template<typename Type>
auto diag(vec_t<2,Type>& v) -> decltype(v(_,0)) {
    assert(v.dims[0] == v.dims[1]);
    decltype(v(_,0)) d(get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i = 0; i < v.dims[0]; ++i) {
        d.data[i] = &v(i,i);
    }

    return d;
}

template<typename Type = double>
auto identity_matrix(uint_t dim) {
    auto m = arr<Type>(dim, dim);
    diag(m) = 1;
    return m;
}

template<typename TX, typename TY>
auto scale_matrix(const TX& sx, const TY& sy) {
    auto m = arr<decltype(sx*sy)>(3, 3);
    m(0,0) = sx;
    m(1,1) = sy;
    m(2,2) = 1;
    return m;
}

template<typename T>
auto scale_matrix(const T& s) {
    auto m = arr<T>(3, 3);
    m(0,0) = s;
    m(1,1) = s;
    m(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto translation_matrix(const TX& tx, const TY& ty) {
    auto m = arr<decltype(tx*ty)>(3, 3);
    diag(m) = 1;
    m(0,2) = tx;
    m(1,2) = ty;
    return m;
}

template<typename A>
auto rotation_matrix(const A& a) {
    auto m = arr<decltype(cos(a))>(3, 3);
    auto ca = cos(a), sa = sin(a);
    m(0,0) = m(1,1) = ca;
    m(0,1) = -sa;
    m(1,0) = sa;
    m(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto point2d(const TX& x, const TY& y) {
    return vec_t<1,decltype(x*y)>{x, y, 1};
}

extern "C" void dgetrf_(int* n, int* m, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);

bool invert(vec2d& i) {
    phypp_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

    int n = i.dims[0];
    int lda = n;
    int info;

    vec_t<1,int> ipiv(n);
    dgetrf_(&n, &n, i.data.data(), &lda, ipiv.data.data(), &info);
    if (info < 0) {
        return false;
    }

    vec1d work(n);
    int lw = n;
    dgetri_(&n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &lw, &info);
    if (info != 0) {
        return false;
    }

    return true;
}

template<typename TypeA>
bool invert(const vec_t<2,TypeA>& a, vec2d& i) {
    i = a;
    return invert(i);
}

extern "C" void dsytrf_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work,
    int* lwork, int* info);
extern "C" void dsytri_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work, int* info);

bool invert_symmetric(vec2d& i) {
    phypp_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

    char uplo = 'U';
    int n = i.dims[0];
    int lda = n;
    int info;

    int lw = n*64;
    // Note: the optimal value for lw is n*nb, where nb is the optimal block size
    // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
    // the Lapack User Guide.

    vec1d work(lw);
    vec_t<1,int> ipiv(n);

    dsytrf_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &lw, &info);
    if (info < 0) {
        return false;
    }

    dsytri_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &info);
    if (info != 0) {
        return false;
    }

    return true;
}

template<typename TypeA>
bool invert_symmetric(const vec_t<2,TypeA>& a, vec2d& i) {
    i = a;
    return invert_symmetric(i);
}

extern "C" void dsysv_(char* uplo, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,
    int* ldb, double* work, int* lwork, int* info);

bool solve_symmetric(vec2d& alpha, vec1d& beta) {
    phypp_check(alpha.dims[0] == alpha.dims[1], "cannot invert a non square matrix (",
        alpha.dims, ")");
    phypp_check(alpha.dims[0] == beta.dims[0], "matrix and vector must have the same dimesions (",
        "got ", alpha.dims[0], " and ", beta.dims[0], ")");

    char uplo = 'U';
    int n = alpha.dims[0];
    int nrhs = 1;
    int lda = n, ldb = n;
    int info;

    int lw = n*64;
    // Note: the optimal value for lw is n*nb, where nb is the optimal block size
    // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
    // the Lapack User Guide.

    vec1d work(lw);
    vec_t<1,int> ipiv(n);

    dsysv_(&uplo, &n, &nrhs, alpha.data.data(), &lda, ipiv.data.data(), beta.data.data(),
        &ldb, work.data.data(), &lw, &info);
    if (info != 0) {
        return false;
    }

    return true;
}

template<typename Type>
void symmetrize(vec_t<2,Type>& alpha) {
    phypp_check(alpha.dims[0] == alpha.dims[1], "cannot symmetrize a non square matrix (",
        alpha.dims, ")");

    for (uint_t i = 0; i < alpha.dims[0]; ++i)
    for (uint_t j = i+1; j < alpha.dims[0]; ++j) {
        alpha(i,j) = alpha(j,i);
    }
}

enum nlfit_status {
    ILL_MODEL,
    MAX_ITER,
    SINGULAR_HESSIAN,
    MIN_REL,
    MIN_ABS,
    PERFECT_FIT
};

struct nlfit_result {
    bool         success;
    nlfit_status status;

    double chi2;
    vec1d  params;
    vec1d  errors;
    vec2d  cov;
    uint_t niter;
};

template<typename Func, typename TypeY, typename TypeE>
nlfit_result nlfit(Func f, const TypeY& y, const TypeE& ye, vec1d params, double e = 0.001,
    double abs = 0.01, const uint_t max_iter = 1000) {

    nlfit_result fr;

    uint_t np = n_elements(params);
    uint_t nm = n_elements(y);

    vec2d alpha = dblarr(np,np);
    vec1d beta  = dblarr(np);

    vec1d best_params = params;
    vec2d best_alpha  = alpha;

    auto make_chi2 = [&](const vec1d& p) {
        return total(sqr((f(p) - y)/ye));
    };

    double best_chi2 = make_chi2(params);
    double prev_chi2 = best_chi2;
    fr.chi2          = best_chi2;

    double lambda = 0.001;
    double factor = 10.0;

    fr.errors = dblarr(np);
    fr.niter = 0;

    vec2d deriv = dblarr(np,nm);

    auto make_alpha = [&]() {
        for (uint_t i = 0; i < np; ++i) {
            deriv(i,_) = flatten(derivate1(f, params, 0.001, i)/ye);
        }

        auto tmp = flatten((y - f(params))/ye);
        for (uint_t i = 0; i < np; ++i) {
            for (uint_t j = 0; j < np; ++j) {
                if (i == j) {
                    alpha(i,i) = (1.0 + lambda)*total(deriv(i,_)*deriv(i,_));
                } else if (i < j) {
                    alpha(i,j) = total(deriv(i,_)*deriv(j,_));
                } else {
                    alpha(i,j) = alpha(j,i);
                }
            }

            beta(i) = total(deriv(i,_)*tmp);
        }

        if (!invert(alpha)) {
            return false;
        } else {
            return true;
        }
    };

    auto fail = [&](nlfit_status s) {
        params = best_params;
        if (make_alpha()) {
            best_alpha = alpha;
        }

        auto d = diag(best_alpha);
        d /= (1.0 + lambda);

        fr.cov     = best_alpha;
        fr.chi2    = best_chi2;
        fr.params  = params;
        fr.errors  = sqrt(d);
        fr.success = false;
        fr.status  = s;
    };

    auto succeed = [&](nlfit_status s) {
        if (make_alpha()) {
            best_alpha = alpha;
        }

        auto d = diag(best_alpha);
        d /= (1.0 + lambda);

        fr.cov     = best_alpha;
        fr.params  = params;
        fr.errors  = sqrt(d);
        fr.success = true;
        fr.status  = s;
    };

    while (max_iter == 0 || fr.niter < max_iter) {
        if (!make_alpha()) {
            fail(SINGULAR_HESSIAN);
            return fr;
        }

        params += mmul(alpha, beta);
        fr.chi2 = make_chi2(params);

        if (!finite(fr.chi2)) {
            fail(ILL_MODEL);
            return fr;
        }

        if (fr.chi2 == 0.0) {
            succeed(PERFECT_FIT);
            return fr;
        }

        if (fr.chi2 >= best_chi2) {
            lambda *= factor;

            prev_chi2 = fr.chi2;
            params = best_params;
        } else {
            if (fr.niter > 10) {
                if ((best_chi2 - fr.chi2)/best_chi2 < e) {
                    succeed(MIN_REL);
                    return fr;
                }

                if (abs != 0.0 && fabs(prev_chi2 - fr.chi2) < abs) {
                    succeed(MIN_ABS);
                    return fr;
                }
            }

            lambda /= factor;

            prev_chi2   = best_chi2 = fr.chi2;
            best_params = params;
            best_alpha  = alpha;
        }

        ++fr.niter;
    }

    fail(MAX_ITER);
    return fr;
}

struct linfit_result {
    bool success;

    double chi2;
    vec1d  params;
    vec1d  errors;
    vec2d  cov;
};

template<typename T, typename TypeE>
void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t) {
    cache(i,_) = flatten(t/ye);
}

template<typename T, typename TypeE, typename ... Args>
void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t, Args&& ... args) {
    cache(i,_) = flatten(t/ye);
    linfit_make_cache_(cache, ye, i+1, args...);
}

template<typename TypeY, typename TypeE>
linfit_result linfit_do_(const TypeY& y, const TypeE& ye, const vec2d& cache) {
    linfit_result fr;

    uint_t np = cache.dims[0];
    uint_t nm = cache.dims[1];

    // Solving 'y +/- e = sum over i of a[i]*x[i]' to get all a[i]'s
    vec2d alpha(np,np);
    vec1d beta(np);
    auto tmp = flatten(y/ye);
    for (uint_t i = 0; i < np; ++i) {
        for (uint_t j = 0; j < np; ++j) {
            if (i <= j) {
                alpha(i,j) = 0.0;
                // alpha(i,j) = sum over all points of x[i]*x[j]/e^2
                for (uint_t m = 0; m < nm; ++m) {
                    alpha(i,j) += cache(i,m)*cache(j,m);
                }
            } else {
                alpha(i,j) = alpha(j,i);
            }
        }

        beta[i] = 0.0;
        // beta[i] = sum over all points of x[i]*y/e^2
        for (uint_t m = 0; m < nm; ++m) {
            beta[i] += cache(i,m)*tmp[m];
        }
    }

    if (!invert_symmetric(alpha)) {
        fr.success = false;
        fr.chi2 = dnan;
        fr.params = replicate(dnan, np);
        fr.errors = replicate(dnan, np);
        symmetrize(alpha);
        fr.cov = alpha;
        return fr;
    }

    symmetrize(alpha);
    fr.success = true;
    fr.params = mmul(alpha, beta);
    fr.errors = sqrt(diag(alpha));

    vec1d model(nm);
    for (uint_t m = 0; m < nm; ++m) {
        model[m] = 0.0;
        for (uint_t i = 0; i < np; ++i) {
            model[m] += fr.params[i]*cache(i,m);
        }
    }

    fr.chi2 = total(sqr(model*flatten(1.0/ye) - tmp));

    return fr;
}

template<typename TypeY, typename TypeE, typename ... Args>
linfit_result linfit(const TypeY& y, const TypeE& ye, Args&&... args) {
    uint_t np = sizeof...(Args);
    uint_t nm = n_elements(y);

    vec2d cache(np,nm);
    linfit_make_cache_(cache, ye, 0, std::forward<Args>(args)...);

    return linfit_do_(y, ye, cache);
}

template<std::size_t Dim, typename TypeY, typename TypeE, typename TypeX>
linfit_result linfit_pack(const vec_t<Dim,TypeY>& y, const vec_t<Dim,TypeE>& ye,
    const vec_t<Dim+1,TypeX>& x) {

    linfit_result fr;

    uint_t np = x.dims[0];
    uint_t nm = n_elements(y);

    vec2d cache(np,nm);
    for (uint_t i = 0; i < np; ++i) {
        for (uint_t j = 0; j < nm; ++j) {
            cache(i,j) = x[i*x.pitch(0) + j]/ye[j];
        }
    }

    return linfit_do_(y, ye, cache);
}

struct affinefit_result {
    bool success;

    double chi2;
    double slope, offset;
    double slope_err, offset_err;
};

// Perform a simple linear fit 'y = offset + slope*x'
template<typename TypeX, typename TypeY, std::size_t Dim, typename TypeE>
affinefit_result affinefit(const TypeX& x, const TypeY& y, const vec_t<Dim,TypeE>& ye) {
    affinefit_result fr;
    fr.success = true;

    auto s = total(1.0/ye);
    auto sx = total(x/ye);
    auto sy = total(y/ye);
    auto sxx = total(x*x/ye);
    auto sxy = total(x*y/ye);

    auto delta = s*sxx - sx*sx;
    fr.offset = (sxx*sy - sx*sxy)/delta;
    fr.slope = (s*sxy - sx*sy)/delta;
    fr.offset_err = sqrt(sxx/delta);
    fr.slope_err = sqrt(s/delta);

    return fr;
}

// Perform a simple linear fit 'y = offset + slope*x'
template<typename TypeX, typename TypeY, typename TypeE>
affinefit_result affinefit(const TypeX& x, const TypeY& y, const TypeE& ye) {
    affinefit_result fr;
    fr.success = true;

    auto s = n_elements(x)/ye;
    auto sx = total(x/ye);
    auto sy = total(y/ye);
    auto sxx = total(x*x/ye);
    auto sxy = total(x*y/ye);

    auto delta = s*sxx - sx*sx;
    fr.offset = (sxx*sy - sx*sxy)/delta;
    fr.slope = (s*sxy - sx*sy)/delta;
    fr.offset_err = sqrt(sxx/delta);
    fr.slope_err = sqrt(s/delta);

    return fr;
}

// Returns the position of the first value in the array that is less than or equal to 'x'.
// Returns 'npos' if no value satisfy this criterium.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, typename Type>
uint_t lower_bound(T x, const vec_t<1,Type>& v) {
    auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
        typename vec_t<1,Type>::comparator());

    if (iter == v.data.begin()) {
        return npos;
    } else {
        return iter - v.data.begin() - 1;
    }
}

// Returns the position of the first value in the array that is greater than 'x'.
// Returns 'npos' if no value satisfy this criterium.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, typename Type>
uint_t upper_bound(T x, const vec_t<1,Type>& v) {
    auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
        typename vec_t<1,Type>::comparator());

    if (iter == v.data.end()) {
        return npos;
    } else {
        return iter - v.data.begin();
    }
}

// Check if a given array is sorted or not
template<typename Type>
bool is_sorted(const vec_t<1,Type>& v) {
    for (uint_t i = 0; i < v.size()-1; ++i) {
        if (v[i] >= v[i+1]) return false;
    }

    return true;
}

// Perform linear interpolation of data 'y' of position 'x' at new positions 'nx'.
// Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
// the arrays contains special values (NaN, inf, ...), all the points that would use these values
// will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
template<typename TypeX2 = double, typename TypeY = double, typename TypeX1 = double>
auto interpolate(const vec_t<1,TypeY>& y, const vec_t<1,TypeX1>& x, const vec_t<1,TypeX2>& nx) {
    using rtypey = rtype_t<TypeY>;
    using rtypex = rtype_t<TypeX1>;

    phypp_check(y.size() == x.size(),
        "interpolate: 'x' and 'y' arrays must contain the same number of elements");
    phypp_check(y.size() >= 2,
        "interpolate: 'x' and 'y' arrays must contain at least 2 elements");

    uint_t nmax = x.size();
    vec_t<1,decltype(y[0]*x[0])> r; r.reserve(nx.size());
    for (auto& tx : nx) {
        uint_t low = lower_bound(tx, x);

        rtypey ylow, yup;
        rtypex xlow, xup;
        if (low != npos) {
            if (low != nmax-1) {
                ylow = dref(y.data[low]); yup = dref(y.data[low+1]);
                xlow = dref(x.data[low]); xup = dref(x.data[low+1]);
            } else {
                ylow = dref(y.data[low-1]); yup = dref(y.data[low]);
                xlow = dref(x.data[low-1]); xup = dref(x.data[low]);
            }
        } else {
            ylow = dref(y.data[0]); yup = dref(y.data[1]);
            xlow = dref(x.data[0]); xup = dref(x.data[1]);
        }

        r.push_back(ylow + (yup - ylow)*(tx - xlow)/(xup - xlow));
    }

    return r;
}

// Perform linear interpolation of data 'y' of position 'x' at new position 'nx'.
// Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
// the arrays contains special values (NaN, inf, ...), all the points that would use these values
// will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
template<typename TypeY = double, typename TypeX = double, typename T = double,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
auto interpolate(const vec_t<1,TypeY>& y, const vec_t<1,TypeX>& x, const T& nx) {
    using rtypey = rtype_t<TypeY>;
    using rtypex = rtype_t<TypeX>;

    phypp_check(n_elements(y) == n_elements(x),
        "interpolate: 'x' and 'y' arrays must contain the same number of elements");
    phypp_check(n_elements(y) >= 2,
        "interpolate: 'x' and 'y' arrays must contain at least 2 elements");

    uint_t nmax = x.size();
    uint_t low = lower_bound(nx, x);

    rtypey ylow, yup;
    rtypex xlow, xup;
    if (low != npos) {
        if (low != nmax-1) {
            ylow = dref(y.data[low]); yup = dref(y.data[low+1]);
            xlow = dref(x.data[low]); xup = dref(x.data[low+1]);
        } else {
            ylow = dref(y.data[low-1]); yup = dref(y.data[low]);
            xlow = dref(x.data[low-1]); xup = dref(x.data[low]);
        }
    } else {
        ylow = dref(y.data[0]); yup = dref(y.data[1]);
        xlow = dref(x.data[0]); xup = dref(x.data[1]);
    }

    return ylow + (yup - ylow)*(nx - xlow)/(xup - xlow);
}

template<typename TypeX, typename TypeY>
auto integrate(const vec_t<1,TypeX>& x, const vec_t<1,TypeY>& y) -> decltype(0.5*y[0]*x[0]) {
    phypp_check(n_elements(x) == n_elements(y),
        "integrate: incompatible x and y array dimensions ("+strn(x.size())+" vs "+
        strn(y.size())+")");

    decltype(0.5*y[0]*x[0]) r = 0;
    for (uint_t i = 0; i < x.size()-1; ++i) {
        r += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i]);
    }

    return r;
}

template<typename F, typename T, typename U>
auto integrate_trap(F f, T x0, U x1, uint_t n) -> decltype(0.5*f(x0)*x0) {
    using rtype = decltype(0.5*f(x0)*x0);
    rtype r = 0;
    rtype y0 = f(x0);
    T x = x0;
    auto dx = (x1 - x0)/double(n);
    for (uint_t i = 0; i < n; ++i) {
        x += dx;
        rtype y1 = f(x);
        r += 0.5*(y1+y0)*dx;
        y0 = y1;
    }

    return r;
}

template<typename F, typename T, typename U>
auto integrate(F f, T x0, U x1, double e = std::numeric_limits<decltype(f(x0))>::epsilon()) ->
    decltype(0.5*f(x0)*x0) {

    using rtype = decltype(0.5*f(x0)*x0);

    std::vector<rtype> buffer;
    buffer.reserve(20);
    buffer.push_back(0.5*(x1 - x0)*(f(x0) + f(x1)));

    uint_t n = 0;
    uint_t oid = 0;

    do {
        ++n;

        rtype tr = 0;
        uint_t tn = 1 << n;
        auto d = (x1 - x0)/double(tn);
        for (uint_t k = 1; k <= tn/2; ++k) {
            tr += f(x0 + (2*k-1)*d);
        }
        buffer.push_back(0.5*buffer[oid] + d*tr);

        for (uint_t m = 1; m <= n; ++m) {
            auto t = 1 << (2*m);
            buffer.push_back((t*buffer.back() - buffer[oid+m-1])/(t - 1));
        }

        oid += n;
    } while (fabs((buffer.back() - *(buffer.end()-2))/buffer.back()) > e);

    return buffer.back();
}

// Perform the convolution of two arrays, assuming that they are based on the same 'x' coordinate
template<typename TypeX, typename TypeY1, typename TypeY2>
auto convolve(const vec_t<1,TypeX>& x, const vec_t<1,TypeY1>& y1, const vec_t<1,TypeY2>& y2) ->
    vec_t<1,decltype(x[0]*y1[0]*y2[0])> {

    phypp_check(x.size() > 3, "convolve needs arrays of at least 3 elements to work");

    vec_t<1,decltype(x[0]*y1[0]*y2[0])> r(x.size());
    for (uint_t i = 0; i < x.size(); ++i) {
        auto dx = x[i];
        if (i == 0) dx = x[i+1] - dx;
        else dx -= x[i-1];

        auto tmp = interpolate(y2, x + x[i], x);
        r += y1[i]*tmp*dx;
    }

    return r;
}

// Build the convex hull of a set of points, returning the indices of the points that form the hull
// in counter-clockwise order.
// Uses the monotone chain algorithm, taken from:
// http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C.2B.2B
template<typename TX, typename TY>
vec1u convex_hull(const TX& x, const TY& y) {
    phypp_check(n_elements(x) == n_elements(y),
        "convex_hull requires same number of elements in 'x' and 'y'");

    uint_t npt = n_elements(x);
    vec1u res = uintarr(2*npt);
    vec1u ids = uindgen(npt);
    std::sort(ids.data.begin(), ids.data.end(), [&](uint_t i, uint_t j) {
        if (x[i] < x[j]) return true;
        if (x[i] > x[j]) return false;
        if (y[i] < y[j]) return true;
        return false;
    });

    auto cross = [&](uint_t i, uint_t j, uint_t k) {
        return (x[j] - x[i])*(y[k] - y[i]) - (y[j] - y[i])*(x[k] - x[i]);
    };

    uint_t k = 0;
    for (uint_t i = 0; i < npt; ++i) {
        while (k >= 2 && cross(res[k-2], res[k-1], ids[i]) <= 0) --k;
        res[k] = ids[i];
        ++k;
    }

    uint_t t = k+1;
    for (uint_t i = npt-2; i != npos; --i) {
        while (k >= t && cross(res[k-2], res[k-1], ids[i]) <= 0) --k;
        res[k] = ids[i];
        ++k;
    }

    res.data.resize(k);
    res.dims[0] = k;

    return res;
}

template<typename TX, typename TY, typename THX, typename THY>
bool in_convex_hull(const TX& x, const TY& y, const vec1u& hull, const THX& hx, const THY& hy) {
    phypp_check(n_elements(hx) == n_elements(hy),
        "is_in_convex_hull requires same number of elements in 'hx' and 'hy'");

    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull[i], p2 = hull[i+1];
        auto cross = (hx[p2] - hx[p1])*(y - hy[p1]) - (hy[p2] - hy[p1])*(x - hx[p1]);
        if (cross < 0) return false;
    }

    return true;
}

template<std::size_t Dim, typename TX, typename TY, typename THX, typename THY>
vec_t<Dim,bool> in_convex_hull(const vec_t<Dim,TX>& x, const vec_t<Dim,TY>& y, const vec1u& hull,
    const THX& hx, const THY& hy) {

    phypp_check(n_elements(x) == n_elements(y),
        "is_in_convex_hull requires same number of elements in 'x' and 'y'");
    phypp_check(n_elements(hx) == n_elements(hy),
        "is_in_convex_hull requires same number of elements in 'hx' and 'hy'");

    vec_t<Dim,bool> res = replicate(true, x.dims);
    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull[i], p2 = hull[i+1];
        for (uint_t p = 0; p < x.size(); ++p) {
            if (!res[p]) continue;
            auto cross = (hx[p2] - hx[p1])*(y[p] - hy[p1]) - (hy[p2] - hy[p1])*(x[p] - hx[p1]);
            if (cross < 0) res[p] = false;
        }
    }

    return res;
}

#endif
