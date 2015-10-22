#ifndef MATH_HPP
#define MATH_HPP

#include "phypp/vec.hpp"
#include "phypp/string.hpp"
#include "phypp/simple_complex.hpp"
#include <random>
#include <cmath>

#ifndef NO_GSL
#include <gsl/gsl_sf_bessel.h>
#endif

#ifndef NO_FFTW
#include <fftw3.h>
#endif

static const double dnan = std::numeric_limits<double>::quiet_NaN();
static const float  fnan = std::numeric_limits<float>::quiet_NaN();
static const double dinf = std::numeric_limits<double>::infinity();
static const float  finf = std::numeric_limits<float>::infinity();
static const double dpi = 3.14159265359;
static const float  fpi = 3.14159265359;

template<typename T>
auto e10(const T& t) -> decltype(pow(10.0, t)) {
    return pow(10.0, t);
}

template<typename T, typename U, typename V,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T clamp(const T& t, const U& mi, const V& ma) {
    return (t < mi ? mi : (t > ma ? ma : t));
}

// Create a range with integer step from 0 to n (exclusive)
template<typename T>
vec1u rgen(T n) {
    phypp_check(n >= 0, "'rgen(n)' needs a positive or null value for 'n' (got ", n, ")");

    vec1u v(n);
    for (uint_t k : range(uint_t(n))) {
        v.safe[k] = k;
    }

    return v;
}

// Create a range with integer step from i to j (inclusive)
template<typename T, typename U>
vec<1,T> rgen(T i, U j) {
    if (i < T(j)) {
        uint_t n = j-i+1;
        vec<1,T> v(n);
        for (uint_t k : range(n)) {
            v.safe[k] = i+k;
        }
        return v;
    } else {
        uint_t n = i-j+1;
        vec<1,T> v(n);
        for (uint_t k : range(n)) {
            v.safe[k] = i-k;
        }
        return v;
    }
}

// Create a range of n steps from i to j (inclusive)
template<typename T, typename U, typename V>
vec1d rgen(T i, U j, V n) {
    phypp_check(n >= 0, "'rgen(a,b,n)' needs a positive or null value for 'n' (got ", n, ")");

    if (n == 1) {
        return {i};
    } else {
        vec1d v(n);
        double dx = (j-i)/double(n-1);
        for (uint_t k : range(uint_t(n))) {
            v.safe[k] = i + k*dx;
        }

        return v;
    }
}

template<typename T, typename U, typename V>
vec1d rgen_log(T i, U j, V n) {
    phypp_check(n >= 0, "'rgen_log(a,b,n)' needs a positive or null value for 'n' (got ", n, ")");

    if (n == 1) {
        return {i};
    } else {
        vec1d v(n);
        double dx = log10(j/i)/double(n-1);
        for (uint_t k : range(uint_t(n))) {
            v.safe[k] = i*e10(k*dx);
        }

        return v;
    }
}

template<typename T>
vec<2,T> make_bins(T mi, T ma) {
    return {{mi}, {ma}};
}

template<typename T>
vec2d make_bins(T mi, T ma, uint_t n) {
    vec2d b(2,n);
    double d = (ma - mi)/double(n);
    for (uint_t i : range(n)) {
        b.safe(0,i) = mi + i*d;
        b.safe(1,i) = mi + (i+1)*d;
    }

    return b;
}

template<typename T = double>
vec<2,T> make_bins(const vec<1,T>& v) {
    vec<2,T> b(2, v.size()-1);
    for (uint_t i : range(b.dims[1])) {
        b.safe(0,i) = v.safe[i];
        b.safe(1,i) = v.safe[i+1];
    }

    return b;
}

template<typename T, typename B = double>
bool in_bin(T t, const vec<1,B>& b) {
    phypp_check(b.dims[0] == 2, "can only be called on a single bin "
        "vector (expected dims=[2], got dims=[", b.dims, "])");

    return t >= b.safe[0] && t < b.safe[1];
}

template<std::size_t Dim, typename Type, typename B = double>
vec<Dim,bool> in_bin(const vec<Dim,Type>& v, const vec<1,B>& b) {
    phypp_check(b.dims[0] == 2, "can only be called on a single bin "
        "vector (expected dims=[2], got dims=[", b.dims, "])");

    auto low = b.safe[0];
    auto up  = b.safe[1];
    vec<Dim,bool> res(v.dims);
    for (uint_t i : range(v)) {
        res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
    }

    return res;
}

template<typename T, typename B = double>
bool in_bin(T t, const vec<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "can only be called on a single bin "
        "vector (expected dims=[2], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    return t >= b.safe(0,ib) && t < b.safe(1,ib);
}

template<std::size_t Dim, typename Type, typename B = double>
vec<Dim,bool> in_bin(const vec<Dim,Type>& v, const vec<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    auto low = b.safe(0,ib);
    auto up  = b.safe(1,ib);
    vec<Dim,bool> res(v.dims);
    for (uint_t i : range(v)) {
        res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
    }

    return res;
}

template<typename T, typename B = double>
bool in_bin_open(T t, const vec<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "can only be called on a single bin "
        "vector (expected dims=[2], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    if (ib == 0) {
        return t < b.safe(1,ib);
    } else if (ib == b.dims[1]-1) {
        return t >= b.safe(0,ib);
    } else {
        return t >= b.safe(0,ib) && t < b.safe(1,ib);
    }
}

template<std::size_t Dim, typename Type, typename B = double>
vec<Dim,bool> in_bin_open(const vec<Dim,Type>& v, const vec<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    auto low = b.safe(0,ib);
    auto up  = b.safe(1,ib);
    vec<Dim,bool> res(v.dims);
    if (ib == 0) {
        for (uint_t i : range(v)) {
            res.safe[i] = v.safe[i] < up;
        }
    } else if (ib == b.dims[1]-1) {
        for (uint_t i : range(v)) {
            res.safe[i] = v.safe[i] >= low;
        }
    } else {
        for (uint_t i : range(v)) {
            res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
        }
    }

    return res;
}

template<typename Type>
auto bin_center(const vec<2,Type>& b) -> vec<1,decltype(0.5*b[0])> {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");

    return 0.5*(b.safe(1,_) + b.safe(0,_));
}

template<typename Type>
auto bin_center(const vec<1,Type>& b) -> decltype(0.5*b[0]) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=2, got dims=", b.dims, ")");

    return 0.5*(b.safe[1] + b.safe[0]);
}

template<typename Type>
vec<1,rtype_t<Type>> bin_width(const vec<2,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");

    return b.safe(1,_) - b.safe(0,_);
}

template<typename Type>
rtype_t<Type> bin_width(const vec<1,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=2, got dims=", b.dims, ")");

    return b.safe[1] - b.safe[0];
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool is_finite(const T& t) {
    return std::isfinite(t);
}

template<std::size_t Dim, typename Type>
vec<Dim,bool> is_finite(const vec<Dim,Type>& v) {
    vec<Dim,bool> r(v.dims);
    for (uint_t i : range(v)) {
        r.safe[i] = std::isfinite(v.safe[i]);
    }

    return r;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool is_nan(const T& t) {
    return std::isnan(t);
}

template<std::size_t Dim, typename Type>
vec<Dim,bool> is_nan(const vec<Dim,Type>& v) {
    vec<Dim,bool> r(v.dims);
    for (uint_t i : range(v)) {
        r.safe[i] = std::isnan(v.safe[i]);
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
vec<dim_total<Args...>::value,double> randomn(T& seed, Args&& ... args) {
    vec<dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed);
    }

    return v;
}

template<typename T>
double randomu(T& seed) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(seed);
}

template<typename T, typename ... Args>
vec<dim_total<Args...>::value,double> randomu(T& seed, Args&& ... args) {
    vec<dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed);
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
auto randomi(T& seed, TMi mi, TMa ma, Args&& ... args) ->
    vec<dim_total<Args...>::value,decltype(mi+ma)> {
    auto v = randomu(seed, std::forward<Args>(args)...);
    using rtype = decltype(mi + ma);
    return vec<vec_dim<decltype(v)>::value,rtype>(v*(ma + 1 - mi) + mi);
}

template<typename T, typename TypeX, typename TypeY, typename ... Args>
rtype_t<TypeX> random_pdf(T& seed, const vec<1,TypeX>& px, const vec<1,TypeY>& py) {
    // TODO: make an alternative version for integers using std::discrete_distribution.

    using rtype = rtype_t<TypeX>;
    std::piecewise_linear_distribution<rtype> distribution(px.begin(), px.end(), py.begin());
    return distribution(seed);
}

template<typename T, typename TypeX, typename TypeY, typename ... Args>
vec<dim_total<Args...>::value,rtype_t<TypeX>> random_pdf(T& seed, const vec<1,TypeX>& px,
    const vec<1,TypeY>& py, Args&& ... args) {

    // TODO: make an alternative version for integers using std::discrete_distribution.

    using rtype = rtype_t<TypeX>;
    vec<dim_total<Args...>::value,rtype> v(std::forward<Args>(args)...);
    std::piecewise_linear_distribution<rtype> distribution(px.begin(), px.end(), py.begin());
    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed);
    }

    return v;
}

template<typename TSeed>
bool random_coin(TSeed& seed, double prob) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(seed) <= prob;
}

template<typename TSeed, uint_t D, typename T>
vec<D,bool> random_coin(TSeed& seed, const vec<D,T>& prob) {
    vec<D,bool> v(prob.dims);

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed) <= prob.safe[i];
    }

    return v;
}

template<typename T, typename ... Args>
vec<dim_total<Args...>::value,bool> random_coin(T& seed, double prob, Args&& ... args) {
    vec<dim_total<Args...>::value,bool> v(std::forward<Args>(args)...);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed) <= prob;
    }

    return v;
}

template<std::size_t Dim, typename Type, typename T>
void inplace_shuffle(T& seed, vec<Dim,Type>& v) {
    std::shuffle(v.begin(), v.end(), seed);
}

template<std::size_t Dim, typename Type, typename T>
vec<Dim,Type> shuffle(T& seed, vec<Dim,Type> v) {
    inplace_shuffle(seed, v);
    return v;
}

template<typename F, F f, std::size_t Dim, typename Type, typename ... Args>
auto run_index_(const vec<Dim,Type>& v, const Args& ... args, uint_t dim) ->
vec<Dim-1,typename return_type<F>::type> {

    vec<Dim-1,typename return_type<F>::type> r;
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

    vec<1,rtype_t<Type>> tmp(v.dims[dim]);
    for (uint_t i : range(r)) {
        uint_t base = (i%mpitch) + (i/mpitch)*v.dims[dim]*mpitch;
        for (uint_t j : range(tmp)) {
            tmp.safe[j] = v.safe[base + j*mpitch];
        }

        r.safe[i] = (*f)(tmp, args...);
    }

    return r;
}

template<std::size_t Dim, typename Type>
vec<1,rtype_t<Type>> run_dim_apply_ids_(const vec1u& ids, const vec<Dim,Type>& v) {
    return v.safe[ids].concretise();
}

template<std::size_t Dim, typename Type, typename ... Args>
std::array<uint_t,Dim> run_dim_get_dim_(const vec<Dim,Type>& v, const Args& ... vs) {
    // TODO: add a check here to make sure that all the other dims are the same
    return v.dims;
}

template<typename F, typename ... Args>
void run_dim_final_(uint_t dim, F&& func, const Args& ... vs) {
    auto ds = run_dim_get_dim_(vs...);
    const uint_t N = array_size<decltype(ds)>::size;

    uint_t nint = ds[dim];
    uint_t mpitch = 1;
    uint_t np = 1;
    for (uint_t i = 0; i < N; ++i) {
        if (i != dim) np *= ds[i];
        if (i > dim) mpitch *= ds[i];
    }

    vec1u ids(nint);
    for (uint_t i : range(np)) {
        uint_t base = (i%mpitch) + (i/mpitch)*nint*mpitch;
        for (uint_t j : range(ids)) {
            ids.safe[j] = base + j*mpitch;
        }

        func(i, run_dim_apply_ids_<N>(ids, vs)...);
    }
}

template<typename ... Args1>
struct run_dim_unroll_ {
    template<typename T>
    static void run(uint_t dim, Args1&&... a1, placeholder_t, T&& t) {
        run_dim_final_(dim, std::forward<T>(t), std::forward<Args1>(a1)...);
    }

    template<typename T, typename ... Args2>
    static void run(uint_t dim, Args1&&... a1, placeholder_t, T&& t, Args2&&... a2) {
        run_dim_unroll_<Args1..., T>::run(dim, std::forward<Args1>(a1)...,
            std::forward<T>(t), _, std::forward<Args2>(a2)...);
    }
};

// Iterate over one dimension of the provided vector and call a function for each slice.
template<typename ... Args>
void run_dim(uint_t dim, Args&& ... args) {
    run_dim_unroll_<>::run(dim, _, std::forward<Args>(args)...);
}

template<typename F, std::size_t Dim, typename Type>
auto reduce(uint_t dim, const vec<Dim,Type>& v, F&& func) ->
    vec<Dim-1,typename return_type<F>::type> {

    vec<Dim-1,typename return_type<F>::type> r;
    for (uint_t i = 0; i < dim; ++i) {
        r.dims[i] = v.dims[i];
    }
    for (uint_t i = dim+1; i < Dim; ++i) {
        r.dims[i-1] = v.dims[i];
    }

    r.resize();

    run_dim(dim, v, [&](uint_t i, vec<1,rtype_t<Type>> tv) {
        r.safe[i] = func(std::move(tv));
    });

    return r;
}

template<typename T>
using total_return_type = typename std::conditional<std::is_integral<T>::value,
    typename std::conditional<std::is_unsigned<T>::value, uint_t, int_t>::type,
    double>::type;

template<std::size_t Dim, typename Type>
total_return_type<rtype_t<Type>> total(const vec<Dim,Type>& v) {
    total_return_type<rtype_t<Type>> total = 0;
    for (auto& t : v) {
        total += t;
    }

    return total;
}


template<std::size_t Dim, typename Type, typename enable =
    typename std::enable_if<std::is_same<rtype_t<Type>, bool>::value>::type>
uint_t count(const vec<Dim,Type>& v) {
    uint_t n = 0u;
    for (bool b : v) {
        if (b) ++n;
    }

    return n;
}

template<std::size_t Dim, typename Type>
double mean(const vec<Dim,Type>& v) {
    double total = 0.0;
    for (auto& t : v) {
        total += t;
    }

    return total/v.size();
}

template<std::size_t Dim, typename T, typename enable =
    typename std::enable_if<std::is_same<rtype_t<T>, bool>::value>::type>
double fraction_of(const vec<Dim,T>& b) {
    return mean(b);
}

template<std::size_t Dim, typename T, typename U, typename enable =
    typename std::enable_if<std::is_same<rtype_t<T>, bool>::value &&
                            std::is_same<rtype_t<U>, bool>::value>::type>
double fraction_of(const vec<Dim,T>& b, const vec<Dim,U>& among) {
    phypp_check(b.dims == among.dims, "incompatible dimensions between count vector and "
        "among vector (", b.dims, " vs. ", among.dims, ")");

    double total = 0.0;
    uint_t count = 0;
    for (uint_t i : range(among)) {
        if (among.safe[i]) {
            total += b.safe[i];
            ++count;
        }
    }

    return total/count;
}

template<typename ... Args>
std::string percent_of(Args&& ... args) {
    return strn(round(100.0*fraction_of(std::forward<Args>(args)...)))+"%";
}

template<std::size_t Dim, typename Type>
rtype_t<Type> inplace_median(vec<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the median of an empty vector");

    using vtype = typename vec<Dim,Type>::vtype;
    using dtype = typename vtype::value_type;

    uint_t nwrong = 0;
    for (auto& t : v) {
        nwrong += is_nan(t);
    }

    if (nwrong == v.size()) return dnan;

    std::ptrdiff_t offset = (v.size()-nwrong)/2;
    std::nth_element(v.data.begin(), v.data.begin() + offset, v.data.end(),
        [](dtype i, dtype j) {
            if (is_nan(dref<Type>(i))) return false;
            if (is_nan(dref<Type>(j))) return true;
            return dref<Type>(i) < dref<Type>(j);
        }
    );

    return *(v.begin() + offset);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> median(vec<Dim,Type> v) {
    return inplace_median(v);
}

template<std::size_t Dim, typename Type, typename U>
rtype_t<Type> percentile(const vec<Dim,Type>& v, const U& u) {
    phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

    vec1u ok = where(is_finite(v));
    if (ok.empty()) return 0;

    // TODO: use same algorithm than median
    typename vec<1,Type>::effective_type t = v.safe[ok];
    std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    return *(t.begin() + offset);
}

template<std::size_t Dim, typename Type>
void percentiles_(vec<1,Type>& r, uint_t i, vec<Dim,Type>& t) {}

template<std::size_t Dim, typename Type, typename U, typename ... Args>
void percentiles_(vec<1,Type>& r, uint_t i, vec<Dim,Type>& t, const U& u, const Args& ... args) {
    std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    r.safe[i] = *(t.begin() + offset);
    ++i;

    percentiles_(r, i, t, args...);
}

template<std::size_t Dim, typename Type, typename ... Args>
typename vec<1,Type>::effective_type percentiles(const vec<Dim,Type>& v, const Args& ... args) {
    phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

    vec1u ok = where(is_finite(v));
    typename vec<1,Type>::effective_type t;
    if (ok.empty()) return t;
    t = v.safe[ok];

    typename vec<1,Type>::effective_type r = arr<rtype_t<Type>>(sizeof...(Args));
    percentiles_(r, 0, t, args...);

    return r;
}

template<std::size_t Dim, typename Type>
vec<Dim,bool> sigma_clip(const vec<Dim,Type>& tv, double sigma) {
    auto v = tv.concretise();
    auto med = inplace_median(v);
    auto mad = median(fabs(v - med));
    // Note: cannot use 'v' below, since the order of the values has changed!
    return fabs(tv - med) < sigma*mad;
}

template<std::size_t Dim, typename Type>
typename vec<Dim,Type>::const_iterator min_(const vec<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the minimum of an empty vector");

    auto iter = std::min_element(v.begin(), v.end(), [](rtype_t<Type> t1, rtype_t<Type> t2){
        if (is_nan(t1)) return false;
        if (is_nan(t2)) return true;
        return t1 < t2;
    });

    if (iter == v.end()) iter = v.begin();
    return iter;
}

template<std::size_t Dim, typename Type>
typename vec<Dim,Type>::const_iterator max_(const vec<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the maximum of an empty vector");

    auto iter = std::max_element(v.begin(), v.end(), [](rtype_t<Type> t1, rtype_t<Type> t2){
        if (is_nan(t1)) return true;
        if (is_nan(t2)) return false;
        return t1 < t2;
    });

    if (iter == v.end()) iter = v.begin();
    return iter;
}

template<std::size_t Dim, typename Type>
rtype_t<Type> min(const vec<Dim,Type>& v) {
    return *min_(v);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> max(const vec<Dim,Type>& v) {
    return *max_(v);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> min(const vec<Dim,Type>& v, uint_t& id) {
    auto iter = min_(v);
    id = iter - v.begin();
    return *iter;
}

template<std::size_t Dim, typename Type>
rtype_t<Type> max(const vec<Dim,Type>& v, uint_t& id) {
    auto iter = max_(v);
    id = iter - v.begin();
    return *iter;
}

template<std::size_t Dim, typename Type>
uint_t min_id(const vec<Dim,Type>& v) {
    return min_(v) - v.begin();
}

template<std::size_t Dim, typename Type>
uint_t max_id(const vec<Dim,Type>& v) {
    return max_(v) - v.begin();
}

template<typename T1, typename T2>
auto min(T1&& t1, T2&& t2) -> decltype((t1 < t2 ? t1 : t2)) {
    return t1 < t2 ? t1 : t2;
}

template<typename T1, typename T2>
auto max(T1&& t1, T2&& t2) -> decltype((t1 > t2 ? t1 : t2)) {
    return t1 > t2 ? t1 : t2;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> min(const vec<Dim,Type1>& v1, const vec<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "min: incompatible vector dimensions "
        "(", v1.dims, " vs. ", v2.dims, ")");

    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::min<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> max(const vec<Dim,Type1>& v1, const vec<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "max: incompatible vector dimensions "
        "(", v1.dims, " vs. ", v2.dims, ")");

    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::max<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> min(const vec<Dim,Type1>& v1, const Type2& v2) {
    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::min<decltype(v1[0]*v2)>(v1.safe[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> max(const vec<Dim,Type1>& v1, const Type2& v2) {
    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::max<decltype(v1[0]*v2)>(v1.safe[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> min(const Type1& v1, const vec<Dim,Type2>& v2) {
    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r.safe[i] = std::min<decltype(v1*v2[0])>(v1, v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec<Dim,rtype_t<Type1>> max(const Type1& v1, const vec<Dim,Type2>& v2) {
    vec<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r.safe[i] = std::max<decltype(v1*v2[0])>(v1, v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type>
double rms(const vec<Dim,Type>& v) {
    double sum = 0;
    for (auto& t : v) {
        sum += t*t;
    }

    return sqrt(sum/v.size());
}

template<std::size_t Dim, typename Type>
double stddev(const vec<Dim,Type>& v) {
    return rms(v - mean(v));
}

template<std::size_t Dim, typename Type>
rtype_t<Type> mad(const vec<Dim,Type>& v) {
    return median(fabs(v - median(v)));
}

#define RUN_INDEX(func) \
    struct func ## _run_index_wrapper_ { \
        template<typename T, typename ... Args> \
        static auto run(const vec<1,T>& v, Args&& ... args) -> \
        decltype(func(v, std::forward<Args>(args)...)) { \
            return func(v, std::forward<Args>(args)...); \
        } \
    }; \
    \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto partial_ ## func (uint_t dim, const vec<Dim,Type>& v, Args&& ... args) -> \
    vec<Dim-1, decltype(func(std::declval<vec<1,rtype_t<Type>>>(), std::forward<Args>(args)...))> { \
        using fptr = decltype(func(std::declval<vec<1,rtype_t<Type>>>(), std::forward<Args>(args)...)) \
            (*)(const vec<1,rtype_t<Type>>&, Args&& ...); \
        using wrapper = func ## _run_index_wrapper_; \
        return run_index_<fptr, &wrapper::run<rtype_t<Type>, Args...>>(v, dim); \
    }

RUN_INDEX(total);
RUN_INDEX(count);
RUN_INDEX(fraction_of);
RUN_INDEX(mean);
RUN_INDEX(median);
RUN_INDEX(min);
RUN_INDEX(max);
RUN_INDEX(percentile);
RUN_INDEX(rms);
RUN_INDEX(stddev);
RUN_INDEX(mad);

#undef RUN_INDEX

template<std::size_t Dim, typename Type, typename TypeB>
vec1u histogram(const vec<Dim,Type>& data, const vec<2,TypeB>& bins) {
    phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2,...], got dims=[", bins.dims, "])");

    using rtype = rtype_t<Type>;
    vec<Dim,rtype> tmp = data;

    uint_t nbin = bins.dims[1];
    vec1u counts(nbin);

    auto first = tmp.data.begin();
    for (uint_t i : range(nbin)) {
        auto last = std::partition(first, tmp.data.end(), [&bins,i](rtype t) {
            return t >= bins.safe(0,i) && t < bins.safe(1,i);
        });

        counts.safe[i] = last - first;
        first = last;

        if (last == tmp.data.end()) break;
    }

    return counts;
}

template<std::size_t Dim, typename Type, typename TypeB, typename TypeW>
vec<1,rtype_t<TypeW>> histogram(const vec<Dim,Type>& data, const vec<Dim,TypeW>& weight,
    const vec<2,TypeB>& bins) {
    phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", bins.dims, "])");
    phypp_check(data.dims == weight.dims, "incompatible dimensions for data and weight "
        "(", data.dims, " vs. ", weight.dims, ")");

    vec1u tmp = uindgen(data.size());

    uint_t nbin = bins.dims[1];
    vec<1,rtype_t<TypeW>> counts(nbin);

    auto first = tmp.data.begin();
    for (uint_t i : range(nbin)) {
        auto last = std::partition(first, tmp.data.end(), [&bins,&data,i](uint_t id) {
            return data.safe[id] >= bins.safe(0,i) && data.safe[id] < bins.safe(1,i);
        });

        for (; first != last; ++first) {
            counts.safe[i] += weight.safe[*first];
        }

        if (last == tmp.data.end()) break;
    }

    return counts;
}

template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX,
    typename TypeBY, typename TypeF>
void histogram2d_impl(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
    const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins, TypeF&& func) {

    phypp_check(xbins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", xbins.dims, "])");
    phypp_check(ybins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", ybins.dims, "])");
    phypp_check(x.dims == y.dims, "incompatible dimensions for x and y (", x.dims, " vs. ",
        y.dims, ")");

    vec1u ids = uindgen(x.size());

    uint_t nxbin = xbins.dims[1];
    uint_t nybin = ybins.dims[1];

    auto firstx = ids.data.begin();
    for (uint_t i : range(nxbin)) {
        auto lastx = std::partition(firstx, ids.data.end(), [&xbins,&x,i](uint_t id) {
            return x.safe[id] >= xbins.safe(0,i) && x.safe[id] < xbins.safe(1,i);
        });

        auto firsty = firstx;
        for (uint_t j : range(nybin)) {
            auto lasty = std::partition(firsty, lastx, [&ybins,&y,j](uint_t id) {
                return y.safe[id] >= ybins.safe(0,j) && y.safe[id] < ybins.safe(1,j);
            });

            func(ids, i, j, firsty, lasty);

            if (lasty == lastx) break;
            firsty = lasty;
        }

        if (lastx == ids.data.end()) break;
        firstx = lastx;
    }
}

template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX, typename TypeBY>
vec2u histogram2d(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
    const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins) {

    vec2u counts(xbins.dims[1], ybins.dims[1]);

    using iterator = vec1u::iterator;
    histogram2d_impl(x, y, xbins, ybins,
        [&](const vec1u& ids, uint_t i, uint_t j, iterator i0, iterator i1) {
            counts.safe(i,j) = i1 - i0;
        }
    );

    return counts;
}

template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX,
    typename TypeBY, typename TypeF>
void histogram2d(const vec<Dim,TypeX>& x, const vec<Dim,TypeY>& y,
    const vec<2,TypeBX>& xbins, const vec<2,TypeBY>& ybins, TypeF&& func) {

    using iterator = vec1u::iterator;
    histogram2d_impl(x, y, xbins, ybins,
        [&](const vec1u& ids, uint_t i, uint_t j, iterator i0, iterator i1) {
            uint_t npts = i1 - i0;
            vec1u tids(npts);
            uint_t k0 = i0 - ids.data.begin();
            for (uint_t k : range(npts)) {
                tids.safe[k] = ids.safe[k0 + k];
            }

            func(i, j, std::move(tids));
        }
    );
}

template<std::size_t Dim, typename Type>
void data_info_(const vec<Dim,Type>& v) {
    vec1u idok = where(is_finite(v));
    print(idok.size(), "/", v.size(), " valid values (dims: ", v.dims, ")");
    if (idok.size() == 0) return;

    vec<1,rtype_t<Type>> tv = v.safe[idok];

    print(" min : ", min(tv));
    print(" 15% : ", percentile(tv, 0.15));
    print(" 50% : ", median(tv));
    print(" mean: ", mean(tv));
    print(" 85% : ", percentile(tv, 0.85));
    print(" max : ", max(tv));
    print(" rms : ", stddev(tv - median(tv)));
}

#define data_info(x) \
    print("data info: ", #x, " (", typeid(x).name(), ")"); \
    data_info_(x);

template<typename T>
auto sign(const T& t) -> decltype(2*(t >= 0) - 1) {
    return 2*(t >= 0) - 1;
}

template<std::size_t Dim, typename Type, typename U,
    typename enable = typename std::enable_if<!is_vec<U>::value>::type>
auto pow(const U& u, const vec<Dim,Type>& v) -> vec<Dim, decltype(pow(u,v(0)))> {
    vec<Dim, decltype(pow(u,v(0)))> r = v;
    for (auto& t : r) {
        t = pow(u, t);
    }
    return r;
}

template<std::size_t Dim, typename Type, typename U>
auto pow(const U& u, vec<Dim,Type>&& v) -> typename std::enable_if<!is_vec<U>::value &&
    !std::is_pointer<Type>::value && std::is_same<decltype(pow(u,v(0))), Type>::value,
    vec<Dim,Type>>::type {
    for (auto& t : v) {
        t = pow(u, t);
    }
    return std::move(v);
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
auto sqr(T t) -> decltype(t*t) {
    return t*t;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
auto invsqr(T t) -> decltype(1.0/(t*t)) {
    return 1.0/(t*t);
}

#define VECTORIZE(name) \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(const vec<Dim,Type>& v, const Args& ... args) -> \
        vec<Dim,decltype(name(v[0], args...))> { \
        using ntype = decltype(name(v[0], args...)); \
        vec<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(name(dref<Type>(t), args...)); \
        } \
        return r; \
    } \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(vec<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<Type>::value && std::is_same<decltype(name(v[0], args...)), Type>::value, \
        vec<Dim,Type>>::type { \
        for (auto& t : v) { \
            t = name(t, args...); \
        } \
        return std::move(v); \
    }

#define VECTORIZE2(name) \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(const vec<D,T1>& v1, const vec<D,T2>& v2, const Args& ... args) -> \
        vec<D,decltype(name(v1[0], v2[0], args...))> { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
            v1.dims, " vs. ", v2.dims, ")"); \
        using ntype = decltype(name(v1[0], v2[0], args...)); \
        vec<D,ntype> r; r.dims = v1.dims; r.data.reserve(v1.size()); \
        for (uint_t i : range(v1)) { \
            r.data.push_back(name(v1.safe[i], v2.safe[i], args...)); \
        } \
        return r; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(vec<D,T1>&& v1, const vec<D,T2>& v2, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<T1>::value && std::is_same<decltype(name(v1[0], v2[0], args...)), T1>::value, \
        vec<D,T1>>::type { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
            v1.dims, " vs. ", v2.dims, ")"); \
        for (uint_t i : range(v1)) { \
            v1.safe[i] = name(v1.safe[i], v2.safe[i], args...); \
        } \
        return v1; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(const vec<D,T1>& v1, vec<D,T2>&& v2, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<T2>::value && std::is_same<decltype(name(v1[0], v2[0], args...)), T2>::value, \
        vec<D,T2>>::type { \
        phypp_check(v1.dims == v2.dims, "incompatible dimensions between V1 and V2 (", \
            v1.dims, " vs. ", v2.dims, ")"); \
        for (uint_t i : range(v1)) { \
            v2.safe[i] = name(v1.safe[i], v2.safe[i], args...); \
        } \
        return v2; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(T1 v1, const vec<D,T2>& v2, const Args& ... args) -> typename std::enable_if<!is_vec<T1>::value, \
        vec<D,decltype(name(v1, v2[0], args...))>>::type { \
        using ntype = decltype(name(v1, v2[0], args...)); \
        vec<D,ntype> r; r.dims = v2.dims; r.data.reserve(v2.size()); \
        for (uint_t i : range(v2)) { \
            r.data.push_back(name(v1, v2.safe[i], args...)); \
        } \
        return r; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(const vec<D,T1>& v1, T2 v2, const Args& ... args) -> typename std::enable_if<!is_vec<T2>::value, \
        vec<D,decltype(name(v1[0], v2, args...))>>::type { \
        using ntype = decltype(name(v1[0], v2, args...)); \
        vec<D,ntype> r; r.dims = v1.dims; r.data.reserve(v1.size()); \
        for (uint_t i : range(v1)) { \
            r.data.push_back(name(v1.safe[i], v2, args...)); \
        } \
        return r; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(T1 v1, vec<D,T2>&& v2, const Args& ... args) -> typename std::enable_if<!is_vec<T1>::value && \
        !std::is_pointer<T2>::value && std::is_same<decltype(name(v1, v2[0], args...)), T2>::value, \
        vec<D,decltype(name(v1, v2[0], args...))>>::type { \
        for (uint_t i : range(v2)) { \
            v2.safe[i] = name(v1, v2.safe[i], args...); \
        } \
        return v2; \
    } \
    template<std::size_t D, typename T1, typename T2, typename ... Args> \
    auto name(vec<D,T1>&& v1, T2 v2, const Args& ... args) -> typename std::enable_if<!is_vec<T2>::value && \
        !std::is_pointer<T1>::value && std::is_same<decltype(name(v1[0], v2, args...)), T1>::value, \
        vec<D,decltype(name(v1[0], v2, args...))>>::type { \
        for (uint_t i : range(v1)) { \
            v1.safe[i] = name(v1.safe[i], v2, args...); \
        } \
        return v1; \
    } \

#define VECTORIZE_REN(name, orig) \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(const vec<Dim,Type>& v, const Args& ... args) -> \
        vec<Dim,decltype(orig(v[0], args...))> { \
        using ntype = decltype(orig(v[0], args...)); \
        vec<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(orig(dref<Type>(t), args...)); \
        } \
        return r; \
    } \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto name(vec<Dim,Type>&& v, const Args& ... args) -> typename std::enable_if< \
        !std::is_pointer<Type>::value && std::is_same<decltype(orig(v[0], args...)), Type>::value, \
        vec<Dim,Type>>::type { \
        for (auto& t : v) { \
            t = orig(t, args...); \
        } \
        return std::move(v); \
    } \
    template<typename ... Args> \
    auto name(Args&& ... args) -> decltype(orig(std::forward<Args>(args)...)) { \
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
VECTORIZE2(atan2);
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
auto derivate1_func(F&& func, double x, double ep) -> decltype(0.5*func(x)) {
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
auto derivate2_func(F&& func, double x, double ep) -> decltype(0.5*func(x)) {
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
auto derivate1_func(F&& func, const vec1d& x, double ep, uint_t ip) -> decltype(0.5*func(x)) {
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
auto derivate2_func(F&& func, const vec1d& x, double ep, uint_t ip) -> decltype(0.5*func(x)) {
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
auto derivate2_func(F&& func, const vec1d& x, double ep, uint_t ip1, uint_t ip2) -> decltype(0.5*func(x)) {
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
auto mmul(const vec<2,TypeA>& a, const vec<2,TypeB>& b) -> vec<2,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying matrices "
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];
    const uint_t m = b.dims[1];

    vec<2,ntype_t> r(n,m);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t j = 0; j < m; ++j)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i,j) += a.safe(i,k)*b.safe(k,j);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec<2,TypeA>& a, const vec<1,TypeB>& b) -> vec<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying matrix by vector "
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];

    vec<1,ntype_t> r(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i) += a.safe(i,k)*b.safe(k);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec<1,TypeB>& b, const vec<2,TypeA>& a) -> vec<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying vector by matrix "
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[0];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[1];

    vec<1,ntype_t> r = arr<ntype_t>(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i) += b.safe(k)*a.safe(k,i);
    }

    return r;
}

template<typename Type>
void mprint(const vec<2,Type>& m) {
    for (uint_t i = 0; i < m.dims[0]; ++i) {
        for (uint_t j = 0; j < m.dims[1]; ++j) {
            if (j != 0) std::cout << ", ";
            std::cout << m.safe(i,j);
        }

        std::cout << "\n";
    }

    std::cout << std::flush;
}

template<typename Type>
vec<2,rtype_t<Type>> transpose(const vec<2,Type>& v) {
    vec<2,rtype_t<Type>> r(v.dims);
    std::swap(r.dims[0], r.dims[1]);

    for (uint_t i : range(r.size())) {
        r[i] = v.safe[(i%r.dims[1])*v.dims[1] + i/r.dims[1]];
    }

    // TODO: see who's faster
    // for (uint_t i : range(r.dims[0]))
    // for (uint_t j : range(r.dims[1])) {
    //     r.safe[j+i*r.dims[1]] = v.safe[i+j*v.dims[1]];
    //     // r.safe(i,j) = v.safe(j,i);
    // }

    return r;
}

template<typename Type>
auto diagonal(vec<2,Type>&& v) -> decltype(v(_,0).concretise()) {
    phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
        v.dims, ")");

    decltype(v(_,0).concretise()) d(v.dims[0]);
    for (uint_t i : range(d)) {
        d.safe[i] = v.safe(i,i);
    }

    return d;
}

template<typename Type>
auto diagonal(const vec<2,Type>& v) -> decltype(v(_,0)) {
    phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
        v.dims, ")");

    decltype(v(_,0)) d(vec_ref_tag, vec_access::get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i : range(d)) {
        d.safe[i] = ptr<Type>(v.safe(i,i));
    }

    return d;
}

template<typename Type>
auto diagonal(vec<2,Type>& v) -> decltype(v(_,0)) {
    phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
        v.dims, ")");

    decltype(v(_,0)) d(vec_ref_tag, vec_access::get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i : range(d)) {
        d.data[i] = ptr<Type>(v.safe(i,i));
    }

    return d;
}

template<typename Type = double>
vec<2,Type> identity_matrix(uint_t dim) {
    vec<2,Type> m(dim, dim);
    diagonal(m) = 1;
    return m;
}

template<typename TX, typename TY>
auto scale_matrix(const TX& sx, const TY& sy) -> vec<2,decltype(sx*sy)> {
    vec<2,decltype(sx*sy)> m(3, 3);
    m.safe(0,0) = sx;
    m.safe(1,1) = sy;
    m.safe(2,2) = 1;
    return m;
}

template<typename T>
vec<2,T> scale_matrix(const T& s) {
    vec<2,T> m(3, 3);
    m.safe(0,0) = s;
    m.safe(1,1) = s;
    m.safe(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto translation_matrix(const TX& tx, const TY& ty) -> vec<2,decltype(tx*ty)> {
    vec<2,decltype(tx*ty)> m(3, 3);
    diagonal(m) = 1;
    m.safe(0,2) = tx;
    m.safe(1,2) = ty;
    return m;
}

template<typename A>
auto rotation_matrix(const A& a) -> vec<2,decltype(cos(a))> {
    vec<2,decltype(cos(a))> m(3, 3);
    auto ca = cos(a), sa = sin(a);
    m.safe(0,0) = m.safe(1,1) = ca;
    m.safe(0,1) = -sa;
    m.safe(1,0) = sa;
    m.safe(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto point2d(const TX& x, const TY& y) -> vec<1,decltype(x*y)> {
    return vec<1,decltype(x*y)>{x, y, 1};
}

// LAPACK functions
// ----------------

#ifndef NO_LAPACK
extern "C" void dgetrf_(int* n, int* m, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info);
extern "C" void dsytrf_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work,
    int* lwork, int* info);
extern "C" void dsytri_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work, int* info);
extern "C" void dsysv_(char* uplo, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,
    int* ldb, double* work, int* lwork, int* info);
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
    double* work, int* lwork, int* info);
#endif

template<typename Dummy = void>
bool inplace_invert(vec2d& i) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
    phypp_check(i.dims[0] == i.dims[1], "cannot invert a non square matrix (", i.dims, ")");

    int n = i.dims[0];
    int lda = n;
    int info;

    vec<1,int> ipiv(n);
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
#endif
}

template<typename TypeA>
bool invert(const vec<2,TypeA>& a, vec2d& i) {
    i = a;
    return inplace_invert<TypeA>(i);
}

template<typename Dummy = void>
bool inplace_invert_symmetric(vec2d& i) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
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
    vec<1,int> ipiv(n);

    dsytrf_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &lw, &info);
    if (info < 0) {
        return false;
    }

    dsytri_(&uplo, &n, i.data.data(), &lda, ipiv.data.data(), work.data.data(), &info);
    if (info != 0) {
        return false;
    }

    return true;
#endif
}

template<typename TypeA>
bool invert_symmetric(const vec<2,TypeA>& a, vec2d& i) {
    i = a;
    return inplace_invert_symmetric<TypeA>(i);
}

template<typename Dummy = void>
bool inplace_solve_symmetric(vec2d& alpha, vec1d& beta) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
    phypp_check(alpha.dims[0] == alpha.dims[1], "cannot invert a non square matrix (",
        alpha.dims, ")");
    phypp_check(alpha.dims[0] == beta.dims[0], "matrix and vector must have the same dimensions (",
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
    vec<1,int> ipiv(n);

    dsysv_(&uplo, &n, &nrhs, alpha.data.data(), &lda, ipiv.data.data(), beta.data.data(),
        &ldb, work.data.data(), &lw, &info);
    if (info != 0) {
        return false;
    }

    return true;
#endif
}

template<typename Dummy = void>
bool solve_symmetric(const vec2d& alpha, const vec1d& beta, vec1d& res) {
    vec2d a = alpha;
    res = beta;
    return inplace_solve_symmetric<Dummy>(a, res);
}

template<typename Dummy = void>
bool inplace_eigen_symmetric(vec2d& a, vec1d& vals) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
    phypp_check(a.dims[0] == a.dims[1], "cannot invert a non square matrix (",
        a.dims, ")");

    char jobz = 'V';
    char uplo = 'U';
    int n = a.dims[0];
    int lda = n;
    int info;

    vals.resize(n);

    int lw = n*64;
    // Note: the optimal value for lw is n*nb, where nb is the optimal block size
    // This value can be obtained using ilaenv_, but 64 should be plenty enough, according to
    // the Lapack User Guide.

    vec1d work(lw);

    dsyev_(&jobz, &uplo, &n, a.data.data(), &lda, vals.data.data(), work.data.data(),
        &lw, &info);
    if (info != 0) {
        return false;
    }

    // Eigen vectors are now stored in 'a' with the following layout:
    // v0 = a(0,_), v1 = a(1,_), ...
    // each corresponding to the eigen values given in 'vals'

    return true;
#endif
}

template<typename Dummy = void>
bool eigen_symmetric(const vec2d& a, vec1d& vals, vec2d& vecs) {
    vecs = a;
    return inplace_eigen_symmetric<Dummy>(vecs, vals);
}


// -----------------------
// end of LAPACK functions


template<typename Type>
void symmetrize(vec<2,Type>& alpha) {
    phypp_check(alpha.dims[0] == alpha.dims[1], "cannot symmetrize a non square matrix (",
        alpha.dims, ")");

    for (uint_t i = 0; i < alpha.dims[0]; ++i)
    for (uint_t j = i+1; j < alpha.dims[0]; ++j) {
        alpha.safe(i,j) = alpha.safe(j,i);
    }
}

struct linfit_result {
    bool success;

    double chi2;
    vec1d  params;
    vec1d  errors;
    vec2d  cov;

    // Reflection data
    MEMBERS1(success, chi2, params, errors, cov);
    MEMBERS2("linfit_result", MAKE_MEMBER(success), MAKE_MEMBER(chi2),
        MAKE_MEMBER(params), MAKE_MEMBER(errors), MAKE_MEMBER(cov));
};

template<typename T, typename TypeE>
void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t) {
    cache.safe(i,_) = flatten(t/ye);
}

template<typename T, typename TypeE, typename ... Args>
void linfit_make_cache_(vec2d& cache, const TypeE& ye, uint_t i, T&& t, Args&& ... args) {
    cache.safe(i,_) = flatten(t/ye);
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
                alpha.safe(i,j) = 0.0;
                // alpha(i,j) = sum over all points of x[i]*x[j]/e^2
                for (uint_t m = 0; m < nm; ++m) {
                    alpha.safe(i,j) += cache.safe(i,m)*cache.safe(j,m);
                }
            } else {
                alpha.safe(i,j) = alpha.safe(j,i);
            }
        }

        beta.safe[i] = 0.0;
        // beta[i] = sum over all points of x[i]*y/e^2
        for (uint_t m = 0; m < nm; ++m) {
            beta.safe[i] += cache.safe(i,m)*tmp.safe[m];
        }
    }

    if (!inplace_invert_symmetric(alpha)) {
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
    fr.errors = sqrt(diagonal(alpha));

    vec1d model(nm);
    for (uint_t m = 0; m < nm; ++m) {
        model.safe[m] = 0.0;
        for (uint_t i = 0; i < np; ++i) {
            model.safe[m] += fr.params.safe[i]*cache.safe(i,m);
        }
    }

    fr.chi2 = total(sqr(model - tmp));

    return fr;
}

template<typename TY>
void linfit_error_dims_(const TY& y, uint_t i) {}

template<typename TY, typename U, typename ... Args>
void linfit_error_dims_(const TY& y, uint_t i, const U& t, const Args& ... args) {
    phypp_check(same_dims_or_scalar(y, t), "incompatible dimensions between Y and X",
        i, " (", dim(y), " vs. ", dim(t), ")");
    linfit_error_dims_(y, i+1, args...);
}

template<typename TypeY, typename TypeE, typename ... Args>
linfit_result linfit(const TypeY& y, const TypeE& ye, Args&&... args) {
    bool bad = !same_dims_or_scalar(y, ye, args...);
    if (bad) {
        phypp_check(same_dims_or_scalar(y, ye), "incompatible dimensions between Y and "
            "YE arrays (", dim(y), " vs. ", dim(ye), ")");
        linfit_error_dims_(y, 0, args...);
    }

    uint_t np = sizeof...(Args);
    uint_t nm = n_elements(y);

    vec2d cache(np,nm);
    linfit_make_cache_(cache, ye, 0, std::forward<Args>(args)...);

    return linfit_do_(y, ye, cache);
}

template<std::size_t Dim, typename TypeY, typename TypeE, typename TypeX>
linfit_result linfit_pack(const vec<Dim,TypeY>& y, const vec<Dim,TypeE>& ye,
    const vec<Dim+1,TypeX>& x) {
    bool good = true;
    for (uint_t i : range(Dim)) {
        if (x.dims[i+1] != ye.dims[i] || x.dims[i+1] != y.dims[i]) {
            good = false;
            break;
        }
    }

    phypp_check(good, "incompatible dimensions between X, Y and YE arrays (", x.dims,
        " vs. ", y.dims, " vs. ", ye.dims, ")")

    linfit_result fr;

    uint_t np = x.dims[0];
    uint_t nm = y.size();

    vec2d cache(np,nm);
    for (uint_t i = 0; i < np; ++i) {
        for (uint_t j = 0; j < nm; ++j) {
            cache.safe(i,j) = x.safe[i*x.pitch(0) + j]/ye.safe[j];
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
affinefit_result affinefit(const TypeY& y, const vec<Dim,TypeE>& ye, const TypeX& x) {
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

template<typename TypeX, typename TypeY, typename TypeE>
affinefit_result affinefit(const TypeY& y, const TypeE& ye, const TypeX& x) {
    return affinefit(y, replicate(ye, y.dims), x);
}

double interpolate(double y1, double y2, double x1, double x2, double x) {
    return y1 + (y2 - y1)*(x - x1)/(x2 - x1);
}

// Perform linear interpolation of data 'y' of position 'x' at new positions 'nx'.
// Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
// the arrays contains special values (NaN, inf, ...), all the points that would use these values
// will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
template<std::size_t DI = 1, std::size_t DX = 1, typename TypeX2 = double, typename TypeY = double,
    typename TypeX1 = double>
auto interpolate(const vec<DI,TypeY>& y, const vec<DI,TypeX1>& x, const vec<DX,TypeX2>& nx) ->
    vec<DX,decltype(y[0]*x[0])> {

    using rtypey = rtype_t<TypeY>;
    using rtypex = rtype_t<TypeX1>;

    phypp_check(y.size() == x.size(),
        "interpolate: 'x' and 'y' arrays must contain the same number of elements");
    phypp_check(y.size() >= 2,
        "interpolate: 'x' and 'y' arrays must contain at least 2 elements");

    uint_t nmax = x.size();
    vec<DX,decltype(y[0]*x[0])> r; r.reserve(nx.size());
    for (auto& tx : nx) {
        uint_t low = lower_bound(tx, x);

        rtypey ylow, yup;
        rtypex xlow, xup;
        if (low != npos) {
            if (low != nmax-1) {
                ylow = y.safe[low]; yup = y.safe[low+1];
                xlow = x.safe[low]; xup = x.safe[low+1];
            } else {
                ylow = y.safe[low-1]; yup = y.safe[low];
                xlow = x.safe[low-1]; xup = x.safe[low];
            }
        } else {
            ylow = y.safe[0]; yup = y.safe[1];
            xlow = x.safe[0]; xup = x.safe[1];
        }

        r.push_back(ylow + (yup - ylow)*(tx - xlow)/(xup - xlow));
    }

    return r;
}

// Perform linear interpolation of data 'y' of position 'x' at new position 'nx'.
// Assumes that the arrays only contain finite elements, and that 'x' is properly sorted. If one of
// the arrays contains special values (NaN, inf, ...), all the points that would use these values
// will be contaminated. If 'x' is not properly sorted, the result will simply be wrong.
template<std::size_t DI, typename TypeY = double, typename TypeX = double, typename T = double,
    typename enable = typename std::enable_if<!is_vec<T>::value>::type>
auto interpolate(const vec<DI,TypeY>& y, const vec<DI,TypeX>& x, const T& nx) ->
    decltype(y[0]*x[0]) {

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
            ylow = y.safe[low]; yup = y.safe[low+1];
            xlow = x.safe[low]; xup = x.safe[low+1];
        } else {
            ylow = y.safe[low-1]; yup = y.safe[low];
            xlow = x.safe[low-1]; xup = x.safe[low];
        }
    } else {
        ylow = y.safe[0]; yup = y.safe[1];
        xlow = x.safe[0]; xup = x.safe[1];
    }

    return ylow + (yup - ylow)*(nx - xlow)/(xup - xlow);
}

template<typename Type>
rtype_t<Type> bilinear(const vec<2,Type>& map, double x, double y) {
    int_t tix = floor(x);
    int_t tiy = floor(y);
    double dx = x - tix;
    double dy = y - tiy;

    if (tix < 0) {
        tix = 0;
        dx = x - tix;
    } else if (uint_t(tix) >= map.dims[0]-1) {
        tix = map.dims[0]-2;
        dx = x - tix;
    }

    if (tiy < 0) {
        tiy = 0;
        dy = y - tiy;
    } else if (uint_t(tiy) >= map.dims[1]-1) {
        tiy = map.dims[1]-2;
        dy = y - tiy;
    }

    uint_t ix = tix;
    uint_t iy = tiy;

    return map.safe(ix,iy)*(1.0 - dx)*(1.0 - dy) + map.safe(ix,iy+1)*(1.0 - dx)*dy
        + map.safe(ix+1,iy)*dx*(1.0 - dy) + map.safe(ix+1,iy+1)*dx*dy;
}

template<typename Type, typename TypeD = double>
rtype_t<Type> bilinear_strict(const vec<2,Type>& map, double x, double y,
    TypeD def = 0.0) {

    int_t tix = floor(x);
    int_t tiy = floor(y);

    if (tix >= map.dims[0]-1 || tix < 0 || tiy >= map.dims[1]-1 || tiy < 0) {
        return def;
    }

    double dx = x - tix;
    double dy = y - tiy;
    uint_t ix = tix;
    uint_t iy = tiy;

    return map.safe(ix,iy)*(1.0 - dx)*(1.0 - dy) + map.safe(ix,iy+1)*(1.0 - dx)*dy
        + map.safe(ix+1,iy)*dx*(1.0 - dy) + map.safe(ix+1,iy+1)*dx*dy;
}

template<typename Type>
vec<2,rtype_t<Type>> rebin(const vec<2,Type>& map, const vec1d& mx,
    const vec1d& my, const vec1d& x, const vec1d& y) {

    phypp_check(map.dims[0] == mx.size(), "incompatible size of MAP and MX (", map.dims,
        " vs. ", mx.size(), ")");
    phypp_check(map.dims[1] == my.size(), "incompatible size of MAP and MY (", map.dims,
        " vs. ", my.size(), ")");

    vec<2,rtype_t<Type>> v(x.size(), y.size());

    vec1d ux = interpolate(dindgen(mx.size()), mx, x);
    vec1d uy = interpolate(dindgen(my.size()), my, y);

    for (uint_t ix : range(ux))
    for (uint_t iy : range(uy)) {
        v.safe(ix,iy) = bilinear(map, ux.safe[ix], uy.safe[iy]);
    }

    return v;
}

// TODO: move this function back into gencat.cpp, it is too specific and will
// likely not be used for anything else
template<typename TX1, typename TX2, typename TY1, typename TY2, typename TX, typename TY>
void merge_add(const vec<1,TX1>& x1, const vec<1,TX2>& x2,
    const vec<1,TY1>& y1, const vec<1,TY2>& y2,
    vec<1,TX>& x, vec<1,TY>& y) {

    phypp_check(x1.dims == y1.dims, "incompatible dimensions between X1 and Y1 (",
        x1.dims, " vs. ", y1.dims, ")");
    phypp_check(x2.dims == y2.dims, "incompatible dimensions between X2 and Y2 (",
        x2.dims, " vs. ", y2.dims, ")");

    uint_t n1 = x1.size(), n2 = x2.size();
    x.clear(); x.reserve(n1+n2);
    y.clear(); y.reserve(n1+n2);

    uint_t i1 = 0, i2 = 0;
    while (i1 < n1 || i2 < n2) {
        if (i1 == n1) {
            x.push_back(x2.safe[i2]);
            y.push_back(y2.safe[i2]);
            ++i2;
        } else if (i2 == n2) {
            x.push_back(x1.safe[i1]);
            y.push_back(y1.safe[i1]);
            ++i1;
        } else {
            if (x1.safe[i1] < x2.safe[i2]) {
                x.push_back(x1.safe[i1]);

                if (i2 == 0) {
                    y.push_back(y1.safe[i1]);
                } else {
                    y.push_back(y1.safe[i1] + interpolate(
                        y2.safe[i2-1], y2.safe[i2], x2.safe[i2-1], x2.safe[i2], x1.safe[i1]
                    ));
                }

                ++i1;
            } else {
                x.push_back(x2.safe[i2]);

                if (i1 == 0) {
                    y.push_back(y2.safe[i2]);
                } else {
                    y.push_back(y2.safe[i2] + interpolate(
                        y1.safe[i1-1], y1.safe[i1], x1.safe[i1-1], x1.safe[i1], x2.safe[i2]
                    ));
                }

                ++i2;
            }
        }
    }
}

// Find the first value in 'x' where 'y' equals zero, using linear interpolation
template<typename TypeX, typename TypeY>
double find_zero(const vec<1,TypeX>& x, const vec<1,TypeY>& y) {
    phypp_check(x.size() == y.size(),
        "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

    for (uint_t i : range(1, x.size())) {
        if (y.safe[i]*y.safe[i-1] <= 0) {
            return interpolate(x.safe[i-1], x.safe[i], y.safe[i-1], y.safe[i], 0.0);
        }
    }

    return dnan;
};

template<typename TypeX, typename TypeY>
auto integrate(const vec<1,TypeX>& x, const vec<1,TypeY>& y)
    -> decltype(0.5*y[0]*(x[1]-x[0])) {

    phypp_check(x.size() == y.size(),
        "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

    decltype(0.5*y[0]*(x[1]-x[0])) r = 0;
    for (uint_t i = 0; i < x.size()-1; ++i) {
        r += 0.5*(y.safe[i+1]+y.safe[i])*(x.safe[i+1]-x.safe[i]);
    }

    return r;
}

template<typename TypeX, typename TypeY>
auto integrate(const vec<2,TypeX>& x, const vec<1,TypeY>& y)
    -> decltype(y[0]*(x[1]-x[0])) {

    phypp_check(x.dims[0] == 2, "x array must be a binned array (expected dim=[2,...], "
        "got dim=[", x.dims[0], ",...]");
    phypp_check(x.dims[1] == y.size(),
        "incompatible x and y array dimensions (", x.dims[1], " vs ", y.size(), ")");

    decltype(y[0]*(x[1]-x[0])) r = 0;
    for (uint_t i = 0; i < x.dims[1]; ++i) {
        r += y.safe[i]*(x.safe[i+x.dims[1]]-x.safe[i]);
    }

    return r;
}

template<typename TypeX, typename TypeY>
auto integrate(const vec<1,TypeX>& x, const vec<1,TypeY>& y, double x0, double x1)
    -> decltype(0.5*y[0]*(x[1]-x[0])) {

    phypp_check(x.size() == y.size(),
        "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

    phypp_check(x.front() <= x0 && x.back() >= x1, "x array must cover the range [x0,x1]");

    uint_t i0 = upper_bound(x0, x);
    uint_t i1 = lower_bound(x1, x);

    if (i0 > i1) {
        return 0.5*(y.safe[i1]+y.safe[i0])*(x1 - x0);
    } else {
        decltype(0.5*y[0]*(x[1]-x[0])) r = 0;
        for (uint_t i = i0; i < i1; ++i) {
            r += 0.5*(y.safe[i+1]+y.safe[i])*(x.safe[i+1]-x.safe[i]);
        }

        if (i0 > 0) {
            double y0 = interpolate(y.safe[i0-1], y.safe[i0], x.safe[i0-1], x.safe[i0], x0);
            r += 0.5*(y.safe[i0]+y0)*(x[i0]-x0);
        }

        if (i1 < x.size()-1) {
            double y1 = interpolate(y.safe[i1], y.safe[i1+1], x.safe[i1], x.safe[i1+1], x1);
            r += 0.5*(y1+y.safe[i1])*(x1-x.safe[i1]);
        }

        return r;
    }
}

template<typename TypeX, typename TypeY>
auto integrate(const vec<2,TypeX>& x, const vec<1,TypeY>& y, double x0, double x1)
    -> decltype(y[0]*(x[1]-x[0])) {

    phypp_check(x.dims[0] == 2, "x array must be a binned array (expected dim=[2,...], "
        "got dim=[", x.dims[0], ",...]");
    phypp_check(x.dims[1] == y.size(),
        "incompatible x and y array dimensions (", x.dims[1], " vs ", y.size(), ")");

    phypp_check(x[0] <= x0 && x[x.size()-1] >= x1, "x array must cover the range [x0,x1]");

    uint_t i0 = upper_bound(x0, x(0,_));
    uint_t i1 = lower_bound(x1, x(0,_));

    if (i0 > i1) {
        return y.safe[i1]*(x1 - x0);
    } else {
        decltype(y[0]*(x[1]-x[0])) r = 0;
        for (uint_t i = i0; i < i1; ++i) {
            r += y.safe[i]*(x.safe[i+x.dims[1]]-x.safe[i]);
        }

        r += y.safe[i0]*(x.safe[i0] - x0);
        r += y.safe[i1]*(x1 - x.safe[i1]);

        return r;
    }
}

template<typename TypeX, typename TypeY>
auto cumul(const vec<1,TypeX>& x, const vec<1,TypeY>& y) -> vec<1,decltype(integrate(x,y))> {
    vec<1,decltype(integrate(x,y))> dr(y.dims);

    phypp_check(x.size() == y.size(),
        "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

    dr[0] = y[0];
    decltype(0.5*y[0]*(x[1]-x[0])) r = y[0];
    for (uint_t i : range(x.size()-1)) {
        r += 0.5*(y.safe[i+1]+y.safe[i])*(x.safe[i+1]-x.safe[i]);

        dr[i+1] = r;
    }

    return dr;
}

template<typename F, typename T, typename U,
    typename enable = typename std::enable_if<!is_vec<F>::value>::type>
auto integrate_func(F f, T x0, U x1,
    double e = std::numeric_limits<typename std::result_of<F(T)>::type>::epsilon()) ->
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

#ifndef NO_FFTW
// Compute the Fast Fourrier Transform (FFT) of the provided 2d array
vec2cd fft(const vec2d& v) {
    vec2cd r(v.dims);

    fftw_plan p;
    p = fftw_plan_dft_r2c_2d(v.dims[0], v.dims[1],
        const_cast<double*>(v.data.data()),
        reinterpret_cast<fftw_complex*>(r.data.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    return r;
}

// Compute the Fast Fourrier Transform (FFT) of the provided 2d array
vec2d ifft(const vec2cd& v) {
    vec2d r(v.dims);

    fftw_plan p;
    p = fftw_plan_dft_c2r_2d(v.dims[0], v.dims[1],
        const_cast<fftw_complex*>(reinterpret_cast<const fftw_complex*>(v.data.data())),
        r.data.data(), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    return r;
}
#endif

// Perform the convolution of two 1D arrays, assuming that they are based on the same 'x' coordinate
template<typename TypeX, typename TypeY1, typename TypeY2>
auto convolve(const vec<1,TypeX>& x, const vec<1,TypeY1>& y1, const vec<1,TypeY2>& y2) ->
    vec<1,decltype(x[0]*y1[0]*y2[0])> {

    phypp_check(x.dims == y1.dims, "incompatible dimensions for X and Y1 "
        "(", x.dims, " vs. ", y1.dims, ")");
    phypp_check(x.dims == y2.dims, "incompatible dimensions for X and Y2 "
        "(", x.dims, " vs. ", y2.dims, ")");
    phypp_check(x.size() > 3, "convolve needs arrays of at least 3 elements to work "
        "(got ", x.size(), ")");

    vec<1,decltype(x[0]*y1[0]*y2[0])> r(x.size());
    for (uint_t i = 0; i < x.size(); ++i) {
        auto dx = x.safe[i];
        if (i == 0) dx = x.safe[i+1] - dx;
        else dx -= x.safe[i-1];

        auto tmp = interpolate(y2, x + x.safe[i], x);
        r += y1.safe[i]*tmp*dx;
    }

    return r;
}

struct no_check {};

template <typename T>
struct convex_hull {
    // Default constructed hull is empty
    convex_hull() {}

    // Build a new convex hull and make sure it is valid
    convex_hull(vec<1,T> tx, vec<1,T> ty) : x(std::move(tx)), y(std::move(ty)) {
        validate();
    }

    // Build a hull without security checks
    // example: convex_hull(x, y, no_check{});
    convex_hull(vec<1,T> tx, vec<1,T> ty, no_check) :
        x(std::move(tx)), y(std::move(ty)), validated(true) {}

    // Number of vertices in the hull
    uint_t size() const {
        return x.size();
    }

    // Check the validity of this convex hull
    // Will only be done once unless 'force' is set to 'true'.
    void validate(bool force = false) const {
        if (!validated || force) {
            phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
                "(", x.dims, " vs. ", y.dims, ")");
            phypp_check(x.size() >= 3, "a hull must have at least 3 elements "
                "(got ", x.size(), ")");

            const auto eps = std::numeric_limits<T>::epsilon();
            const uint_t hend = size()-1;
            closed = (fabs(x.safe[0] - x.safe[hend]) < eps &&
                 fabs(y.safe[0] - y.safe[hend]) < eps);

            orient = ((x.safe[1] - x.safe[0])*(y.safe[2] - y.safe[1]) -
                     (y.safe[1] - y.safe[0])*(x.safe[2] - x.safe[1]) > 0) ? 1 : -1;

            validated = true;
        }
    }

    // Create the (ux, uy, nx, ny) vectors if they do not exist yet
    void build_vectors() const {
        if (!vectors_built) {
            ux.resize(size()-1); uy.resize(size()-1);
            nx.resize(size()-1); ny.resize(size()-1);
            length.resize(size()-1);

            for (uint_t i : range(size()-1)) {
                // Get unit vector
                ux.safe[i] = x.safe[i+1] - x.safe[i];
                uy.safe[i] = y.safe[i+1] - y.safe[i];
                // Get perpendicular vector, pointing inside the hull
                nx.safe[i] = orient == 1 ? y.safe[i+1] - y.safe[i] : y.safe[i] - y.safe[i+1];
                ny.safe[i] = orient == 1 ? x.safe[i] - x.safe[i+1] : x.safe[i+1] - x.safe[i];
                // Normalize
                auto l = sqrt(sqr(nx.safe[i]) + sqr(ny.safe[i]));
                length.safe[i] = l;
                ux.safe[i] /= l; uy.safe[i] /= l; nx.safe[i] /= l; ny.safe[i] /= l;
            }

            vectors_built = true;
        }
    }

    vec<1,T> x, y; // vertices

    // Cached computations
    mutable bool validated = false; // flag if the hull has been validated or not
    mutable bool closed = false;    // is first vertex = last vertex ?
    mutable int_t orient = 0;       // orientation of the hull, 1: ccw, -1: cw

    mutable bool vectors_built = false; // flag if segment data has been computed
    mutable vec<1,T> ux, uy;            // unit vector
    mutable vec<1,T> nx, ny;            // normal vector
    mutable vec<1,T> length;            // vector length
};

// Build the convex hull of a set of points, returning the indices of the points that form the hull
// in counter-clockwise order.
// Uses the monotone chain algorithm, taken from:
// http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C.2B.2B
template<typename T>
convex_hull<T> build_convex_hull(const vec<1,T>& x, const vec<1,T>& y) {
    phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
        "(", x.dims, " vs. ", y.dims, ")");

    uint_t npt = n_elements(x);
    vec1u res = uintarr(2*npt);
    vec1u ids = uindgen(npt);
    std::sort(ids.data.begin(), ids.data.end(), [&](uint_t i, uint_t j) {
        if (x.safe[i] < x.safe[j]) return true;
        if (x.safe[i] > x.safe[j]) return false;
        if (y.safe[i] < y.safe[j]) return true;
        return false;
    });

    auto cross = [&](uint_t i, uint_t j, uint_t k) {
        return (x.safe[j] - x.safe[i])*(y.safe[k] - y.safe[i])
             - (y.safe[j] - y.safe[i])*(x.safe[k] - x.safe[i]);
    };

    uint_t k = 0;
    for (uint_t i = 0; i < npt; ++i) {
        while (k >= 2 && cross(res.safe[k-2], res.safe[k-1], ids.safe[i]) <= 0) --k;
        res.safe[k] = ids.safe[i];
        ++k;
    }

    uint_t t = k+1;
    for (uint_t i = npt-2; i != npos; --i) {
        while (k >= t && cross(res.safe[k-2], res.safe[k-1], ids.safe[i]) <= 0) --k;
        res.safe[k] = ids.safe[i];
        ++k;
    }

    res.data.resize(k);
    res.dims[0] = k;

    convex_hull<T> hull(x.safe[res], y.safe[res], no_check{});
    hull.orient = 1; // counter-clockwise by construction
    hull.closed = true; // closed by construction

    return hull;
}

template<typename TX, typename TY, typename H>
bool in_convex_hull(const TX& x, const TY& y, const convex_hull<H>& hull) {
    hull.validate();
    phypp_check(hull.closed, "the provided hull must be closed");

    // Find out if the hull is built counter-clockwise or not
    for (uint_t i : range(hull.size()-1)) {
        auto cross = (hull.x.safe[i+1] - hull.x.safe[i])*(y - hull.y.safe[i]) -
                     (hull.y.safe[i+1] - hull.y.safe[i])*(x - hull.x.safe[i]);
        if (cross*hull.orient < 0) return false;
    }

    return true;
}

template<std::size_t Dim, typename TX, typename TY, typename H>
vec<Dim,bool> in_convex_hull(const vec<Dim,TX>& x, const vec<Dim,TY>& y,
    const convex_hull<H>& hull) {

    vec<Dim,bool> res = replicate(true, x.dims);

    hull.validate();
    phypp_check(hull.closed, "the provided hull must be closed");
    phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
        "(", x.dims, " vs. ", y.dims, ")");

    for (uint_t i : range(hull.size()-1)) {
        for (uint_t p : range(x)) {
            if (!res.safe[p]) continue;

            auto cross = (hull.x.safe[i+1] - hull.x.safe[i])*(y.safe[p] - hull.y.safe[i]) -
                         (hull.y.safe[i+1] - hull.y.safe[i])*(x.safe[p] - hull.x.safe[i]);

            if (cross*hull.orient < 0) res.safe[p] = false;
        }
    }

    return res;
}

// Compute the signed distance of a set of points with respect to the provided convex
// hull. Positive distances mean that the point lies inside the hull.
template<typename TX, typename TY, typename H, typename enable =
    typename std::enable_if<!is_vec<TX>::value && !is_vec<TY>::value>::type>
auto convex_hull_distance(const TX& x, const TY& y, const convex_hull<H>& hull)
    -> decltype(sqrt(x*y)) {

    hull.validate();
    phypp_check(hull.closed, "the provided hull must be closed");
    hull.build_vectors();

    decltype(sqrt(x*y)) res = finf;
    bool inhull = true;

    for (uint_t i : range(hull.size()-1)) {
        // Compute signed distance to current hull face line
        auto dx = x - hull.x.safe[i];
        auto dy = y - hull.y.safe[i];
        auto d  = dx*hull.nx.safe[i] + dy*hull.ny.safe[i];

        // Check if the point is inside the hull or not
        if (d > 0) {
            inhull = false;
        }

        d = fabs(d);
        if (res > d) {
            // Find if the projection of the point on the face lies on the segment
            auto proj = dx*hull.ux.safe[i] + dy*hull.uy.safe[i];
            if (proj < 0) {
                // Projection lies before starting point
                d = sqrt(sqr(dx) + sqr(dy));
            } else if (proj > hull.length.safe[i]) {
                // Projection lies after ending point
                d = sqrt(sqr(x - hull.x.safe[i+1]) + sqr(y - hull.y.safe[i+1]));
            }

            // Keep this distance if it is minimal
            if (res > d) res = d;
        }
    }

    // Give negative distance to points outside the hull
    if (!inhull) {
        res = -res;
    }

    return res;
}

// Compute the signed distance of a set of points with respect to the provided convex
// hull. Positive distances mean that the point lies inside the hull.
template<std::size_t Dim, typename TX, typename TY, typename H>
auto convex_hull_distance(const vec<Dim,TX>& x, const vec<Dim,TY>& y,
    const convex_hull<H>& hull) -> vec<Dim, decltype(sqrt(x[0]*y[0]))> {

    phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
        "(", x.dims, " vs. ", y.dims, ")");
    hull.validate();
    phypp_check(hull.closed, "the provided must be closed");
    hull.build_vectors();

    vec<Dim,decltype(sqrt(x[0]*y[0]))> res = replicate(finf, x.dims);
    vec<Dim,bool> inhull = replicate(true, x.dims);

    for (uint_t i : range(hull.size()-1)) {
        for (uint_t p : range(x)) {
            // Compute signed distance to current hull face line
            auto dx = x.safe[p] - hull.x.safe[i];
            auto dy = y.safe[p] - hull.y.safe[i];
            auto d  = dx*hull.nx.safe[i] + dy*hull.ny.safe[i];

            // Check if the point is inside the hull or not
            if (d > 0) {
                inhull.safe[p] = false;
            }

            d = fabs(d);
            if (res.safe[p] > d) {
                // Find if the projection of the point on the face lies on the segment
                auto proj = dx*hull.ux.safe[i] + dy*hull.uy.safe[i];
                if (proj < 0) {
                    // Projection lies before starting point
                    d = sqrt(sqr(dx) + sqr(dy));
                } else if (proj > hull.length.safe[i]) {
                    // Projection lies after ending point
                    d = sqrt(sqr(x.safe[p] - hull.x.safe[i+1]) + sqr(y.safe[p] - hull.y.safe[i+1]));
                }

                // Keep this distance if it is minimal
                if (res.safe[p] > d) res.safe[p] = d;
            }
        }
    }

    // Give negative distance to points outside the hull
    for (uint_t p : range(inhull)) {
        if (!inhull.safe[p]) {
            res.safe[p] *= -1.0;
        }
    }

    return res;
}

// Compute the angular distance between two RA/Dec positions [radian].
// Assumes that RA & Dec coordinates are in radian.
double angdistr(double ra1, double dec1, double ra2, double dec2) {
    double sra = sin(0.5*(ra2 - ra1));
    double sde = sin(0.5*(dec2 - dec1));
    return 2.0*asin(sqrt(sde*sde + sra*sra*cos(dec2)*cos(dec1)));
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
double angdist(double tra1, double tdec1, double tra2, double tdec2) {
    const double d2r = dpi/180.0;
    double ra1 = d2r*tra1, ra2 = d2r*tra2, dec1 = d2r*tdec1, dec2 = d2r*tdec2;
    double sra = sin(0.5*(ra2 - ra1));
    double sde = sin(0.5*(dec2 - dec1));
    return 3600.0*2.0*asin(sqrt(sde*sde + sra*sra*cos(dec2)*cos(dec1)))/d2r;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1, typename TR2, typename TD2>
vec<N,double> angdist(const vec<N,TR1>& tra1, const vec<N,TD1>& tdec1,
    const vec<N,TR2>& tra2, const vec<N,TD2>& tdec2) {
    phypp_check(tra1.dims == tdec1.dims, "first RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");
    phypp_check(tra2.dims == tdec2.dims, "second RA and Dec dimensions do not match (",
        tra2.dims, " vs ", tdec2.dims, ")");
    phypp_check(tra1.dims == tra2.dims, "position sets dimensions do not match (",
        tra1.dims, " vs ", tra2.dims, ")");

    vec<N,double> res(tra1.dims);
    for (uint_t i : range(tra1)) {
        res.safe[i] = angdist(tra1.safe[i], tdec1.safe[i], tra2.safe[i], tdec2.safe[i]);
    }

    return res;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1>
vec<N,double> angdist(const vec<N,TR1>& tra1, const vec<N,TD1>& tdec1,
    double tra2, double tdec2) {
    phypp_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");

    vec<N,double> res(tra1.dims);
    for (uint_t i : range(tra1)) {
        res.safe[i] = angdist(tra1.safe[i], tdec1.safe[i], tra2, tdec2);
    }

    return res;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1>
vec<N,bool> angdist_less(const vec<N,TR1>& tra1, const vec<N,TD1>& tdec1,
    double tra2, double tdec2, double radius) {
    phypp_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");

    const double d2r = dpi/180.0;
    const double rrad = d2r*radius/3600.0;
    const double crad = sqr(sin(rrad/2.0));

    vec<N,bool> res(tra1.dims);
    for (uint_t i : range(tra1)) {
        double ra1 = d2r*tra1.safe[i], ra2 = d2r*tra2;
        double dec1 = d2r*tdec1.safe[i], dec2 = d2r*tdec2;

        double sra = sin(0.5*(ra2 - ra1));
        double sde = sin(0.5*(dec2 - dec1));

        res.safe[i] = sqr(sde) + sqr(sra)*cos(dec2)*cos(dec1) < crad;
    }

    return res;
}

// Move a point in RA/Dec [degree] by an increment in [arcsec]
// This is an approximation for small movements of the order of a degree or less.
void move_ra_dec(double& ra, double& dec, double dra, double ddec) {
    ra += dra/cos(dec*dpi/180.0)/3600.0;
    dec += ddec/3600.0;

    // Note: another, in principle more accurate implementation is
    // ra *= cos(dec*dpi/180.0);
    // dec += ddec/3600.0;
    // ra += dra/3600.0;
    // ra /= cos(dec*dpi/180.0);
    // But this implementation suffers from numerical instability!
    // The culprit is the cos(dec_old)/cos(dec_new) ratio.
}

template<std::size_t D, typename T, typename U>
void move_ra_dec(vec<D,T>& ra, vec<D,U>& dec, double dra, double ddec) {
    phypp_check(ra.dims == dec.dims, "RA and Dec dimensions do not match (",
        ra.dims, " vs ", dec.dims, ")");

    dra /= 3600.0;
    ddec /= 3600.0;
    const double d2r = dpi/180.0;

    for (uint_t i : range(ra)) {
        ra.safe[i] += dra/cos(d2r*dec.safe[i]);
        dec.safe[i] += ddec;
    }
}

// Adjust RA and Dec coordinates to remain within the range [0,360] and [-90,90].
// The geometry is not modified, this is just a convention for the standard
// representation of the coordinates.
template<typename T, typename U, typename enable = typename std::enable_if<
    !is_vec<T>::value && !is_vec<U>::value>::type>
void normalize_coordinates(T& ra, U& dec) {
    while (dec > 90.0) {
        dec = 180.0 - dec;
        ra += 180.0;
    }
    while (dec < -90.0) {
        dec = -180.0 - dec;
        ra += 180.0;
    }

    while (ra > 360.0) {
        ra -= 360.0;
    }
    while (ra < 0.0) {
        ra += 360.0;
    }
}

template<std::size_t D, typename T, typename U>
void normalize_coordinates(vec<D,T>& ra, vec<D,U>& dec) {
    phypp_check(ra.dims == dec.dims, "RA and Dec dimensions do not match (",
        ra.dims, " vs ", dec.dims, ")");

    for (auto i : range(ra)) {
        normalize_coordinates(ra.safe[i], dec.safe[i]);
    }
}

// Return the rotation angle in degrees between the equator and the geodesic that passes
// through the celestial origin (i.e., the point of coordinates RA=0 and Dec=0) and the
// provided position.
template <typename T=void>
double vernal_angle(double ra, double dec) {
    ra *= dpi/180;
    dec *= dpi/180;

    double p = sqr(sin(dec/2)) + cos(dec)*sqr(sin(ra/2));
    double e = (180/dpi)*asin(sin(dec)/(2*(sqrt(p-p*p))));

    if (ra > 180.0) e = -e;

    return e;
}

// Rotate one set of RA/Dec positions around the vernal equinox by an angle 'e' given in
// degrees. NB: the vernal equinox is the line that goes from the center of Earth and
// passes through the celestial origin (i.e., the point of coordinates RA=0 and Dec=0).
// It is the 'x' axis in cartesian coordinates.
template <typename T1, typename T2, typename T3, typename T4>
void angrot_vernal(const T1& tr1, const T2& td1, double e, T3& r2, T4& d2) {
    phypp_check(tr1.dims == td1.dims, "incompatible dimensions between RA and Dec "
        "(", tr1.dims, " vs. ", td1.dims, ")");

    e *= dpi/180.0;
    double ce = cos(e);
    double se = sin(e);

    auto r1 = tr1*dpi/180.0;
    auto d1 = td1*dpi/180.0;

    auto cr1 = cos(r1), sr1 = sin(r1);
    auto cd1 = cos(d1), sd1 = sin(d1);

    auto x1 = cd1*cr1, y1 = cd1*sr1, z1 = sd1;

    auto y2 = ce*y1 - se*z1;
    auto z2 = ce*z1 + se*y1;

    d2 = asin(z2);
    r2 = atan2(y2, x1);

    // Other solutions
    // d2 = asin(sin(d1)*ce + cos(d1)*sin(r1)*se);
    // r2 = acos(cos(r1)*cos(d1)/cos(d2));
    // r2 = acos(x1/sqrt(1-sqr(z2)));

    r2 *= 180/dpi;
    d2 *= 180/dpi;

    // Come back to [0, 360]
    for (auto& r : r2) {
        if (r < 0) r += 360.0;
    }
}

#endif
