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
vec_t<1,T> rgen(T i, U j) {
    if (i < T(j)) {
        uint_t n = j-i+1;
        vec_t<1,T> v(n);
        for (uint_t k : range(n)) {
            v.safe[k] = i+k;
        }
        return v;
    } else {
        uint_t n = i-j+1;
        vec_t<1,T> v(n);
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
vec_t<2,T> make_bins(T mi, T ma) {
    return {{mi}, {ma}};
}

template<typename T>
vec_t<2,T> make_bins(T mi, T ma, uint_t n) {
    vec_t<2,T> b(2,n);
    T d = (ma - mi)/n;
    for (uint_t i : range(n)) {
        b.safe(0,i) = mi + i*d;
        b.safe(1,i) = mi + (i+1)*d;
    }

    return b;
}

template<typename T = double>
vec_t<2,T> make_bins(const vec_t<1,T>& v) {
    vec_t<2,T> b(2, v.size()-1);
    for (uint_t i : range(b.dims[1])) {
        b.safe(0,i) = v.safe[i];
        b.safe(1,i) = v.safe[i+1];
    }

    return b;
}

template<std::size_t Dim, typename Type, typename B = double>
vec_t<Dim,bool> in_bin(const vec_t<Dim,Type>& v, const vec_t<1,B>& b) {
    phypp_check(b.dims[0] == 2, "can only be called on a single bin "
        "vector (expected dims=[2], got dims=[", b.dims, "])");

    auto low = b.safe[0];
    auto up  = b.safe[1];
    vec_t<Dim,bool> res(v.dims);
    for (uint_t i : range(v)) {
        res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
    }

    return res;
}

template<std::size_t Dim, typename Type, typename B>
vec_t<Dim,bool> in_bin(const vec_t<Dim,Type>& v, const vec_t<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    auto low = b.safe(0,ib);
    auto up  = b.safe(1,ib);
    vec_t<Dim,bool> res(v.dims);
    for (uint_t i : range(v)) {
        res.safe[i] = v.safe[i] >= low && v.safe[i] < up;
    }

    return res;
}

template<std::size_t Dim, typename Type, typename B>
vec_t<Dim,bool> in_bin_open(const vec_t<Dim,Type>& v, const vec_t<2,B>& b, uint_t ib) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");
    phypp_check(ib < b.dims[1], "bin index is out of bounds ",
        "(", ib, " vs. ", b.dims[1], ")");

    auto low = b.safe(0,ib);
    auto up  = b.safe(1,ib);
    vec_t<Dim,bool> res(v.dims);
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
vec_t<1,rtype_t<Type>> bin_center(const vec_t<2,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");

    return 0.5*(b.safe(1,_) + b.safe(0,_));
}

template<typename Type>
vec_t<1,rtype_t<Type>> bin_width(const vec_t<2,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=[2, ...], got dims=[", b.dims, "])");

    return b.safe(1,_) - b.safe(0,_);
}

template<typename Type>
rtype_t<Type> bin_center(const vec_t<1,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=2, got dims=", b.dims, ")");

    return 0.5*(b.safe[1] + b.safe[0]);
}

template<typename Type>
rtype_t<Type> bin_width(const vec_t<1,Type>& b) {
    phypp_check(b.dims[0] == 2, "B is not a bin vector "
        "(expected dims=2, got dims=", b.dims, ")");

    return b.safe[1] - b.safe[0];
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool finite(const T& t) {
    return std::isfinite(t);
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> finite(const vec_t<Dim,Type>& v) {
    vec_t<Dim,bool> r(v.dims);
    for (uint_t i : range(v)) {
        r.safe[i] = std::isfinite(v.safe[i]);
    }

    return r;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
bool nan(const T& t) {
    return std::isnan(t);
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> nan(const vec_t<Dim,Type>& v) {
    vec_t<Dim,bool> r(v.dims);
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
vec_t<dim_total<Args...>::value,double> randomn(T& seed, Args&& ... args) {
    vec_t<dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
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
vec_t<dim_total<Args...>::value,double> randomu(T& seed, Args&& ... args) {
    vec_t<dim_total<Args...>::value,double> v(std::forward<Args>(args)...);
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
    vec_t<dim_total<Args...>::value,decltype(mi+ma)> {
    auto v = randomu(seed, std::forward<Args>(args)...);
    using rtype = decltype(mi + ma);
    return vec_t<vec_dim<decltype(v)>::value,rtype>(v*(ma + 1 - mi) + mi);
}

template<typename T, typename TypeX, typename TypeY, typename ... Args>
vec_t<dim_total<Args...>::value,rtype_t<TypeX>> random_pdf(T& seed, const vec_t<1,TypeX>& px,
    const vec_t<1,TypeY>& py, Args&& ... args) {

    // TODO: make an alternative version for integers using std::discrete_distribution.

    using rtype = rtype_t<TypeX>;
    vec_t<dim_total<Args...>::value,rtype> v(std::forward<Args>(args)...);
    std::piecewise_linear_distribution<rtype> distribution(px.begin(), px.end(), py.begin());
    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed);
    }

    return v;
}

template<typename T>
bool random_coin(T& seed, double prob) {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(seed) <= prob;
}

template<typename T, typename ... Args>
vec_t<dim_total<Args...>::value,bool> random_coin(T& seed, double prob, Args&& ... args) {
    vec_t<dim_total<Args...>::value,bool> v(std::forward<Args>(args)...);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (uint_t i : range(v)) {
        v.safe[i] = distribution(seed) <= prob;
    }

    return v;
}

template<std::size_t Dim, typename Type, typename T>
vec_t<Dim,Type> shuffle(vec_t<Dim,Type> v, T& seed) {
    std::shuffle(v.begin(), v.end(), seed);
    return v;
}

template<typename F, F f, std::size_t Dim, typename Type, typename ... Args>
auto run_index_(const vec_t<Dim,Type>& v, const Args& ... args, uint_t dim) ->
vec_t<Dim-1,typename return_type<F>::type> {

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
            tmp[j] = v.safe[base + j*mpitch];
        }

        r[i] = (*f)(tmp, args...);
    }

    return r;
}

// Iterate over one dimension of the provided vector and call a function for each slice.
// Note: the slice provided to the function is mutable. It is allowed to modify its elements or
// swap them, but the slice should never be resized. No check will be performed.
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
            tmp[j] = v.safe[base + j*mpitch];
        }

        func(i, tmp);
    }
}

template<std::size_t Dim, typename Type>
vec_t<1,rtype_t<Type>> run_dim_idx_apply_ids_(const vec1u& ids, const vec_t<Dim,Type>& v) {
    return v.safe[ids].concretise();
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
auto run_dim(const vec_t<Dim,Type>& v, uint_t dim, F&& func) ->
    vec_t<Dim-1,typename return_type<F>::type> {

    vec_t<Dim-1,typename return_type<F>::type> r;
    for (uint_t i = 0; i < dim; ++i) {
        r.dims[i] = v.dims[i];
    }
    for (uint_t i = dim+1; i < Dim; ++i) {
        r.dims[i-1] = v.dims[i];
    }

    r.resize();

    run_dim_idx(v, dim, [&](uint_t i, vec_t<1, rtype_t<Type>>& tv) {
        r[i] = func(tv);
    });

    return r;
}

template<typename T>
using total_return_type = typename std::conditional<std::is_integral<T>::value,
    typename std::conditional<std::is_unsigned<T>::value, uint_t, int_t>::type,
    double>::type;

template<std::size_t Dim, typename Type>
total_return_type<rtype_t<Type>> total(const vec_t<Dim,Type>& v) {
    total_return_type<rtype_t<Type>> total = 0;
    for (auto& t : v) {
        total += t;
    }

    return total;
}


template<std::size_t Dim, typename Type, typename enable =
    typename std::enable_if<std::is_same<rtype_t<Type>, bool>::value>::type>
uint_t count(const vec_t<Dim,Type>& v) {
    uint_t n = 0u;
    for (bool b : v) {
        if (b) ++n;
    }

    return n;
}


template<std::size_t Dim, typename T>
double fraction_of(const vec_t<Dim,T>& b) {
    return mean(b);
}

template<std::size_t Dim, typename Type>
double mean(const vec_t<Dim,Type>& v) {
    double total = 0.0;
    for (auto& t : v) {
        total += t;
    }

    return total/v.size();
}

template<std::size_t Dim, typename Type>
rtype_t<Type> inplace_median(vec_t<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the median of an empty vector");

    using vtype = typename vec_t<Dim,Type>::vtype;
    using dtype = typename vtype::value_type;

    uint_t nwrong = 0;
    for (auto& t : v) {
        nwrong += nan(t);
    }

    if (nwrong == v.size()) return dnan;

    std::ptrdiff_t offset = (v.size()-nwrong)/2;
    std::nth_element(v.data.begin(), v.data.begin() + offset, v.data.end(),
        [](dtype i, dtype j) {
            if (nan(dref<Type>(i))) return false;
            if (nan(dref<Type>(j))) return true;
            return dref<Type>(i) < dref<Type>(j);
        }
    );

    return *(v.begin() + offset);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> median(vec_t<Dim,Type> v) {
    return inplace_median(v);
}

template<std::size_t Dim, typename Type, typename U>
rtype_t<Type> percentile(const vec_t<Dim,Type>& v, const U& u) {
    phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

    vec1u ok = where(finite(v));
    if (ok.empty()) return 0;

    // TODO: use same algorithm than median
    typename vec_t<1,Type>::effective_type t = v.safe[ok];
    std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    return *(t.begin() + offset);
}

template<std::size_t Dim, typename Type>
void percentiles_(vec_t<1,Type>& r, uint_t i, vec_t<Dim,Type>& t) {}

template<std::size_t Dim, typename Type, typename U, typename ... Args>
void percentiles_(vec_t<1,Type>& r, uint_t i, vec_t<Dim,Type>& t, const U& u, const Args& ... args) {
    std::ptrdiff_t offset = clamp(t.size()*u, 0u, t.size()-1);
    std::nth_element(t.begin(), t.begin() + offset, t.end());
    r.safe[i] = *(t.begin() + offset);
    ++i;

    percentiles_(r, i, t, args...);
}

template<std::size_t Dim, typename Type, typename ... Args>
typename vec_t<1,Type>::effective_type percentiles(const vec_t<Dim,Type>& v, const Args& ... args) {
    phypp_check(!v.empty(), "cannot find the percentiles of an empty vector");

    vec1u ok = where(finite(v));
    typename vec_t<1,Type>::effective_type t;
    if (ok.empty()) return t;
    t = v.safe[ok];

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
    return v > p.safe[0] && v < p.safe[1];
}

template<std::size_t Dim, typename Type>
vec_t<Dim,bool> mad_clip(const vec_t<Dim,Type>& tv, double sigma) {
    auto v = tv.concretise();
    auto med = inplace_median(v);
    auto mad = median(fabs(v - med));
    return fabs(tv - med) < sigma*mad;
}

template<std::size_t Dim, typename Type>
typename vec_t<Dim,Type>::const_iterator min_(const vec_t<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the minimum of an empty vector");

    auto iter = std::min_element(v.begin(), v.end(), [](rtype_t<Type> t1, rtype_t<Type> t2){
        if (nan(t1)) return false;
        if (nan(t2)) return true;
        return t1 < t2;
    });

    if (iter == v.end()) iter = v.begin();
    return iter;
}

template<std::size_t Dim, typename Type>
typename vec_t<Dim,Type>::const_iterator max_(const vec_t<Dim,Type>& v) {
    phypp_check(!v.empty(), "cannot find the maximum of an empty vector");

    auto iter = std::max_element(v.begin(), v.end(), [](rtype_t<Type> t1, rtype_t<Type> t2){
        if (nan(t1)) return true;
        if (nan(t2)) return false;
        return t1 < t2;
    });

    if (iter == v.end()) iter = v.begin();
    return iter;
}

template<std::size_t Dim, typename Type>
rtype_t<Type> min(const vec_t<Dim,Type>& v) {
    return *min_(v);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> max(const vec_t<Dim,Type>& v) {
    return *max_(v);
}

template<std::size_t Dim, typename Type>
rtype_t<Type> min(const vec_t<Dim,Type>& v, uint_t& id) {
    auto iter = min_(v);
    id = iter - v.begin();
    return *iter;
}

template<std::size_t Dim, typename Type>
rtype_t<Type> max(const vec_t<Dim,Type>& v, uint_t& id) {
    auto iter = max_(v);
    id = iter - v.begin();
    return *iter;
}

template<std::size_t Dim, typename Type>
uint_t min_id(const vec_t<Dim,Type>& v) {
    return min_(v) - v.begin();
}

template<std::size_t Dim, typename Type>
uint_t max_id(const vec_t<Dim,Type>& v) {
    return max_(v) - v.begin();
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const vec_t<Dim,Type1>& v1, const vec_t<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "min: incompatible vector dimensions "
        "(", v1.dims, " vs. ", v2.dims, ")");

    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::min<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const vec_t<Dim,Type1>& v1, const vec_t<Dim,Type2>& v2) {
    phypp_check(v1.dims == v2.dims, "max: incompatible vector dimensions "
        "(", v1.dims, " vs. ", v2.dims, ")");

    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::max<decltype(v1[0]*v2[0])>(v1.safe[i], v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const vec_t<Dim,Type1>& v1, const Type2& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::min<decltype(v1[0]*v2)>(v1.safe[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const vec_t<Dim,Type1>& v1, const Type2& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v1.dims);
    for (uint_t i = 0; i < v1.size(); ++i) {
        r.safe[i] = std::max<decltype(v1[0]*v2)>(v1.safe[i], v2);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> min(const Type1& v1, const vec_t<Dim,Type2>& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r.safe[i] = std::min<decltype(v1*v2[0])>(v1, v2.safe[i]);
    }
    return r;
}

template<std::size_t Dim, typename Type1, typename Type2>
vec_t<Dim,rtype_t<Type1>> max(const Type1& v1, const vec_t<Dim,Type2>& v2) {
    vec_t<Dim,rtype_t<Type1>> r = arr<rtype_t<Type1>>(v2.dims);
    for (uint_t i = 0; i < v2.size(); ++i) {
        r.safe[i] = std::max<decltype(v1*v2[0])>(v1, v2.safe[i]);
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
double stddev(const vec_t<Dim,Type>& v) {
    return rms(v - mean(v));
}

template<std::size_t Dim, typename Type>
double mad(const vec_t<Dim,Type>& v) {
    return median(fabs(v - median(v)));
}

#define RUN_INDEX(func) \
    struct func ## _run_index_wrapper_ { \
        template<typename T, typename ... Args> \
        static auto run(const vec_t<1,T>& v, Args&& ... args) -> \
        decltype(func(v, std::forward<Args>(args)...)) { \
            return func(v, std::forward<Args>(args)...); \
        } \
    }; \
    \
    template<std::size_t Dim, typename Type, typename ... Args> \
    auto func(const vec_t<Dim,Type>& v, Args&& ... args, uint_t dim) -> \
    vec_t<Dim-1, decltype(func(std::declval<vec_t<1,rtype_t<Type>>>(), std::forward<Args>(args)...))> { \
        using fptr = decltype(func(std::declval<vec_t<1,rtype_t<Type>>>(), std::forward<Args>(args)...)) \
            (*)(const vec_t<1,rtype_t<Type>>&, Args&& ...); \
        using wrapper = func ## _run_index_wrapper_; \
        return run_index_<fptr, &wrapper::run<rtype_t<Type>, Args...>>(v, dim); \
    }

RUN_INDEX(total);
RUN_INDEX(count);
RUN_INDEX(fraction_of);
RUN_INDEX(mean);
RUN_INDEX(median);
RUN_INDEX(percentile);
RUN_INDEX(rms);
RUN_INDEX(stddev);
RUN_INDEX(mad);

#undef RUN_INDEX

template<std::size_t Dim, typename Type, typename TypeB>
vec1u histogram(const vec_t<Dim,Type>& data, const vec_t<2,TypeB>& bins) {
    phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2,...], got dims=[", bins.dims, "])");

    using rtype = rtype_t<Type>;
    vec_t<Dim,rtype> tmp = data;

    uint_t nbin = bins.dims[1];
    vec1u counts(nbin);

    auto first = tmp.data.begin();
    for (uint_t i : range(nbin)) {
        auto last = std::partition(first, tmp.data.end(), [&bins,i](rtype t) {
            return t >= bins.safe(0,i) && t < bins.safe(1,i);
        });

        if (last == tmp.data.end()) break;

        counts.safe[i] = last - first;
        first = last;
    }

    return counts;
}

template<std::size_t Dim, typename Type, typename TypeB, typename TypeW>
vec_t<1,rtype_t<TypeW>> histogram(const vec_t<Dim,Type>& data, const vec_t<Dim,TypeW>& weight,
    const vec_t<2,TypeB>& bins) {
    phypp_check(bins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", bins.dims, "])");
    phypp_check(data.dims == weight.dims, "incompatible dimensions for data and weight "
        "(", data.dims, " vs. ", weight.dims, ")");

    vec1u tmp = uindgen(data.size());

    uint_t nbin = bins.dims[1];
    vec_t<1,rtype_t<TypeW>> counts(nbin);

    auto first = tmp.data.begin();
    for (uint_t i : range(nbin)) {
        auto last = std::partition(first, tmp.data.end(), [&bins,&data,i](uint_t id) {
            return data.safe[id] >= bins.safe(0,i) && data.safe[id] < bins.safe(1,i);
        });

        if (last == tmp.data.end()) break;

        for (; first != last; ++first) {
            counts.safe[i] += weight.safe[*first];
        }
    }

    return counts;
}

template<std::size_t Dim, typename TypeX, typename TypeY, typename TypeBX, typename TypeBY>
vec2u histogram2d(const vec_t<Dim,TypeX>& x, const vec_t<Dim,TypeY>& y,
    const vec_t<2,TypeBX>& xbins, const vec_t<2,TypeBY>& ybins) {

    phypp_check(xbins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", xbins.dims, "])");
    phypp_check(ybins.dims[0] == 2, "can only be called with a bin vector (expected "
        "dims=[2, ...], got dims=[", ybins.dims, "])");
    phypp_check(x.dims == y.dims, "incompatible dimensions for x and y (", x.dims, " vs. ",
        y.dims, ")");

    vec1u ids = uindgen(x.size());

    uint_t nxbin = xbins.dims[1];
    uint_t nybin = ybins.dims[1];
    vec2u counts(nxbin, nybin);

    auto firstx = ids.data.begin();
    for (uint_t i : range(nxbin)) {
        auto lastx = std::partition(firstx, ids.data.end(), [&xbins,&x,i](uint_t id) {
            return x.safe[id] >= xbins.safe(0,i) && x.safe[id] < xbins.safe(1,i);
        });

        if (lastx == ids.data.end()) break;

        auto firsty = firstx;
        for (uint_t j : range(nybin)) {
            auto lasty = std::partition(firsty, lastx, [&ybins,&y,j](uint_t id) {
                return y.safe[id] >= ybins.safe(0,j) && y.safe[id] < ybins.safe(1,j);
            });

            if (lasty == lastx) break;

            counts.safe(i,j) = lasty - firsty;

            firsty = lasty;
        }

        firstx = lastx;
    }

    return counts;
}

template<std::size_t Dim, typename Type>
void data_info_(const vec_t<Dim,Type>& v) {
    vec1u idok = where(finite(v));
    print(idok.size(), "/", v.size(), " valid values (dims: ", v.dims, ")");
    if (idok.size() == 0) return;

    vec_t<1,rtype_t<Type>> tv = v.safe[idok];

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
auto sign(const T& t) -> decltype(2*(t > 0) - 1) {
    return 2*(t > 0) - 1;
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
    auto name(const vec_t<Dim,Type>& v, const Args& ... args) -> \
        vec_t<Dim,decltype(name(v[0], args...))> { \
        using ntype = decltype(name(v[0], args...)); \
        vec_t<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(name(dref<Type>(t), args...)); \
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
    auto name(const vec_t<Dim,Type>& v, const Args& ... args) -> \
        vec_t<Dim,decltype(orig(v[0], args...))> { \
        using ntype = decltype(orig(v[0], args...)); \
        vec_t<Dim,ntype> r; r.dims = v.dims; r.data.reserve(v.size()); \
        for (auto& t : v.data) { \
            r.data.push_back(orig(dref<Type>(t), args...)); \
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
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];
    const uint_t m = b.dims[1];

    vec_t<2,ntype_t> r(n,m);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t j = 0; j < m; ++j)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i,j) += a.safe(i,k)*b.safe(k,j);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec_t<2,TypeA>& a, const vec_t<1,TypeB>& b) -> vec_t<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying matrix by vector "
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[1];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[0];

    vec_t<1,ntype_t> r(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i) += a.safe(i,k)*b.safe(k);
    }

    return r;
}

template<typename TypeA, typename TypeB>
auto mmul(const vec_t<1,TypeB>& b, const vec_t<2,TypeA>& a) -> vec_t<1,decltype(a(0,0)*b(0,0))> {
    phypp_check(a.dims[1] == b.dims[0], "wrong dimensions multiplying vector by matrix "
        "(", a.dims, " x ", b.dims, ")");

    const uint_t o = a.dims[0];

    using ntype_t = decltype(a(0,0)*b(0,0));
    const uint_t n = a.dims[1];

    vec_t<1,ntype_t> r = arr<ntype_t>(n);
    for (uint_t i = 0; i < n; ++i)
    for (uint_t k = 0; k < o; ++k) {
        r.safe(i) += b.safe(k)*a.safe(k,i);
    }

    return r;
}

template<typename Type>
void mprint(const vec_t<2,Type>& m) {
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
vec_t<2,rtype_t<Type>> transpose(const vec_t<2,Type>& v) {
    vec_t<2,rtype_t<Type>> r(v.dims);
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
auto diagonal(const vec_t<2,Type>& v) -> decltype(v(_,0)) {
    phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
        v.dims, ")");

    decltype(v(_,0)) d(vec_ref_tag, vec_access::get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i = 0; i < v.dims[0]; ++i) {
        d.safe[i] = ptr<Type>(v.safe(i,i));
    }

    return d;
}

template<typename Type>
auto diagonal(vec_t<2,Type>& v) -> decltype(v(_,0)) {
    phypp_check(v.dims[0] == v.dims[1], "can only be called on square matrix (got ",
        v.dims, ")");

    decltype(v(_,0)) d(vec_ref_tag, vec_access::get_parent(v));
    d.dims[0] = v.dims[0];
    d.resize();
    for (uint_t i = 0; i < v.dims[0]; ++i) {
        d.data[i] = ptr<Type>(v.safe(i,i));
    }

    return d;
}

template<typename Type = double>
vec_t<2,Type> identity_matrix(uint_t dim) {
    vec_t<2,Type> m(dim, dim);
    diagonal(m) = 1;
    return m;
}

template<typename TX, typename TY>
auto scale_matrix(const TX& sx, const TY& sy) -> vec_t<2,decltype(sx*sy)> {
    vec_t<2,decltype(sx*sy)> m(3, 3);
    m.safe(0,0) = sx;
    m.safe(1,1) = sy;
    m.safe(2,2) = 1;
    return m;
}

template<typename T>
vec_t<2,T> scale_matrix(const T& s) {
    vec_t<2,T> m(3, 3);
    m.safe(0,0) = s;
    m.safe(1,1) = s;
    m.safe(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto translation_matrix(const TX& tx, const TY& ty) -> vec_t<2,decltype(tx*ty)> {
    vec_t<2,decltype(tx*ty)> m(3, 3);
    diagonal(m) = 1;
    m.safe(0,2) = tx;
    m.safe(1,2) = ty;
    return m;
}

template<typename A>
auto rotation_matrix(const A& a) -> vec_t<2,decltype(cos(a))> {
    vec_t<2,decltype(cos(a))> m(3, 3);
    auto ca = cos(a), sa = sin(a);
    m.safe(0,0) = m.safe(1,1) = ca;
    m.safe(0,1) = -sa;
    m.safe(1,0) = sa;
    m.safe(2,2) = 1;
    return m;
}

template<typename TX, typename TY>
auto point2d(const TX& x, const TY& y) -> vec_t<1,decltype(x*y)> {
    return vec_t<1,decltype(x*y)>{x, y, 1};
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
bool invert(vec2d& i) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
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
#endif
}

template<typename TypeA>
bool invert(const vec_t<2,TypeA>& a, vec2d& i) {
    i = a;
    return invert<TypeA>(i);
}

template<typename Dummy = void>
bool invert_symmetric(vec2d& i) {
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
#endif
}

template<typename TypeA>
bool invert_symmetric(const vec_t<2,TypeA>& a, vec2d& i) {
    i = a;
    return invert_symmetric<TypeA>(i);
}

template<typename Dummy = void>
bool solve_symmetric(vec2d& alpha, vec1d& beta) {
#ifdef NO_LAPACK
    static_assert(!std::is_same<Dummy,Dummy>::value, "LAPACK support has been disabled, "
        "please enable LAPACK to use this function");
#else
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
#endif
}

template<typename Dummy = void>
bool eigen_symmetric(vec2d& a, vec1d& vals) {
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
    return eigen_symmetric<Dummy>(vecs, vals);
}


// -----------------------
// end of LAPACK functions


template<typename Type>
void symmetrize(vec_t<2,Type>& alpha) {
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
    fr.errors = sqrt(diagonal(alpha));

    vec1d model(nm);
    for (uint_t m = 0; m < nm; ++m) {
        model.safe[m] = 0.0;
        for (uint_t i = 0; i < np; ++i) {
            model.safe[m] += fr.params.safe[i]*cache.safe(i,m);
        }
    }

    fr.chi2 = total(sqr(model*flatten(1.0/ye) - tmp));

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
linfit_result linfit_pack(const vec_t<Dim,TypeY>& y, const vec_t<Dim,TypeE>& ye,
    const vec_t<Dim+1,TypeX>& x) {
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

// Returns the position of the last value in the array that is less than or equal to 'x'.
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

// Return the position of the last value in 'v' that is less than or equal to 'x1' and
// the position of the first value in 'v' that is greater than 'x2'.
// Note: assumes that:
//  1) 'v' is sorted and does not contain NaN values,
//  2) 'x2' is greater than or equal to 'x1'.
template<typename T, typename U, typename Type>
std::array<uint_t,2> bounds(T x1, U x2, const vec_t<1,Type>& v) {
    auto iter = std::upper_bound(v.data.begin(), v.data.end(), x1,
        typename vec_t<1,Type>::comparator());

    std::array<uint_t,2> res;
    if (iter == v.data.begin()) {
        res[0] = npos;
    } else {
        iter--;
        res[0] = iter - v.data.begin();
    }

    iter = std::upper_bound(iter, v.data.end(), x2,
        typename vec_t<1,Type>::comparator());

    if (iter == v.data.end()) {
        res[1] = npos;
    } else {
        res[1] = iter - v.data.begin();
    }

    return res;
}

// Return the indices of all the values in the array that are equal to 'x'.
// Note: assumes that 'v' is sorted and does not contain NaN values.
template<typename T, std::size_t Dim, typename Type>
vec1u equal_range(T x, const vec_t<Dim,Type>& v) {
    auto res = std::equal_range(v.data.begin(), v.data.end(), x,
        typename vec_t<Dim,Type>::comparator());

    return uindgen(1 + (res.second - res.first)) + (res.first - v.data.begin());
}

// Check if a given array is sorted or not
template<typename Type>
bool is_sorted(const vec_t<1,Type>& v) {
    for (uint_t i = 0; i < v.size()-1; ++i) {
        if (v.safe[i] >= v.safe[i+1]) return false;
    }

    return true;
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
auto interpolate(const vec_t<DI,TypeY>& y, const vec_t<DI,TypeX1>& x, const vec_t<DX,TypeX2>& nx) ->
    vec_t<DX,decltype(y[0]*x[0])> {

    using rtypey = rtype_t<TypeY>;
    using rtypex = rtype_t<TypeX1>;

    phypp_check(y.size() == x.size(),
        "interpolate: 'x' and 'y' arrays must contain the same number of elements");
    phypp_check(y.size() >= 2,
        "interpolate: 'x' and 'y' arrays must contain at least 2 elements");

    uint_t nmax = x.size();
    vec_t<DX,decltype(y[0]*x[0])> r; r.reserve(nx.size());
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
auto interpolate(const vec_t<DI,TypeY>& y, const vec_t<DI,TypeX>& x, const T& nx) ->
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
rtype_t<Type> bilinear(const vec_t<2,Type>& map, double x, double y) {
    int_t tix = floor(x);
    int_t tiy = floor(y);
    double dx = x - tix;
    double dy = y - tiy;

    if (tix >= map.dims[0]-1) {
        tix = map.dims[0]-2;
        dx = x - tix;
    } else if (tix < 0) {
        tix = 0;
        dx = x - tix;
    }

    if (tiy >= map.dims[1]-1) {
        tiy = map.dims[1]-2;
        dy = y - tiy;
    } else if (tiy < 0) {
        tiy = 0;
        dy = y - tiy;
    }

    uint_t ix = tix;
    uint_t iy = tiy;

    return map.safe(ix,iy)*(1.0 - dx)*(1.0 - dy) + map.safe(ix,iy+1)*(1.0 - dx)*dy
        + map.safe(ix+1,iy)*dx*(1.0 - dy) + map.safe(ix+1,iy+1)*dx*dy;
}

template<typename Type>
vec_t<2,rtype_t<Type>> rebin(const vec_t<2,Type>& map, const vec1d& mx,
    const vec1d& my, const vec1d& x, const vec1d& y) {

    phypp_check(map.dims[0] == mx.size(), "incompatible size of MAP and MX (", map.dims,
        " vs. ", mx.size(), ")");
    phypp_check(map.dims[1] == my.size(), "incompatible size of MAP and MY (", map.dims,
        " vs. ", my.size(), ")");

    vec_t<2,rtype_t<Type>> v(x.size(), y.size());

    vec1d ux = interpolate(dindgen(mx.size()), mx, x);
    vec1d uy = interpolate(dindgen(my.size()), my, y);

    for (uint_t ix : range(ux))
    for (uint_t iy : range(uy)) {
        v.safe(ix,iy) = bilinear(map, ux.safe[ix], uy.safe[iy]);
    }

    return v;
}

template<typename TypeX, typename TypeY>
auto integrate(const vec_t<1,TypeX>& x, const vec_t<1,TypeY>& y)
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
auto integrate(const vec_t<2,TypeX>& x, const vec_t<1,TypeY>& y)
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
auto integrate(const vec_t<1,TypeX>& x, const vec_t<1,TypeY>& y, double x0, double x1)
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
auto integrate(const vec_t<2,TypeX>& x, const vec_t<1,TypeY>& y, double x0, double x1)
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

template<typename F, typename T, typename U,
    typename enable = typename std::enable_if<!is_vec<F>::value>::type>
auto integrate(F f, T x0, U x1,
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
template<typename T>
using complex_type = phypp::complex<T>;

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
auto convolve(const vec_t<1,TypeX>& x, const vec_t<1,TypeY1>& y1, const vec_t<1,TypeY2>& y2) ->
    vec_t<1,decltype(x[0]*y1[0]*y2[0])> {

    phypp_check(x.dims == y1.dims, "incompatible dimensions for X and Y1 "
        "(", x.dims, " vs. ", y1.dims, ")");
    phypp_check(x.dims == y2.dims, "incompatible dimensions for X and Y2 "
        "(", x.dims, " vs. ", y2.dims, ")");
    phypp_check(x.size() > 3, "convolve needs arrays of at least 3 elements to work "
        "(got ", x.size(), ")");

    vec_t<1,decltype(x[0]*y1[0]*y2[0])> r(x.size());
    for (uint_t i = 0; i < x.size(); ++i) {
        auto dx = x.safe[i];
        if (i == 0) dx = x.safe[i+1] - dx;
        else dx -= x.safe[i-1];

        auto tmp = interpolate(y2, x + x.safe[i], x);
        r += y1.safe[i]*tmp*dx;
    }

    return r;
}

// Build the convex hull of a set of points, returning the indices of the points that form the hull
// in counter-clockwise order.
// Uses the monotone chain algorithm, taken from:
// http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C.2B.2B
template<typename TX, typename TY>
vec1u convex_hull(const TX& x, const TY& y) {
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

    return res;
}

template<typename THX, typename THY>
bool is_hull_closed(const vec1u& hull, const THX& hx, const THY& hy) {
    phypp_check(hull.size() >= 3, "a hull must have at least 3 elements "
        "(got ", hull.size(), ")");
    phypp_check(hx.dims == hy.dims, "incompatible dimensions between HX and HY "
        "(", hx.dims, " vs. ", hy.dims, ")");

    using rtype = typename THX::rtype;
    const auto eps = std::numeric_limits<rtype>::epsilon();
    const uint_t hend = hull.size()-1;
    return hull.safe[0] == hull.safe[hend] ||
        (fabs(hx[hull.safe[0]] - hx[hull.safe[hend]]) < eps &&
         fabs(hy[hull.safe[0]] - hy[hull.safe[hend]]) < eps);
}

template<typename TX, typename TY, typename THX, typename THY>
bool in_convex_hull(const TX& x, const TY& y, const vec1u& hull, const THX& hx, const THY& hy) {
    phypp_check(is_hull_closed(hull, hx, hy), "the provided hull is not closed");

    // Find out if the hull is built counter-clockwise or not
    uint_t i0 = hull.safe[0], i1 = hull.safe[1], i2 = hull.safe[2];
    bool sign = (hx[i1] - hx[i0])*(hy[i2] - hy[i1]) - (hy[i1] - hy[i0])*(hx[i2] - hx[i1]) > 0;

    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull.safe[i], p2 = hull.safe[i+1];
        auto cross = (hx[p2] - hx[p1])*(y - hy[p1]) - (hy[p2] - hy[p1])*(x - hx[p1]);
        if ((cross < 0) == sign) return false;
    }

    return true;
}

template<std::size_t Dim, typename TX, typename TY, typename THX, typename THY>
vec_t<Dim,bool> in_convex_hull(const vec_t<Dim,TX>& x, const vec_t<Dim,TY>& y, const vec1u& hull,
    const THX& hx, const THY& hy) {

    phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
        "(", x.dims, " vs. ", y.dims, ")");
    phypp_check(is_hull_closed(hull, hx, hy), "the provided hull is not closed");

    // Find out if the hull is built counter-clockwise or not
    uint_t i0 = hull.safe[0], i1 = hull.safe[1], i2 = hull.safe[2];
    bool sign = (hx[i1] - hx[i0])*(hy[i2] - hy[i1]) - (hy[i1] - hy[i0])*(hx[i2] - hx[i1]) > 0;

    vec_t<Dim,bool> res = replicate(true, x.dims);
    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull.safe[i], p2 = hull.safe[i+1];
        for (uint_t p = 0; p < x.size(); ++p) {
            if (!res.safe[p]) continue;
            auto cross = (hx[p2] - hx[p1])*(y.safe[p] - hy[p1])
                       - (hy[p2] - hy[p1])*(x.safe[p] - hx[p1]);
            if ((cross < 0) == sign) res.safe[p] = false;
        }
    }

    return res;
}

// Compute the signed distance of a set of points with respect to the provided convex
// hull. Positive distances mean that the point lies inside the hull.
template<typename TX, typename TY, typename THX, typename THY, typename enable =
    typename std::enable_if<!is_vec<TX>::value && !is_vec<TY>::value>::type>
auto convex_hull_distance(const TX& x, const TY& y, const vec1u& hull,
    const THX& hx, const THY& hy) -> decltype(sqrt(x*y)) {

    phypp_check(is_hull_closed(hull, hx, hy), "the provided hull is not closed");

    // Find out if the hull is built counter-clockwise or not
    uint_t i0 = hull.safe[0], i1 = hull.safe[1], i2 = hull.safe[2];
    bool sign = (hx[i1] - hx[i0])*(hy[i2] - hy[i1]) - (hy[i1] - hy[i0])*(hx[i2] - hx[i1]) > 0;

    decltype(sqrt(x*y)) res = finf;
    bool inhull = true;

    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull.safe[i], p2 = hull.safe[i+1];

        // Get unit vector
        auto ux = hx[p2] - hx[p1];
        auto uy = hy[p2] - hy[p1];
        // Get perpendicular vector, pointing inside the hull
        auto nx = sign ? hy[p2] - hy[p1] : hy[p1] - hy[p2];
        auto ny = sign ? hx[p1] - hx[p2] : hx[p2] - hx[p1];
        // Normalize hull segment
        auto l = sqrt(sqr(nx) + sqr(ny));
        ux /= l; uy /= l; nx /= l; ny /= l;

        // Compute signed distance to current hull face line
        auto dx = x - hx[p1];
        auto dy = y - hy[p1];
        auto d = dx*nx + dy*ny;

        // Check if the point is inside the hull or not
        if (d > 0) {
            inhull = false;
        }

        d = fabs(d);
        if (res > d) {
            // Find if the projection of the point on the face lies on the segment
            auto proj = dx*ux + dy*uy;
            if (proj < 0) {
                // Projection lies before starting point
                d = sqrt(sqr(dx) + sqr(dy));
            } else if (proj > l) {
                // Projection lies after ending point
                d = sqrt(sqr(x - hx[p2]) + sqr(y - hy[p2]));
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
template<std::size_t Dim, typename TX, typename TY, typename THX, typename THY>
auto convex_hull_distance(const vec_t<Dim,TX>& x, const vec_t<Dim,TY>& y, const vec1u& hull,
    const THX& hx, const THY& hy) -> vec_t<Dim, decltype(sqrt(x[0]*y[0]))> {

    phypp_check(x.dims == y.dims, "incompatible dimensions between X and Y "
        "(", x.dims, " vs. ", y.dims, ")");
    phypp_check(is_hull_closed(hull, hx, hy), "the provided hull is not closed");

    // Find out if the hull is built counter-clockwise or not
    uint_t i0 = hull.safe[0], i1 = hull.safe[1], i2 = hull.safe[2];
    bool sign = (hx[i1] - hx[i0])*(hy[i2] - hy[i1]) - (hy[i1] - hy[i0])*(hx[i2] - hx[i1]) > 0;

    vec_t<Dim,decltype(sqrt(x[0]*y[0]))> res = replicate(finf, x.dims);
    vec_t<Dim,bool> inhull = replicate(true, x.dims);

    for (uint_t i = 0; i < hull.size()-1; ++i) {
        uint_t p1 = hull.safe[i], p2 = hull.safe[i+1];

        // Get unit vector
        auto ux = hx[p2] - hx[p1];
        auto uy = hy[p2] - hy[p1];
        // Get perpendicular vector, pointing inside the hull
        auto nx = sign ? hy[p2] - hy[p1] : hy[p1] - hy[p2];
        auto ny = sign ? hx[p1] - hx[p2] : hx[p2] - hx[p1];
        // Normalize hull segment
        auto l = sqrt(sqr(nx) + sqr(ny));
        ux /= l; uy /= l; nx /= l; ny /= l;

        for (uint_t p = 0; p < x.size(); ++p) {
            // Compute signed distance to current hull face line
            auto dx = x.safe[p] - hx[p1];
            auto dy = y.safe[p] - hy[p1];
            auto d = dx*nx + dy*ny;

            // Check if the point is inside the hull or not
            if (d > 0) {
                inhull.safe[p] = false;
            }

            d = fabs(d);
            if (res.safe[p] > d) {
                // Find if the projection of the point on the face lies on the segment
                auto proj = dx*ux + dy*uy;
                if (proj < 0) {
                    // Projection lies before starting point
                    d = sqrt(sqr(dx) + sqr(dy));
                } else if (proj > l) {
                    // Projection lies after ending point
                    d = sqrt(sqr(x.safe[p] - hx[p2]) + sqr(y.safe[p] - hy[p2]));
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
    const double d2r = 3.14159265359/180.0;
    double ra1 = d2r*tra1, ra2 = d2r*tra2, dec1 = d2r*tdec1, dec2 = d2r*tdec2;
    double sra = sin(0.5*(ra2 - ra1));
    double sde = sin(0.5*(dec2 - dec1));
    return 3600.0*2.0*asin(sqrt(sde*sde + sra*sra*cos(dec2)*cos(dec1)))/d2r;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1, typename TR2, typename TD2>
vec_t<N,double> angdist(const vec_t<N,TR1>& tra1, const vec_t<N,TD1>& tdec1,
    const vec_t<N,TR2>& tra2, const vec_t<N,TD2>& tdec2) {
    phypp_check(tra1.dims == tdec1.dims, "first RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");
    phypp_check(tra2.dims == tdec2.dims, "second RA and Dec dimensions do not match (",
        tra2.dims, " vs ", tdec2.dims, ")");
    phypp_check(tra1.dims == tra2.dims, "position sets dimensions do not match (",
        tra1.dims, " vs ", tra2.dims, ")");

    vec_t<N,double> res(tra1.dims);
    for (uint_t i : range(tra1)) {
        res.safe[i] = angdist(tra1.safe[i], tdec1.safe[i], tra2.safe[i], tdec2.safe[i]);
    }

    return res;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1>
vec_t<N,double> angdist(const vec_t<N,TR1>& tra1, const vec_t<N,TD1>& tdec1,
    double tra2, double tdec2) {
    phypp_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");

    vec_t<N,double> res(tra1.dims);
    for (uint_t i : range(tra1)) {
        res.safe[i] = angdist(tra1.safe[i], tdec1.safe[i], tra2, tdec2);
    }

    return res;
}

// Compute the angular distance between two RA/Dec positions [arcsec].
// Assumes that RA & Dec coordinates are in degrees.
template<std::size_t N, typename TR1, typename TD1>
vec_t<N,bool> angdist_less(const vec_t<N,TR1>& tra1, const vec_t<N,TD1>& tdec1,
    double tra2, double tdec2, double radius) {
    phypp_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
        tra1.dims, " vs ", tdec1.dims, ")");

    const double d2r = 3.14159265359/180.0;
    const double rrad = d2r*radius/3600.0;
    const double crad = sqr(sin(rrad/2.0));

    vec_t<N,bool> res(tra1.dims);
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
void move_ra_dec(double& ra, double& dec, double dra, double ddec) {
    ra *= cos(dec*dpi/180.0);
    dec += ddec/3600.0;
    ra += dra/3600.0;
    ra /= cos(dec*dpi/180.0);
}

// Find the closest point in a 2D array that satisfies a given criterium
bool astar_find(const vec2b& map, uint_t& x, uint_t& y) {
    phypp_check(!map.empty(), "this algorithm requires a non empty 2D vector");

    if (x >= map.dims[0]) x = map.dims[0]-1;
    if (y >= map.dims[1]) y = map.dims[1]-1;

    if (map.safe(x,y)) return true;


    using vec_pair = vec_t<1,std::pair<uint_t,uint_t>>;
    vec_pair open;
    open.push_back(std::make_pair(x,y));

    vec2b visit(map.dims);
    visit.safe(x,y) = true;

    bool found = false;
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

#endif
