#ifndef ASTRO_HPP
#define ASTRO_HPP

#include "phypp/vec.hpp"
#include "phypp/math.hpp"
#include "phypp/print.hpp"
#include "phypp/image.hpp"
#include "phypp/thread.hpp"
#include <map>

static std::string data_dir = system_var("PHYPP_DATA_DIR", "./");
static std::string temporary_dir = system_var("PHYPP_TMP_DIR", "./");

struct psffit_result {
    double bg;
    double bg_err;
    double flux;
    double flux_err;
    double chi2;
};

template<typename TypeM, typename TypeE = TypeM, typename TypeP = TypeM>
psffit_result psffit(const vec_t<2,TypeM>& img, const vec_t<2,TypeE>& terr,
    const vec_t<2,TypeP>& psf, const vec1i& pos) {

    using img_type = typename vec_t<2,TypeM>::effective_type;
    using err_type = typename vec_t<2,TypeE>::effective_type;
    static const auto max_error = std::numeric_limits<rtype_t<TypeE>>::max();

    int_t hx = psf.dims[1]/2, hy = psf.dims[0]/2;

    err_type err = terr;
    if (err.dims[0] == 1 && err.dims[1] == 1) {
        err = arr<rtype_t<TypeE>>(img.dims)*0+terr[0];
    } else {
        assert(err.dims[0] == img.dims[0] && err.dims[1] == img.dims[1]);
    }

    img_type icut = subregion(img, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, 0.0);
    err_type ecut = subregion(err, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, max_error);

    vec1u idnan = where(!finite(icut) || !finite(ecut));
    if (!idnan.empty()) {
        icut[idnan] = 0.0;
        ecut[idnan] = max_error;
    }

    auto fr = linfit(icut, ecut, 1.0, psf);
    if (fr.success) {
        return {fr.params[0], fr.errors[0], fr.params[1], fr.errors[1], fr.chi2};
    } else {
        return {dnan, dnan, dnan, dnan, fr.chi2};
    }
}

template<typename TypeM, typename TypeE = TypeM, typename TypeP = TypeM>
psffit_result psffit(const vec_t<2,TypeM>& img, const vec_t<2,TypeE>& err, const vec_t<2,TypeP>& psf) {
    return psffit(img, err, psf, {img.dims[0]/2, img.dims[1]/2});
}

template<typename TypeM, typename TypeE, typename TypeP = TypeM>
psffit_result psffit(const vec_t<2,TypeM>& img, const TypeE& err, const vec_t<2,TypeP>& psf) {
    vec_t<2,TypeE> terr = arr<TypeE>(1,1) + err;
    return psffit(img, terr, psf, {img.dims[0]/2, img.dims[1]/2});
}

template<typename TypeM, typename TypeP = TypeM>
void srcsub(vec_t<2,TypeM>& img, const vec_t<2,TypeP>& psf, double flux, const vec1i& pos) {
    int_t hx = psf.dims[1]/2, hy = psf.dims[0]/2;

    vec1i rr, rs;
    subregion(img, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, rr, rs);
    img[rr] -= flux*psf[rs];
}

template<typename TypeM, typename TypeP = TypeM>
void srcsub(vec_t<2,TypeM>& img, const vec_t<2,TypeP>& psf, double flux) {
    srcsub(img, psf, flux, {img.dims[0]/2, img.dims[1]/2});
}

struct cosmo_t {
    double H0 = 70.0;
    double wL = 0.73;
    double wm = 0.27;
    double wk = 0.0;
};

cosmo_t cosmo_wmap() {
    cosmo_t c;
    c.H0 = 70.0;
    c.wL = 0.73;
    c.wm = 0.27;
    c.wk = 0.0;
    return c;
}

cosmo_t cosmo_std() {
    cosmo_t c;
    c.H0 = 70.0;
    c.wL = 0.7;
    c.wm = 0.3;
    c.wk = 0.0;
    return c;
}

cosmo_t get_cosmo(const std::string& name) {
    if (name == "wmap") {
        return cosmo_wmap();
    } else if (name == "std") {
        return cosmo_std();
    } else {
        cosmo_t c;
        warning("unknown cosmology '"+name+"', using default (H0="+strn(c.H0)+", omega_lambda="+
            strn(c.wL)+", omega_matter="+strn(c.wm)+", omega_rad="+strn(c.wk)+")");
        return c;
    }
}

// Luminosity distance [Mpc] as a function of redshift 'z'.
// Note: assumes that cosmo.wk = 0.
// There is no analytic form for this function, hence it must be numerically integrated.
// Note that this function is costly. If you need to compute it for many redshift values, you are
// advised to use the interpolated version below.
template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T lumdist(const T& z, const cosmo_t& cosmo) {
    if (z <= 0) return 0.0;
    return (1.0+z)*(2.99792458e5/cosmo.H0)*integrate([&](T t) {
        return pow(pow((1.0+t),3)*cosmo.wm + cosmo.wL, -0.5);
    }, 0, z);
}

// Lookback time [Gyr] as a function of redshift 'z'.
// There is no analytic form for this function, hence it must be numerically integrated.
// Note that this function is costly. If you need to compute it for many redshift values, you are
// advised to use the interpolated version below.
template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T lookback_time(const T& z, const cosmo_t& cosmo) {
    if (z <= 0) return 0.0;
    return (3.09/(cosmo.H0*3.155e-3))*integrate([&](T t) {
        return pow(pow((1+t),3)*cosmo.wm + cosmo.wL + pow(1+t,2)*cosmo.wk, -0.5)/(1+t);
    }, 0, z);
}

// Volume of the universe [Mpc^3] within a sphere of redshift 'z'.
// Note: assumes that cosmo.wk = 0;
// There is no analytic form for this function, hence it must be numerically integrated.
// Note that this function is costly. If you need to compute it for many redshift values, you are
// advised to use the interpolated version below.
template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T vuniverse(const T& z, const cosmo_t& cosmo) {
    if (z <= 0) return 0.0;
    return (4.0/3.0)*dpi*pow(lumdist(z, cosmo)/(1.0+z), 3);
}

// The above functions have a vectorized version that, when trying to compute lots of elements,
// only evaluates the function on a logarithmic grid (dependent on the provided sample), and then
// interpolates the result for the full sample, effectively reducing the number of function calls.
// The quality of the approximation depends on the sparsity of the sample (sparse samples will be
// less good than compact samples), as well as the number of interpolation points, defined as an
// argument to the vectorization macro.
#define VECTORIZE_INTERPOL(name, npt) \
    template<std::size_t Dim, typename Type, typename ... Args> \
    typename vec_t<Dim,Type>::effective_type name(const vec_t<Dim,Type>& z, const Args& ... args) { \
        if (z.size() < 2*npt) { \
            using rtype = rtype_t<Type>; \
            vec_t<Dim,rtype> r = z; \
            for (auto& t : r) { \
                t = name(t, args...); \
            } \
            return r; \
        } else { \
            auto idz = where(z > 0); \
            auto mi = min(z[idz]); \
            auto ma = max(z[idz]); \
            auto tz = e10(rgen(log10(mi), log10(ma), npt)); \
            prepend(tz, {0.0}); \
            using rtype = rtype_t<Type>; \
            vec_t<1,rtype> td = arr<rtype>(tz.dims); \
            for (uint_t i : range(tz)) { \
                td[i] = name(tz[i], args...); \
            } \
            vec_t<Dim,rtype> r(z.dims); \
            r[idz] = interpolate(td, tz, z[idz]); \
            idz = where(z <= 0); \
            r[idz] = name(0.0, args...); \
            return r; \
        } \
    }

VECTORIZE_INTERPOL(lumdist, 800);
VECTORIZE_INTERPOL(lookback_time, 800);
VECTORIZE_INTERPOL(vuniverse, 800);

#undef VECTORIZE_INTERPOL

// Absolute luminosity [Lsun] to observed flux [uJy], using luminosity distance 'd' [Mpc],
// redshift 'z' [1], and rest-frame wavelength 'lam' [um]
template<typename T, typename U, typename V, typename W>
auto lsun2uJy(const T& z, const U& d, const V& lam, const W& lum) ->
    decltype((1.0 + z)*lam*lum/(d*d)) {

    const double Mpc = 3.0856e22; // [m/Mpc]
    const double Lsol = 3.839e26; // [W/Lsol]
    const double uJy = 1.0e32;    // [uJy/(W.m-2.Hz-1)]
    const double c = 2.9979e14;   // [um.s-1]
    const double factor = uJy*Lsol/(c*4.0*dpi*Mpc*Mpc);

    return factor*(1.0 + z)*lam*lum/(d*d);
}

// Observed flux 'flx' [uJy] to absolute luminosity [Lsun], using luminosity distance 'd' [Mpc],
// redshift 'z' [1], and observed wavelength 'lam' [um]
template<typename T, typename U, typename V, typename W>
auto uJy2lsun(const T& z, const U& d, const V& lam, const W& flx) ->
    decltype(1.0*flx*d*d/lam) {

    const double Mpc = 3.0856e22; // [m/Mpc]
    const double Lsol = 3.839e26; // [W/Lsol]
    const double uJy = 1.0e32;    // [uJy/(W.m-2.Hz-1)]
    const double c = 2.9979e14;   // [um.s-1]
    const double factor = c*4.0*dpi*Mpc*Mpc/(uJy*Lsol);

    return factor*flx*d*d/lam;
}

// Absolute luminosity [Lsun] to absolute magnitude [AB] using rest-frame wavelength
// 'lam' [um]
template<typename T, typename U>
auto lsun2mag(const T& lam, const U& lum, double zp = 23.9) -> decltype(1.0*lam*lum) {
    const double Mpc = 3.0856e22; // [m/Mpc]
    const double Lsol = 3.839e26; // [W/Lsol]
    const double uJy = 1.0e32;    // [uJy/(W.m-2.Hz-1)]
    const double c = 2.9979e14;   // [um.s-1]
    const double d = 10.0/1e6;    // [Mpc]
    const double factor = uJy*Lsol/(c*4.0*dpi*Mpc*Mpc*d*d);
    return -2.5*log10(factor*lam*lum) + zp;
}

// Absolute magnitude [AB] to absolute luminosity [Lsun] using rest-frame wavelength
// 'lam' [um]
template<typename T, typename U>
auto mag2lsun(const T& mag, const U& lum, double zp = 23.9) -> decltype(e10(mag)/lam) {
    const double Mpc = 3.0856e22; // [m/Mpc]
    const double Lsol = 3.839e26; // [W/Lsol]
    const double uJy = 1.0e32;    // [uJy/(W.m-2.Hz-1)]
    const double c = 2.9979e14;   // [um.s-1]
    const double d = 10.0/1e6;    // [Mpc]
    const double factor = uJy*Lsol/(c*4.0*dpi*Mpc*Mpc*d*d);
    return e10(0.4*(zp - mag))/factor/lam;
}

// Flux in uJy to AB magnitude
template<typename T>
auto uJy2mag(const T& x, double zp = 23.9) -> decltype(-2.5*log10(x) + zp) {
    return -2.5*log10(x) + zp;
}

// AB magnitude to flux in uJy
template<typename T>
auto mag2uJy(const T& x, double zp = 23.9) -> decltype(e10(0.4*(zp - x))) {
    return e10(0.4*(zp - x));
}

// Compute the area covered by a field given a set of source coordinates [deg^2] and pre-computed
// convex hull (as obtained from convex_hull() with the same coordinates).
// Coordinates are assumed to be given in degrees.
// NOTE: this algorithm assumes that the field is well characterized by a *convex* hull, meaning
// that is has no hole (masked stars, ...) and the borders are not concave (no "zigzag" shape, ...).
// If these hypotheses do not hold, use field_area_h2d instead.
template<typename TX, typename TY, typename TH>
double field_area_hull(const TH& hull, const TX& ra, const TY& dec) {
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    decltype(1.0*ra[0]*dec[0]) area = 0;

    auto d2r = dpi/180.0;
    typename TX::effective_type hx = ra[hull]*d2r;
    typename TY::effective_type hy = dec[hull]*d2r;
    uint_t nh = hull.size();
    for (uint_t i = 2; i < nh; ++i) {
        double e1 = angdistr(hx[0],   hy[0],   hx[i-1], hy[i-1]);
        double e2 = angdistr(hx[i-1], hy[i-1], hx[i],   hy[i]);
        double e3 = angdistr(hx[i],   hy[i],   hx[0],   hy[0]);
        double p = 0.5*(e1 + e2 + e3);
        area += sqrt(p*(p-e1)*(p-e2)*(p-e3));
    }

    area /= d2r*d2r;

    return area;
}

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
double field_area_hull(const TX& ra, const TY& dec) {
    return field_area_hull(convex_hull(ra, dec), ra, dec);
}

// Compute the area covered by a field given a set of source coordinates [deg^2] by iteratively
// building a 2d histogram of point sources.
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
double field_area_h2d(const TX& ra, const TY& dec) {
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    // Optimal bin size
    // If field is perfectly uniform, then all cells of the field will contain 'optinum' objects
    const uint_t optinum = 10;

    if (ra.size() < 2*optinum) {
        // Too few sources, fallback to hull based approach
        return field_area_hull(ra, dec);
    }

    // Get rough extents of the field
    double min_ra  = min(ra),  max_ra  = max(ra);
    double min_dec = min(dec), max_dec = max(dec);
    double rough_area = (angdist(min_ra, min_dec, max_ra, min_dec)/3600.0)*
                        (angdist(min_ra, min_dec, min_ra, max_dec)/3600.0);

    double dr = sqrt(optinum)*sqrt(rough_area/ra.size());

    // Since the field is possibly not uniform, we have to introduce a tolerence threshold.
    uint_t threshold = 0.5*optinum;

    // Now define the grid
    uint_t nra  = ceil((max_ra  - min_ra)/dr);
    uint_t ndec = ceil((max_dec - min_dec)/dr);
    vec2d rb = make_bins(min_ra,  max_ra,  nra);
    vec2d db = make_bins(min_dec, max_dec, ndec);

    // Compute the area occupied by a cell
    double cell_area = (rb(1,0) - rb(0,0))*(db(1,0) - db(0,0));

    // Build the 2d histogram
    vec2u counts = histogram2d(ra, dec, rb, db);

    // Sum up the areas of all cells with counts above the threshold
    return cell_area*total(counts > threshold);
}

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
double field_area(const TX& ra, const TY& dec) {
    return field_area_h2d(ra, dec);
}

// Compute 2 point angular correlation function of a data set with positions 'ra' and 'dec'
// against a set of random positions uniformly drawn in the same region of space 'rra' and
// 'rdec'. For good results, there must be at least as many random positions as there are
// input positions, and results get better the more random positions are given.
// Compute the correlation in given bins of angular separation (in arcseconds).
// Uses the Landy-Szalay estimator, computed with a brute force approach.
template<std::size_t N1, typename TR1, typename TD1,
    std::size_t N2, typename TR2, typename TD2, typename TB>
vec1d angcorrel(const vec_t<N1,TR1>& ra, const vec_t<N1,TD1>& dec,
    const vec_t<N2,TR2>& rra, const vec_t<N2,TD2>& rdec, const vec_t<2,TB>& bins) {
    phypp_check(ra.dims == dec.dims, "RA and Dec dimensions do not match for the "
        "input catalog (", ra.dims, " vs ", dec.dims, ")");
    phypp_check(rra.dims == rdec.dims, "RA and Dec dimensions do not match for the "
        "random catalog (", rra.dims, " vs ", rdec.dims, ")");

    uint_t nbin = bins.dims[1];

    vec1d dd(nbin);
    vec1d dr(nbin);
    vec1d rr(nbin);

    for (uint_t i : range(ra)) {
        vec1d d = angdist(ra, dec, ra[i], dec[i]);
        dd += histogram(d, bins);

        d = angdist(rra, rdec, ra[i], dec[i]);
        dr += histogram(d, bins);
    }

    for (uint_t i : range(rra)) {
        vec1d d = angdist(rra, rdec, rra[i], rdec[i]);
        rr += histogram(d, bins);
    }

    vec1d dbin = bins(1,_) - bins(0,_);
    double norm = rra.size()/double(ra.size());
    return ((dd*sqr(norm) - dr*norm) + (rr - dr*norm))/(rr*dbin);
}

struct randpos_status {
    bool success = true;
    std::string failure;
};

struct randpos_uniform_options {
    uint_t nsrc = 1000;
    uint_t max_iter = 1000;
};

template<std::size_t N1, typename TR1, typename TD1, typename TSeed, typename F>
randpos_status randpos_uniform(TSeed& seed, vec1d rra, vec1d rdec, F&& in_region,
    vec_t<N1,TR1>& ra, vec_t<N1,TD1>& dec, randpos_uniform_options options) {

    randpos_status status;

    if (options.nsrc == 0) {
        // Nothing to do...
        ra.clear(); dec.clear();
        return status;
    }

    ra = (randomu(seed, options.nsrc) - 0.5)*(rra[1] - rra[0]) + mean(rra);
    dec = (randomu(seed, options.nsrc) - 0.5)*(rdec[1] - rdec[0]) + mean(rdec);

    vec1u bid = where(!in_region(ra, dec));
    uint_t iter = 0;
    while (!bid.empty() && iter < options.max_iter) {
        ra[bid] = (randomu(seed, bid.size()) - 0.5)*(rra[1] - rra[0]) + mean(rra);
        dec[bid] = (randomu(seed, bid.size()) - 0.5)*(rdec[1] - rdec[0]) + mean(rdec);
        bid = bid[where(!in_region(ra[bid], dec[bid]))];
        ++iter;
    }

    if (iter == options.max_iter) {
        status.success = false;
        status.failure = "maximum number of iterations reached "
            "("+strn(options.max_iter)+"): try increasing the max_iter value or "
            "check that the provided ranges in RA and Dec overlap the requested region";
        return status;
    }

    return status;
}

struct randpos_power_options {
    double power = 1.0;
    uint_t levels = 5;
    double inflate = 1.0;
    double overload = 1.0;
    uint_t nsrc = 1000;
    uint_t max_iter = 1000;
};

template<std::size_t N1, typename TR1, typename TD1, typename TSeed>
randpos_status randpos_power_circle(TSeed& seed, double ra0, double dec0, double r0,
    vec_t<N1,TR1>& ra, vec_t<N1,TD1>& dec, randpos_power_options options) {

    randpos_status status;

    if (options.nsrc == 0) {
        // Nothing to do...
        ra.clear(); dec.clear();
        return status;
    }

    // If asked, inflate the initial generation radius
    double rr0 = options.inflate*r0;

    // Choose the number of objects to generate
    uint_t nsim = options.overload*options.nsrc;
    // Increase it to compensate area inflation
    double inflate_fudge_factor = 1.1;
    nsim = ceil(nsim*sqr(options.inflate)*inflate_fudge_factor);

    // If the wanted number of object is too low, the algorithm will not work
    // well, so we generate more and will randomly pick among those afterwards.
    if (nsim < 1000) nsim = 1000;

    // The number of object per cell
    uint_t eta = ceil(pow(nsim, 1.0/options.levels));
    if (eta <= 1) {
        eta = 2;
        options.levels = ceil(log10(nsim)/log10(eta));
    }

    // The radius shrinking factor
    double lambda = pow(eta, 1.0/(2.0 - options.power));

    auto in_circle = [&ra0, &dec0, rr0](vec1d ra, vec1d dec, double threshold) {
        return angdist_less(ra, dec, ra0, dec0, (rr0+threshold)*3600.0);
    };

    auto get_fill = [&in_circle, &seed](double ra1, double dec1,
        double r1, uint_t num, double threshold) mutable {

        vec1d tdec = (randomu(seed, num) - 0.5)*2*r1 + dec1;
        vec1d tra  = (randomu(seed, num) - 0.5)*2*r1/cos(tdec*dpi/180.0) + ra1;

        vec1u id = where(angdist_less(tra, tdec, ra1, dec1, r1*3600.0));

        return fraction_of(in_circle(tra[id], tdec[id], threshold));
    };

    auto gen = [&in_circle, &seed, &options](double ra1, double dec1,
        double r1, uint_t num, double threshold, vec1d& tra, vec1d& tdec) mutable {

        tdec = (randomu(seed, num) - 0.5)*2*r1 + dec1;
        tra  = (randomu(seed, num) - 0.5)*2*r1/cos(tdec*dpi/180.0) + ra1;

        vec1u bid = where(!in_circle(tra, tdec, threshold) ||
            !angdist_less(tra, tdec, ra1, dec1, r1*3600.0));

        uint_t iter = 0;
        while (!bid.empty() && iter < options.max_iter) {
            tdec[bid] = (randomu(seed, bid.size()) - 0.5)*2*r1 + dec1;
            tra[bid]  = (randomu(seed, bid.size()) - 0.5)*2*r1/cos(tdec[bid]*dpi/180.0) + ra1;
            bid = bid[where(!in_circle(tra[bid], tdec[bid], threshold) ||
                !angdist_less(tra[bid], tdec[bid], ra1, dec1, r1*3600.0))];
            ++iter;
        }

        return bid.size();
    };

    // Initialize the algorithm with random positions in the whole area
    uint_t bad = gen(ra0, dec0, rr0, eta, (options.levels <= 1 ? 0.0 : rr0), ra, dec);
    if (bad != 0) {
        status.success = false;
        status.failure = "could not place all random positions "
            "("+strn(bad)+"/"+strn(eta)+") in initial radius "
            "("+strn(rr0*3600)+"\", "+strn(ra0)+", "+strn(dec0)+")";
        return status;
    }

    for (int_t l : range(options.levels-1)) {
        // For each position, generate new objects within a "cell" of radius 'r1'
        double r1 = rr0*pow(lambda, -l-1);
        double rth = (l == options.levels-2 ? 0.0 : 0.8*r1/lambda);

        // First, compute the filling factor of each cell, i.e. what fraction of the
        // area of the cell is inside the convex hull. Each cell will then only
        // generate this fraction of objects, in order not to introduce higher densities
        // close to the borders of the field.
        vec1d fill(ra.size());
        for (uint_t i : range(ra)) {
            fill[i] = get_fill(ra[i], dec[i], r1, 100u, rth);
        }

        // Normalize filling factors so that the total number of generated object
        // is kept constant.
        fill /= mean(fill);

        // Assign a number of objects to each positions
        vec1u ngal = floor(eta*fill);
        vec1d dngal = eta*fill - ngal;
        uint_t miss0 = pow(eta, l+2) - total(ngal);
        uint_t miss = miss0;

        // Randomly assign fractional number of objects to make sure that the total
        // number of object is preserved.
        uint_t iter = 0;
        while (miss != 0 && iter < options.max_iter) {
            for (uint_t i : range(ra)) {
                if (dngal[i] != 0.0 && randomu(seed) < dngal[i]) {
                    ++ngal[i];
                    dngal[i] = 0.0;

                    --miss;
                    if (miss == 0) break;
                }
            }

            ++iter;
        }

        if (miss != 0) {
            status.success = false;
            status.failure = "could not reassign some fractional random positions "
                "("+strn(miss)+"/"+strn(miss0)+") over the available positions "
                "("+strn(ra.size())+", level="+strn(l+1)+")";
            return status;
        }

        // Generate new positions
        vec1d ra1, dec1;
        for (uint_t i : range(ra)) {
            if (ngal[i] == 0) continue;

            vec1d tra1, tdec1;
            bad = gen(ra[i], dec[i], r1, ngal[i], rth, tra1, tdec1);
            if (bad != 0) {
                status.success = false;
                status.failure = "could not place all random positions "
                    "("+strn(bad)+"/"+strn(eta)+") in level "+strn(l+1)+" "
                    "radius ("+strn(r1*3600)+"\", "+strn(ra[i])+", "+strn(dec[i])+") "
                    "and requested circle ("+strn(rr0*3600)+"\", "+strn(ra0)+", "+strn(dec0)+")";
                return status;
            }

            append(ra1, tra1);
            append(dec1, tdec1);
        }

        // Use the newly generated positions as input for the next level
        std::swap(ra, ra1);
        std::swap(dec, dec1);
    }

    // Remove sources outside of asked radius
    vec1u idi = where(angdist_less(ra, dec, ra0, dec0, r0*3600.0));
    if (idi.size() < options.nsrc) {
        status.success = false;
        status.failure = "too few generated positions end up in the requested circle "
            "("+strn(idi.size())+" vs "+strn(options.nsrc)+"): try increasing the "
            "'overload'";
        return status;
    }

    // Trim catalog to match the input number of sources
    vec1u fids = shuffle(idi, seed)[uindgen(options.nsrc)];
    ra = ra[fids];
    dec = dec[fids];

    return status;
}

template<std::size_t N1, typename TR1, typename TD1,
    std::size_t N2, typename TR2, typename TD2, typename TSeed>
randpos_status randpos_power(TSeed& seed,
    const vec1u& hull, const vec_t<N2,TR2>& hra, const vec_t<N2,TD2>& hdec,
    vec_t<N1,TR1>& ra, vec_t<N1,TD1>& dec, randpos_power_options options) {

    randpos_status status;

    if (options.nsrc == 0) {
        // Nothing to do...
        ra.clear(); dec.clear();
        return status;
    }

    // Compute bounding box of the provided hull
    vec1d rra = {min(hra[hull]), max(hra[hull])};
    vec1d rdec = {min(hdec[hull]), max(hdec[hull])};

    // The initial radius is taken as the maximum distance from the field center.
    double r0 = max(angdist(hra[hull], hdec[hull], mean(rra), mean(rdec)))/3600.0;

    // Compute the ratio between the area of the generation circle and the final
    // request convex hull to keep the source density constant
    double filling_factor = dpi*sqr(r0)/field_area_hull(hull, hra, hdec);

    // Add some overload to make sure that enough objects are created
    double circle_overload = 4.0;

    auto circle_options = options;
    circle_options.nsrc *= circle_overload*filling_factor*sqr(options.inflate);
    if (circle_options.nsrc < 3000) circle_options.nsrc = 3000;
    circle_options.inflate = 1.0;

    status = randpos_power_circle(
        seed, mean(rra), mean(rdec), options.inflate*r0, ra, dec, circle_options
    );

    if (!status.success) {
        return status;
    }

    vec1u idi = where(in_convex_hull(ra, dec, hull, hra, hdec));
    if (idi.size() < options.nsrc) {
        status.success = false;
        status.failure = "too few generated positions end up in the requested hull "
            "("+strn(idi.size())+" vs "+strn(options.nsrc)+"): try increasing the "
            "'overload'";
        return status;
    }

    vec1u fids = shuffle(idi, seed)[uindgen(options.nsrc)];
    ra = ra[fids];
    dec = dec[fids];

    return status;
}

// Convert a set of sexagesimal coordinates ('hh:mm:ss.ms') into degrees
bool sex2deg(const std::string& sra, const std::string& sdec, double& ra, double& dec) {
    vec1s starr1 = split(sra, ":");
    vec1s starr2 = split(sdec, ":");

    if (starr1.size() != 3 || starr2.size() != 3) {
        ra = dnan;
        dec = dnan;
        return false;
    }

    int_t rah, ram, dech, decm;
    double ras, decs;

    if (!from_string(starr1[0], rah)  || !from_string(starr1[1], ram) ||
        !from_string(starr1[2], ras)  || !from_string(starr2[0], dech) ||
        !from_string(starr2[1], decm) || !from_string(starr2[2], decs)) {
        ra = dnan;
        dec = dnan;
        return false;
    }

    double signr = sign(rah);
    double signd = sign(dech);

    ra = (rah + ram*signr/60.0 + ras*signr/3600.0)*15.0;
    dec = dech + decm*signd/60.0 + decs*signd/3600.0;

    return true;
}

template<std::size_t Dim, typename TSR, typename TSD, typename TR, typename TD>
vec_t<Dim,bool> sex2deg(const vec_t<Dim,TSR>& sra, const vec_t<Dim,TSD>& sdec,
    vec_t<Dim,TR>& ra, vec_t<Dim,TD>& dec) {

    phypp_check(sra.size() == sdec.size(), "RA and Dec dimensions do not match (",
        sra.dims, " vs ", sdec.dims, ")");

    ra.resize(sra.dims);
    dec.resize(sra.dims);
    vec_t<Dim,bool> res(sra.dims);
    for (uint_t i : range(sra)) {
        res[i] = sex2deg(sra[i], sdec[i], ra[i], dec[i]);
    }

    return res;
}

// Convert a set of degree coordinates into sexagesimal format ('hh:mm:ss.ms')
void deg2sex(double ra, double dec, std::string& sra, std::string& sdec) {
    double signr = sign(ra);
    ra /= 15.0*signr;
    int_t  rah = ra;
    int_t  ram = (ra - rah)*60.0;
    double ras = ((ra - rah)*60 - ram)*60;
    rah *= signr;

    double signd = sign(dec);
    dec *= signd;
    int_t  dech = dec;
    int_t  decm = (dec - dech)*60.0;
    double decs = ((dec - dech)*60 - decm)*60;
    dech *= signd;

    auto format_sec = [](double sec) {
        std::string s = strn(sec);
        auto p = s.find_first_of('.');
        if (p == s.npos) {
            if (s.size() != 2) {
                return "0" + s + ".0";
            } else {
                return s + ".0";
            }
        } else {
            if (p != 2u) {
                return "0" + s;
            } else {
                return s;
            }
        }
    };

    sra = strn(rah)+':'+strn(ram,2)+':'+format_sec(ras);
    sdec = strn(dech)+':'+strn(decm,2)+':'+format_sec(decs);
}

template<std::size_t Dim, typename TSR, typename TSD, typename TR, typename TD>
void deg2sex(const vec_t<Dim,TR>& ra, const vec_t<Dim,TD>& dec, vec_t<Dim,TSR>& sra, vec_t<Dim,TSD>& sdec) {
    phypp_check(ra.size() == dec.size(), "RA and Dec dimensions do not match (",
        ra.dims, " vs ", dec.dims, ")");

    sra.resize(ra.dims);
    sdec.resize(ra.dims);
    for (uint_t i : range(ra)) {
        deg2sex(ra[i], dec[i], sra[i], sdec[i]);
    }
}

template<std::size_t Dim, typename TR, typename TD>
void print_radec(const std::string& file, const vec_t<Dim,TR>& ra, const vec_t<Dim,TD>& dec) {
    std::ofstream f(file);

    for (uint_t i : range(ra)) {
        f << ra[i] << "\t" << dec[i] << "\n";
    }
}

template<typename Cat>
vec1s get_band_get_notes_(const Cat& cat, std::false_type) {
    return replicate("?", cat.bands.size());
}

template<typename Cat>
vec1s get_band_get_notes_(const Cat& cat, std::true_type) {
    return cat.notes;
}

template<typename Cat>
struct get_band_has_notes_impl_ {
    template <typename U> static std::true_type dummy(typename std::decay<
        decltype(std::declval<U&>().notes)>::type*);
    template <typename U> static std::false_type dummy(...);
    using type = decltype(dummy<Cat>(0));
};

template<typename Cat>
using get_band_has_notes_ = typename get_band_has_notes_impl_<typename std::decay<Cat>::type>::type;

template<typename Cat>
bool get_band(const Cat& cat, const std::string& band, uint_t& bid) {
    vec1s notes = get_band_get_notes_(cat, get_band_has_notes_<Cat>{});
    vec1u id = where(match(cat.bands, band));
    if (id.empty()) {
        error("no band matching '"+band+"'");
        return false;
    } else if (id.size() > 1) {
        error("multiple bands matching '"+band+"'");
        for (uint_t b : id) {
            note("  candidate: band='"+cat.bands[b]+"', note='"+notes[b]+"'");
        }
        return false;
    }

    bid = id[0];
    return true;
}

template<typename Cat>
bool get_band(const Cat& cat, const std::string& band, const std::string& note_, uint_t& bid) {
    vec1u id = where(match(cat.bands, band) && match(cat.notes, note_));
    if (id.empty()) {
        error("no band matching '"+band+"' & '"+note_+"'");
        return false;
    } else if (id.size() > 1) {
        error("multiple bands matching '"+band+"' & '"+note_+"'");
        for (uint_t b : id) {
            note("  candidate: band='"+cat.bands[b]+"', note='"+cat.notes[b]+"'");
        }
        return false;
    }

    bid = id[0];
    return true;
}

template<typename Type>
void pick_sources(const vec_t<2,Type>& img, const vec1d& x, const vec1d& y,
    int_t hsize, vec_t<3,rtype_t<Type>>& cube, vec1u& ids) {

    uint_t width = img.dims[0];
    uint_t height = img.dims[1];
    cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*x.size());
    ids.reserve(ids.size() + x.size());

    // Loop over all sources
    vec1i r = rgen(-hsize, hsize);
    for (uint_t i = 0; i < x.size(); ++i) {
        // Discard any source that falls out of the boundaries of the image
        if (x[i]-hsize < 1 || x[i]+hsize >= width || y[i] - hsize < 1 || y[i] + hsize >= height) {
            continue;
        }

        // TODO: optimize this using subsrc?
        // TODO: fix this (will probably not compile)
        auto cut = img(x[i]+r, y[i]+r).concretise();

        // Discard any source that contains a bad pixel (either infinite or NaN)
        if (total(!finite(cut)) != 0) {
            continue;
        }

        ids.push_back(i);
        cube.push_back(cut);
    }
}

struct filter_t {
    double rlam = dnan; // central wavelenght (micron, integrated)
    vec1d lam;          // array of wavelength positions (micron)
    vec1d res;          // array of response for each wavelength (normalized to unit area)
};

// Note: filters are defined so that the measured flux is always obtained by:
//                  flux = integrate(lam, res*sed)
// In other words, it must be converted to an "energy counter" filter (also called RSR), and
// normalized to unit integral (i.e. integrate(lam, res) == 1).

template<typename TypeL, typename TypeS>
double sed2flux(const filter_t& filter, const vec_t<1,TypeL>& lam, const vec_t<1,TypeS>& sed) {
    uint_t nflam = filter.lam.size();

    auto bnd = bounds(filter.lam.safe[0], filter.lam.safe[nflam-1], lam);
    if (bnd[0] == npos || bnd[1] == npos) {
        return dnan;
    }

    vec1d nlam; nlam.reserve(nflam + bnd[1] - bnd[0] + 1);
    vec1d nrs;  nrs.reserve( nflam + bnd[1] - bnd[0] + 1);

    uint_t j = bnd[0];
    for (uint_t i : range(filter.lam)) {
        nlam.push_back(filter.lam.safe[i]);
        nrs.push_back(filter.res.safe[i]*interpolate(
            sed.safe[j-1], sed.safe[j], lam.safe[j-1], lam.safe[j], filter.lam.safe[i]));

        if (i != nflam - 1) {
            while (lam[j] < filter.lam.safe[i+1]) {
                nlam.push_back(lam.safe[j]);
                nrs.push_back(sed.safe[j]*interpolate(
                    filter.res.safe[i], filter.res.safe[i+1],
                    filter.lam.safe[i], filter.lam.safe[i+1], lam.safe[j]
                ));
                ++j;
            }
        }
    }

    return integrate(nlam, nrs);
}

template<typename TypeL, typename TypeS>
vec1d sed2flux(const filter_t& filter, const vec_t<2,TypeL>& lam, const vec_t<2,TypeS>& sed) {
    vec1d r;
    const uint_t nsed = sed.dims[0];
    r.reserve(nsed);

    for (uint_t s : range(nsed)) {
        r.push_back(sed2flux(filter, lam.safe(s,_).concretise(), sed.safe(s,_).concretise()));
    }

    return r;
}

template<typename TypeL, typename TypeS>
double sed_convert(const filter_t& from, const filter_t& to, double z, double d,
    const vec_t<1,TypeL>& lam, const vec_t<1,TypeS>& sed) {

    vec1d rflam = lam*(1.0 + z);
    vec1d rfsed = lsun2uJy(z, d, lam, sed);
    double flx = sed2flux(from, rflam, rfsed);
    double ref = sed2flux(to, lam, sed);

    return ref/flx;
}

// Compute the LIR luminosity of a rest-frame SED (8um to 1000um).
template<typename TL, typename TS>
double lir_8_1000(const vec_t<1,TL>& lam, const vec_t<1,TS>& sed) {
    uint_t s = lower_bound(8.0, lam);
    uint_t e = upper_bound(1000.0, lam);
    if (s == npos || e == npos) return dnan;

    return integrate(lam[s-_-e], sed[s-_-e]/lam[s-_-e]);
}

bool make_psf(const std::array<uint_t,2>& dims, double x0, double y0,
    const std::string& psf_model, vec2d& psf) {

    vec1s params = trim(split(psf_model, ","));
    if (params[0] == "gaussian") {
        uint_t narg = 2;
        if (params.size() < narg+1) {
            error("'gaussian' PSF model requires ", narg, " arguments (peak amplitude "
                "and width in pixels), but only ", params.size()-1, " are provided");
            return false;
        } else if (params.size() > narg+1) {
            warning("'gaussian' PSF model requires ", narg, " arguments (peak amplitude "
                "and width in pixels), but ", params.size()-1, " are provided");
        }

        double amp, width;
        from_string(params[1], amp);
        from_string(params[2], width);

        psf = generate_img(dims, [=](uint_t x, uint_t y) {
            return amp*exp(-(sqr(double(x) - x0) + sqr(double(y) - y0))/(2*sqr(width)));
        });
    } else if (params[0] == "gaussring") {
        uint_t narg = 5;
        if (params.size() < narg+1) {
            error("'gaussring' PSF model requires ", narg, " arguments (peak amplitude "
                ", peak width in pixels, ring distance in pixels, ring fraction and ring width "
                "in pixels), but only ", params.size()-1, " are provided");
            return false;
        } else if (params.size() > narg+1) {
            warning("'gaussring' PSF model requires ", narg, " arguments (peak amplitude "
                ", peak width in pixels, ring distance in pixels, ring fraction and ring width "
                "in pixels), but ", params.size()-1, " are provided");
        }

        double amp, width, rpos, rfrac, rwidth;
        from_string(params[1], amp);
        from_string(params[2], width);
        from_string(params[3], rpos);
        from_string(params[4], rfrac);
        from_string(params[5], rwidth);

        psf = generate_img(dims, [=](uint_t x, uint_t y) {
            double r2 = sqr(double(x) - x0) + sqr(double(y) - y0);
            return amp*(
                exp(-r2/(2*sqr(width)))
                + rfrac*exp(-sqr(sqrt(r2) - rpos)/(2.0*sqr(rwidth)))
            );
        });
    } else {
        error("unknown PSF model '", params[0], "'");
        return false;
    }

    return true;
}

using filter_bank_t = vec_t<1, filter_t>;
using filter_db_t = std::map<std::string, std::string>;

filter_db_t read_filter_db(const std::string& filename) {
    std::ifstream file(filename);

    filter_db_t db;
    if (!file.is_open()) {
        error("read_filter_db: cannot find '"+filename+"'");
        return db;
    }

    std::string db_dir = file::get_directory(filename);

    while (!file.eof()) {
        std::string line;
        std::getline(file, line);
        vec1s slice = split(line, "=");
        if (slice.size() != 2) continue;

        db.insert(std::make_pair(slice[0], db_dir+slice[1]));
    }

    return db;
}

filter_t get_filter(const filter_db_t& db, const std::string& str) {
    auto iter = db.find(str);
    if (iter == db.end()) {
        warning("get_filter: unknown filter '", str,"'");
        uint_t best_d = -1;
        vec1s candidates;
        for (auto& i : db) {
            uint_t d = distance(str, i.first);
            if (d < best_d) {
                best_d = d;
                candidates.data.clear();
                candidates.data.push_back(i.first);
            } else if (d == best_d) {
                candidates.data.push_back(i.first);
            }
        }

        if (candidates.size() == 1) {
            note("get_filter: did you mean '", candidates[0], "'? ");
        } else if (!candidates.empty()) {
            note("get_filter: did you mean one of ", candidates, "?");
        }

        return filter_t();
    } else {
        filter_t f;
        fits::read_table(iter->second, ftable(f.lam, f.res));
        f.rlam = integrate(f.lam, f.res*f.lam);
        return f;
    }
}

template<typename Type = std::string>
filter_bank_t get_filters(const filter_db_t& db, const vec_t<1,Type>& str) {
    filter_bank_t fils;
    for (auto& s : str) {
        fils.push_back(get_filter(db, s));
    }

    return fils;
}

void print_filters(const filter_db_t& db) {
    for (auto& sf : db) {
        filter_t f;
        fits::read_table(sf.second, ftable(f.lam, f.res));
        f.rlam = integrate(f.lam, f.res*f.lam);
        print(sf.first, ": ", f.rlam);
    }
}

#endif
