#ifndef VIF_ASTRO_ASTRO_HPP
#define VIF_ASTRO_ASTRO_HPP

#include <map>
#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/utility/thread.hpp"
#include "vif/math/base.hpp"
#include "vif/math/linfit.hpp"
#include "vif/math/convex_hull.hpp"
#include "vif/math/histogram.hpp"
#include "vif/math/reduce.hpp"
#include "vif/math/random.hpp"
#include "vif/io/fits.hpp"
#include "vif/io/ascii.hpp"
#include "vif/astro/image.hpp"

namespace vif {
namespace astro {
    static std::string data_dir = system_var("VIF_DATA_DIR", "./");
    static std::string temporary_dir = system_var("VIF_TMP_DIR", "./");

    // Compute the angular distance between two RA/Dec positions [radian].
    // Assumes that RA & Dec coordinates are in radian.
    inline double angdistr(double ra1, double dec1, double ra2, double dec2) {
        double sra = sin(0.5*(ra2 - ra1));
        double sde = sin(0.5*(dec2 - dec1));
        return 2.0*asin(sqrt(sde*sde + sra*sra*cos(dec2)*cos(dec1)));
    }

    // Compute the angular distance between two RA/Dec positions [arcsec].
    // Assumes that RA & Dec coordinates are in degrees.
    inline double angdist(double tra1, double tdec1, double tra2, double tdec2) {
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
        vif_check(tra1.dims == tdec1.dims, "first RA and Dec dimensions do not match (",
            tra1.dims, " vs ", tdec1.dims, ")");
        vif_check(tra2.dims == tdec2.dims, "second RA and Dec dimensions do not match (",
            tra2.dims, " vs ", tdec2.dims, ")");
        vif_check(tra1.dims == tra2.dims, "position sets dimensions do not match (",
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
        vif_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
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
        vif_check(tra1.dims == tdec1.dims, "RA and Dec dimensions do not match (",
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
    inline void move_ra_dec(double& ra, double& dec, double dra, double ddec) {
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
        vif_check(ra.dims == dec.dims, "RA and Dec dimensions do not match (",
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
        !meta::is_vec<T>::value && !meta::is_vec<U>::value>::type>
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
        vif_check(ra.dims == dec.dims, "RA and Dec dimensions do not match (",
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
        vif_check(tr1.dims == td1.dims, "incompatible dimensions between RA and Dec "
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

    struct psffit_result {
        double bg;
        double bg_err;
        double flux;
        double flux_err;
        double chi2;
    };

    template<typename TypeM, typename TypeE = TypeM, typename TypeP = TypeM>
    psffit_result psffit(const vec<2,TypeM>& img, const vec<2,TypeE>& terr,
        const vec<2,TypeP>& psf, const vec1i& pos) {

        using img_type = typename vec<2,TypeM>::effective_type;
        using err_type = typename vec<2,TypeE>::effective_type;
        static const auto max_error = std::numeric_limits<meta::rtype_t<TypeE>>::max();

        int_t hx = psf.dims[1]/2, hy = psf.dims[0]/2;

        err_type err = terr;
        if (err.dims[0] == 1 && err.dims[1] == 1) {
            err = arr<meta::rtype_t<TypeE>>(img.dims)*0+terr[0];
        } else {
            assert(err.dims[0] == img.dims[0] && err.dims[1] == img.dims[1]);
        }

        img_type icut = subregion(img, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, 0.0);
        err_type ecut = subregion(err, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, max_error);

        vec1u idnan = where(!is_finite(icut) || !is_finite(ecut));
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
    psffit_result psffit(const vec<2,TypeM>& img, const vec<2,TypeE>& err, const vec<2,TypeP>& psf) {
        return psffit(img, err, psf, {img.dims[0]/2, img.dims[1]/2});
    }

    template<typename TypeM, typename TypeE, typename TypeP = TypeM>
    psffit_result psffit(const vec<2,TypeM>& img, const TypeE& err, const vec<2,TypeP>& psf) {
        vec<2,TypeE> terr = arr<TypeE>(1,1) + err;
        return psffit(img, terr, psf, {img.dims[0]/2, img.dims[1]/2});
    }

    template<typename TypeM, typename TypeP = TypeM>
    void srcsub(vec<2,TypeM>& img, const vec<2,TypeP>& psf, double flux, const vec1i& pos) {
        int_t hx = psf.dims[1]/2, hy = psf.dims[0]/2;

        vec1i rr, rs;
        subregion(img, {pos(0)-hx, pos(1)-hy, pos(0)+hx, pos(1)+hy}, rr, rs);
        img[rr] -= flux*psf[rs];
    }

    template<typename TypeM, typename TypeP = TypeM>
    void srcsub(vec<2,TypeM>& img, const vec<2,TypeP>& psf, double flux) {
        srcsub(img, psf, flux, {img.dims[0]/2, img.dims[1]/2});
    }

    struct cosmo_t {
        double H0 = 70.0;
        double wL = 0.73;
        double wm = 0.27;
        double wk = 0.0;
    };

    inline cosmo_t cosmo_wmap() {
        cosmo_t c;
        c.H0 = 70.0;
        c.wL = 0.73;
        c.wm = 0.27;
        c.wk = 0.0;
        return c;
    }

    inline cosmo_t cosmo_plank() {
        cosmo_t c;
        c.H0 = 67.8;
        c.wL = 0.692;
        c.wm = 0.308;
        c.wk = 0.0;
        return c;
    }

    inline cosmo_t cosmo_std() {
        cosmo_t c;
        c.H0 = 70.0;
        c.wL = 0.7;
        c.wm = 0.3;
        c.wk = 0.0;
        return c;
    }

    inline cosmo_t get_cosmo(const std::string& name) {
        if (name == "wmap") {
            return cosmo_wmap();
        } else if (name == "plank") {
            return cosmo_plank();
        } else if (name == "std") {
            return cosmo_std();
        } else {
            cosmo_t c;
            warning("unknown cosmology '"+name+"', using default (H0=", c.H0, ", omega_lambda=",
                c.wL, ", omega_matter=", c.wm, ", omega_rad=", c.wk, ")");
            return c;
        }
    }

    inline vec1s cosmo_list() {
        return {"std", "wmap", "plank"};
    }

    // Proper size of an object in [kpc/arcsec] as a function of redshift 'z'.
    // Uses lumdist() internally.
    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    T propsize(const T& z, const cosmo_t& cosmo) {
        if (z <= 0) return dinf;

        return (1.0/3.6)*(dpi/180.0)*lumdist(z, cosmo)/pow(1.0+z, 2.0);
    }

    // Luminosity distance [Mpc] as a function of redshift 'z'.
    // Note: assumes that cosmo.wk = 0.
    // There is no analytic form for this function, hence it must be numerically integrated.
    // Note that this function is costly. If you need to compute it for many redshift values, you are
    // advised to use the interpolated version below.
    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    T lumdist(const T& z, const cosmo_t& cosmo) {
        if (z <= 0) return 0.0;
        return (1.0+z)*(2.99792458e5/cosmo.H0)*integrate_func([&](T t) {
            return pow(pow((1.0+t),3)*cosmo.wm + cosmo.wL, -0.5);
        }, 0, z);
    }

    // Lookback time [Gyr] as a function of redshift 'z'.
    // There is no analytic form for this function, hence it must be numerically integrated.
    // Note that this function is costly. If you need to compute it for many redshift values, you are
    // advised to use the interpolated version below.
    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
    T lookback_time(const T& z, const cosmo_t& cosmo) {
        if (z <= 0) return 0.0;
        return (3.09/(cosmo.H0*3.155e-3))*integrate_func([&](T t) {
            return pow(pow((1+t),3)*cosmo.wm + cosmo.wL + pow(1+t,2)*cosmo.wk, -0.5)/(1+t);
        }, 0, z);
    }

    // Volume of the universe [Mpc^3] within a sphere of redshift 'z'.
    // Note: assumes that cosmo.wk = 0;
    // There is no analytic form for this function, hence it must be numerically integrated.
    // Note that this function is costly. If you need to compute it for many redshift values, you are
    // advised to use the interpolated version below.
    template<typename T, typename enable = typename std::enable_if<!meta::is_vec<T>::value>::type>
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
        typename vec<Dim,Type>::effective_type name(const vec<Dim,Type>& z, const Args& ... args) { \
            if (z.size() < 2*npt) { \
                using rtype = meta::rtype_t<Type>; \
                vec<Dim,rtype> r = z; \
                for (auto& t : r) { \
                    t = name(t, args...); \
                } \
                return r; \
            } else { \
                auto idz = where(z > 0); \
                auto mi = min(z[idz]); \
                auto ma = max(z[idz]); \
                auto tz = rgen_log(mi, ma, npt); \
                using rtype = meta::rtype_t<Type>; \
                vec<1,rtype> td = arr<rtype>(tz.dims); \
                for (uint_t i : range(tz)) { \
                    td[i] = name(tz[i], args...); \
                } \
                vec<Dim,rtype> r(z.dims); \
                r[idz] = interpolate(td, tz, z[idz]); \
                idz = where(z <= 0); \
                r[idz] = name(0.0, args...); \
                return r; \
            } \
        }

    VECTORIZE_INTERPOL(lumdist, 800);
    VECTORIZE_INTERPOL(lookback_time, 800);
    VECTORIZE_INTERPOL(propsize, 800);
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
    auto mag2lsun(const T& lam, const T& mag, double zp = 23.9) -> decltype(e10(mag)/lam) {
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

    // Flux in uJy to erg/s/cm2/A (cgs), lambda in [um]
    template<typename T, typename U>
    auto uJy2cgs(const T& lam, const U& flx) -> decltype(1.0*lam*flx) {
        const double toerg = 1.0e-29;   // 1 uJy in erg/s/cm2/Hz
        const double c = 2.9979e8;      // speed of light in m/s
        const double tometer = 1e-6;    // 1 um in meter
        const double toangstrom = 1e10; // 1 meter in Angstrom
        const double factor = toerg*c/sqr(tometer)/toangstrom; // 1 uJy to erg/s/cm2/A
        return factor*flx/lam/lam;
    }

    // Flux in erg/s/cm2/A (cgs) to uJy, lambda in [A]
    template<typename T, typename U>
    auto cgs2uJy(const T& lam, const U& flx) -> decltype(1.0*lam*flx) {
        const double toerg = 1.0e-29; // 1 uJy in erg/s/cm2/Hz
        const double c = 2.9979e8;    // speed of light in m/s
        const double tometer = 1e-10; // 1 Angstrom in meter
        const double factor = tometer/(toerg*c); // 1 uJy to erg/s/cm2/A
        return factor*flx*lam*lam;
    }

    // Compute the area covered by a field given a set of source coordinates [deg^2] and pre-computed
    // convex hull (as obtained from convex_hull() with the same coordinates).
    // Coordinates are assumed to be given in degrees.
    // NOTE: this algorithm assumes that the field is well characterized by a *convex* hull, meaning
    // that is has no hole (masked stars, ...) and the borders are not concave (no "zigzag" shape, ...).
    // If these hypotheses do not hold, the area might be overestimated. Use field_area_h2d instead.
    template<typename H>
    double field_area_hull(const convex_hull<H>& hull) {
        hull.validate();
        vif_check(hull.closed, "the provided must be closed");

        double area = 0;

        auto d2r = dpi/180.0;
        auto hx = hull.x*d2r;
        auto hy = hull.y*d2r;
        for (uint_t i : range(2, hull.size())) {
            double e1 = angdistr(hx.safe[0],   hy.safe[0],   hx.safe[i-1], hy.safe[i-1]);
            double e2 = angdistr(hx.safe[i-1], hy.safe[i-1], hx.safe[i],   hy.safe[i]);
            double e3 = angdistr(hx.safe[i],   hy.safe[i],   hx.safe[0],   hy.safe[0]);
            double p = 0.5*(e1 + e2 + e3);
            area += sqrt(p*(p-e1)*(p-e2)*(p-e3));
        }

        area /= sqr(d2r);

        return area;
    }

    // Compute the area covered by a field given a set of source coordinates [deg^2].
    // Coordinates are assumed to be given in degrees.
    template<typename TX, typename TY>
    double field_area_hull(const TX& ra, const TY& dec) {
        return field_area_hull(build_convex_hull(ra, dec));
    }

    // Compute the area covered by a field given a set of source coordinates [deg^2] by iteratively
    // building a 2d histogram of point sources.
    // Coordinates are assumed to be given in degrees.
    template<typename TX, typename TY>
    double field_area_h2d(const TX& ra, const TY& dec) {
        vif_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

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

        // Now define the grid
        uint_t inset = 2;
        uint_t nra  = ceil((max_ra  - min_ra)/dr + 2*inset);
        uint_t ndec = ceil((max_dec - min_dec)/dr + 2*inset);
        vec2d rb = make_bins(min_ra-inset*dr,  max_ra+inset*dr,  nra);
        vec2d db = make_bins(min_dec-inset*dr, max_dec+inset*dr, ndec);

        // Build the 2d histograms
        // In each cell, we compute both:
        // - the exact total area of of the cell
        vec2d carea(nra, ndec);
        // - the filling factor of the cell, i.e., what fraction is occupied by sources
        vec2d cfill(nra, ndec);

        histogram2d(ra, dec, rb, db, [&](uint_t r, uint_t d, vec1u ids) {
            if (ids.size() < 2) return;

            // Compute the area of this cell
            carea(r,d) = (angdist(rb(0,r), db(0,d), rb(1,r), db(0,d))/3600.0)*
                         (angdist(rb(0,r), db(0,d), rb(0,r), db(1,d))/3600.0);

            // Compute the average distance between sources
            const double d2r = dpi/180.0;
            auto dra  = ra[ids]*d2r;
            auto ddec = dec[ids]*d2r;
            auto dcdec = cos(ddec);

            const uint_t n = ra.size();

            auto distance_proxy = [&](uint_t i, uint_t j) {
                double sra = sin(0.5*(dra.safe[j] - dra.safe[i]));
                double sde = sin(0.5*(ddec.safe[j] - ddec.safe[i]));
                return sde*sde + sra*sra*dcdec.safe[j]*dcdec.safe[i];
            };

            auto distance = [](double td) {
                return (180.0/dpi)*2*asin(sqrt(td));
            };

            uint_t m = ids.size();
            vec1d dmin = replicate(dinf, m);
            for (uint_t i : range(m))
            for (uint_t j : range(i+1, m)) {
                double td = distance_proxy(i,j);
                if (td < dmin.safe[i]) dmin.safe[i] = td;
                if (td < dmin.safe[j]) dmin.safe[j] = td;
            }

            double dist = distance(mean(dmin));

            // Compute the filling factor
            // We make a second 2D histogram in which we only find cells where there is
            // at least one source
            const uint_t tnb = 20;
            vec2d trb = make_bins(rb(0,r), rb(1,r), tnb);
            vec2d tdb = make_bins(db(0,d), db(1,d), tnb);

            vec2b mask = histogram2d(ra[ids], dec[ids], trb, tdb) >= 1;

            // ... then inflate the corresponding mask by a radius that is equal to
            // the average minimal distance between sources, i.e., the average "radius"
            // of each source (this is not the physical size of the galaxy).
            mask = mask_inflate(mask, ceil(dist/(sqrt(carea(r,d))/tnb)));

            // The filling factor is the fraction of cells in the mask
            cfill(r,d) = fraction_of(mask);
        });

        // Find cells which are well filled
        vec2b mask = cfill > 0.5;
        // Inflate and deflate this mask to fill the holes due to statiscial fluctations
        mask = !mask_inflate(!mask_inflate(mask, 1), 1);
        // Attribute filling factor of 1 to all cells in the mask
        cfill[where(mask)] = 1.0;

        // Compute total area
        return total(cfill*carea);
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
    vec1d angcorrel(const vec<N1,TR1>& ra, const vec<N1,TD1>& dec,
        const vec<N2,TR2>& rra, const vec<N2,TD2>& rdec, const vec<2,TB>& bins) {
        vif_check(ra.dims == dec.dims, "RA and Dec dimensions do not match for the "
            "input catalog (", ra.dims, " vs ", dec.dims, ")");
        vif_check(rra.dims == rdec.dims, "RA and Dec dimensions do not match for the "
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

        double norm1 = rra.size()/double(ra.size());
        double norm2 = norm1*((rra.size() - 1.0)/(ra.size() - 1.0));
        return ((dd*norm2 - dr*norm1) + (rr - dr*norm1))/rr;
    }

    struct randpos_status {
        bool success = true;
        std::string failure;
    };

    struct randpos_uniform_options {
        uint_t max_iter = 1000;
    };


    template<typename TX, typename TY, typename F1, typename F2>
    bool rejection_sampling(vec<1,TX>& x, vec<1,TY>& y, uint_t nsrc, uint_t max_iter,
        F1&& genpos, F2&& in_region) {

        using rtx_t = meta::rtype_t<TX>;
        using rty_t = meta::rtype_t<TY>;

        x.resize(nsrc);
        y.resize(nsrc);

        for (uint_t i : range(nsrc)) {
            rtx_t tx; rty_t ty;
            genpos(tx, ty);
            uint_t iter = 0;
            while (!in_region(tx, ty)) {
                genpos(tx, ty);
                ++iter;

                if (iter == max_iter) {
                    return false;
                }
            }

            x.safe[i] = tx;
            y.safe[i] = ty;
        }

        return true;
    }

    template<typename TX, typename TY, typename TSeed, typename F>
    void randpos_uniform_box(TSeed& seed, uint_t nsrc, vec1d rx, vec1d ry,
        vec<1,TX>& x, vec<1,TY>& y) {

        if (nsrc == 0) {
            // Nothing to do...
            x.clear(); y.clear();
            return;
        }

        // Generate uniform positions in a box
        x = randomu(seed, nsrc)*(rx[1] - rx[0]) + rx[0];
        y = randomu(seed, nsrc)*(ry[1] - ry[0]) + ry[0];
    }

    template<typename TX, typename TY, typename TSeed, typename F>
    randpos_status randpos_uniform_box(TSeed& seed, uint_t nsrc, vec1d rx, vec1d ry,
        vec<1,TX>& x, vec<1,TY>& y, F&& in_region,
        randpos_uniform_options options = randpos_uniform_options{}) {

        randpos_status status;

        if (nsrc == 0) {
            // Nothing to do...
            x.clear(); y.clear();
            return status;
        }

        // Generate uniform positions in a box
        auto genpos = [&seed,rx,ry](double& tx, double& ty) mutable {
            tx = randomu(seed)*(rx[1] - rx[0]) + rx[0];
            ty = randomu(seed)*(ry[1] - ry[0]) + ry[0];
        };

        // Use rejection sampling to only fill the requested region
        bool good = rejection_sampling(x, y, nsrc, options.max_iter,
            genpos, std::forward<F>(in_region)
        );

        if (!good) {
            status.success = false;
            status.failure = "maximum number of iterations reached "
                "("+to_string(options.max_iter)+"): try increasing the max_iter value in the "
                "options or check that the provided ranges in X and Y overlap the requested "
                "region";
            return status;
        }

        return status;
    }

    template<typename TX, typename TY, typename TSeed>
    void randpos_uniform_circle(TSeed& seed, uint_t nsrc, double x0, double y0, double r0,
        vec<1,TX>& x, vec<1,TY>& y) {

        vec1d theta = randomu(seed, nsrc)*2*dpi;
        vec1d r = r0*sqrt(randomu(seed, nsrc));

        x = x0 + r*cos(theta);
        y = y0 + r*sin(theta);
    };

    template<typename TX, typename TY, typename TSeed, typename F>
    randpos_status randpos_uniform_circle(TSeed& seed, uint_t nsrc,
        double x0, double y0, double r0,
        vec<1,TX>& x, vec<1,TY>& y, F&& in_region,
        randpos_uniform_options options = randpos_uniform_options{}) {

        randpos_status status;

        if (nsrc == 0) {
            // Nothing to do...
            x.clear(); y.clear();
            return status;
        }

        // Generate uniform positions in a box
        auto genpos = [&seed,x0,y0,r0](double& tx, double& ty) mutable {
            double theta = randomu(seed)*2*dpi;
            double r = r0*sqrt(randomu(seed));

            tx = x0 + r*cos(theta);
            ty = y0 + r*sin(theta);
        };

        // Use rejection sampling to only fill the requested region
        bool good = rejection_sampling(x, y, nsrc, options.max_iter,
            genpos, std::forward<F>(in_region)
        );

        if (!good) {
            status.success = false;
            status.failure = "maximum number of iterations reached "
                "("+to_string(options.max_iter)+"): try increasing the max_iter value in the "
                "options or check that the provided circle overlaps the requested "
                "region";
            return status;
        }

        return status;
    }

    struct randpos_power_options {
        double lambda = 1.5;    // circle shrinking factor
        uint_t eta = 4;         // number of positions per circle
        uint_t eta_end = npos;  // number of positions per circle in the last level
                                // (default: same as 'eta')
        uint_t levels = 6;      // number of levels (ignored in randpos_power)
        uint_t max_iter = 1000; // maximum number of internal iterations for rejection
                                //   sampling (ignored in randpos_power_base)
    };

    // This is the Soneira & Pebbles algorithm
    template<typename TX, typename TY, typename TSeed>
    void randpos_power_base(TSeed& seed, double x0, double y0, double r0,
        vec<1,TX>& x, vec<1,TY>& y, randpos_power_options options) {

        // Initial conditions
        x = {x0};
        y = {y0};
        double r = r0;

        // Loop over the levels
        for (uint_t l : range(options.levels)) {
            // Generate positions of this level
            vec1d nx, ny;
            uint_t nprod = x.size()*options.eta;
            nx.reserve(nprod);
            ny.reserve(nprod);

            for (uint_t i : range(x)) {
                vec1d tx, ty;
                randpos_uniform_circle(seed, options.eta, x[i], y[i], r, tx, ty);
                append(nx, tx);
                append(ny, ty);
            }

            // Keep these positions for the next level
            std::swap(nx, x);
            std::swap(ny, y);

            // Shrink for next level
            r /= options.lambda;
        }
    }

    // Generate correlated positions within a circle centered at (x0,y0) and radius r0,
    // ensuring that these positions also satisfy the criterion given by in_region. The
    // algorithm assumes that r0 is not significantly larger than the extent of this region.
    // It is adapted from the Soneira & Pebbles algorithm.
    template<typename TX, typename TY, typename TSeed, typename F>
    randpos_status randpos_power_circle(TSeed& seed, double x0, double y0, double r0,
        vec<1,TX>& x, vec<1,TY>& y, F&& in_region, randpos_power_options options) {

        randpos_status status;

        // Compute the filling factor of the circle (x1,y1,r1) and its intersection with
        // the region (1: the circle lies completely inside the region, 0: the circle
        // is completely out of the region).
        auto get_fill = [&in_region, &seed](double x1, double y1, double r1,
            uint_t nfil, double threshold) mutable {

            vec1d tx, ty;
            randpos_uniform_circle(seed, nfil, x1, y1, r1, tx, ty);
            // TODO: (c++14) simplify this once generic lambdas are in the game
            vec1b inside(tx.size());
            for (uint_t i : range(tx)) {
                inside.safe[i] = in_region(tx.safe[i], ty.safe[i], threshold);
            }
            return fraction_of(inside);
        };

        // Generate uniform positions inside a circle of radius r1 at (x1,y1) that also
        // lies inside the region defined by in_region (with a tolerance margin of
        // equal to threshold, in arcsec).
        // The tolerance margin is used so that cells can be produced with a center outside
        // of the region, provided there is a sufficient filling factor for the next level.
        // This prevents holes in the position distribution close to the borders of the region
        auto genpos = [&in_region, &seed, &options](double x1, double y1, double r1,
            uint_t num, double threshold, vec1d& tx, vec1d& ty) mutable {

            return randpos_uniform_circle(seed, num, x1, y1, r1, tx, ty,
                [&](double nx, double ny) {
                    return in_region(nx, ny, threshold);
                }
            );
        };

        if (options.eta_end == npos) options.eta_end = options.eta;

        // Initialize the algorithm with random positions in the whole area
        double rth = (options.levels <= 1 ? 0.0 : 0.8*r0/options.lambda);
        uint_t eta_level = (options.levels <= 1 ? options.eta_end : options.eta);
        status = genpos(x0, y0, r0, eta_level, rth, x, y);

        if (!status.success) {
            status.failure += "\n"
                "context: placing "+to_string(options.eta)+" random positions "
                "in initial radius r="+to_string(r0)+", x="+to_string(x0)+", y="+to_string(y0)+")";
            return status;
        }

        double r1 = r0;
        for (uint_t l : range(options.levels-1)) {
            // For each position, generate new objects within a "cell" of radius 'r1'
            r1 /= options.lambda;
            rth = (l+2 == options.levels ? 0.0 : 0.8*r1/options.lambda);

            // First, compute the filling factor of each cell, i.e. what fraction of the
            // area of the cell is inside the convex hull. Each cell will then only
            // generate this fraction of objects, in order not to introduce higher densities
            // close to the borders of the field.
            vec1d fill(x.size());
            for (uint_t i : range(x)) {
                fill[i] = get_fill(x[i], y[i], r1, 200u, rth);
            }

            // Normalize filling factors so that the total number of generated object
            // is kept constant (i.e.: re-assign the missing objects to other, fully
            // filled cells)
            double mfill = mean(fill);
            vif_check(mfill != 0.0, "mean filling factor is zero, not good...");

            fill /= mfill;

            // Assign a number of object to each positions
            eta_level = (l+2 == options.levels ? options.eta_end : options.eta);
            vec1u ngal = floor(eta_level*fill);
            vec1d dngal = eta_level*fill - ngal;
            uint_t miss0 = eta_level*x.size() - total(ngal);
            uint_t miss = miss0;

            // Randomly assign fractional number of objects to make sure that the total
            // number of object is preserved.
            uint_t iter = 0;
            while (miss != 0 && iter < options.max_iter) {
                for (uint_t i : range(x)) {
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
                    "("+to_string(miss)+"/"+to_string(miss0)+") over the available positions "
                    "("+to_string(x.size())+", level="+to_string(l+1)+")";
                return status;
            }

            // Generate new positions
            vec1d x1, y1;
            for (uint_t i : range(x)) {
                if (ngal[i] == 0) continue;

                vec1d tx1, ty1;
                status = genpos(x[i], y[i], r1, ngal[i], rth, tx1, ty1);

                if (!status.success) {
                    status.failure += "\n"
                        "context: placing "+to_string(ngal[i])+" random positions "
                        "in level "+to_string(l+1)+", radius r="+to_string(r1)+", x="+to_string(x[i])+
                        ", y="+to_string(y[i])+")";
                    return status;
                }

                append(x1, tx1);
                append(y1, ty1);
            }

            // Use the newly generated positions as input for the next level
            std::swap(x, x1);
            std::swap(y, y1);
        }

        return status;
    }

    // Convert a set of sexagesimal coordinates ('hh:mm:ss.ms') into degrees
    inline bool sex2deg(const std::string& sra, const std::string& sdec, double& ra, double& dec) {
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

    template<std::size_t Dim = 1, typename TSR = std::string, typename TSD = std::string,
        typename TR, typename TD, typename enable = typename std::enable_if<
            meta::is_compatible_output_type<vec<Dim,TSR>,TR>::value &&
            meta::is_compatible_output_type<vec<Dim,TSD>,TD>::value>::type>
    vec<Dim,bool> sex2deg(const vec<Dim,TSR>& sra, const vec<Dim,TSD>& sdec, TR&& ra, TD&& dec) {
        vif_check(sra.size() == sdec.size(), "RA and Dec dimensions do not match (",
            sra.dims, " vs ", sdec.dims, ")");

        meta::resize_or_check(ra, sra.dims);
        meta::resize_or_check(dec, sra.dims);
        vec<Dim,bool> res(sra.dims);
        for (uint_t i : range(sra)) {
            res[i] = sex2deg(sra[i], sdec[i], ra[i], dec[i]);
        }

        return res;
    }

    // Convert a set of degree coordinates into sexagesimal format ('hh:mm:ss.ms')
    inline void deg2sex(double ra, double dec, std::string& sra, std::string& sdec) {
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
            std::string s = to_string(sec);
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

        sra = to_string(rah)+':'+align_right(to_string(ram),2,'0')+':'+format_sec(ras);
        sdec = to_string(dech)+':'+align_right(to_string(decm),2,'0')+':'+format_sec(decs);
    }

    template<std::size_t Dim = 1, typename TSR = double, typename TSD = double,
        typename TR, typename TD, typename enable = typename std::enable_if<
            meta::is_compatible_output_type<vec<Dim,TR>,TSR>::value &&
            meta::is_compatible_output_type<vec<Dim,TR>,TSD>::value>::type>
    void deg2sex(const vec<Dim,TR>& ra, const vec<Dim,TD>& dec, TSR&& sra, TSD&& sdec) {
        vif_check(ra.size() == dec.size(), "RA and Dec dimensions do not match (",
            ra.dims, " vs ", dec.dims, ")");

        meta::resize_or_check(sra, ra.dims);
        meta::resize_or_check(sdec, ra.dims);
        for (uint_t i : range(ra)) {
            deg2sex(ra[i], dec[i], sra[i], sdec[i]);
        }
    }

    template<std::size_t Dim, typename TR, typename TD>
    void print_radec(const std::string& file, const vec<Dim,TR>& ra, const vec<Dim,TD>& dec) {
        std::ofstream f(file);

        for (uint_t i : range(ra)) {
            f << ra[i] << "\t" << dec[i] << "\n";
        }
    }

    inline bool get_band(const vec1s& bands, const std::string& band, uint_t& bid) {
        vec1u id = where(regex_match(bands, band));
        if (id.empty()) {
            error("no band matching '"+band+"'");
            return false;
        } else if (id.size() > 1) {
            error("multiple bands matching '"+band+"'");
            for (uint_t b : id) {
                note("  candidate: band='"+bands[b]);
            }
            return false;
        }

        bid = id[0];
        return true;
    }

    template<typename Cat>
    bool get_band(const vec1s& bands, const vec1s& notes, const std::string& band,
        const std::string& note_, uint_t& bid) {
        vec1u id = where(regex_match(bands, band) && regex_match(notes, note_));
        if (id.empty()) {
            error("no band matching '"+band+"' & '"+note_+"'");
            return false;
        } else if (id.size() > 1) {
            error("multiple bands matching '"+band+"' & '"+note_+"'");
            for (uint_t b : id) {
                note("  candidate: band='"+bands[b]+"', note='"+notes[b]+"'");
            }
            return false;
        }

        bid = id[0];
        return true;
    }

}

namespace impl {
    namespace astro_impl {
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
    }
}

namespace astro {
    template<typename Cat>
    bool get_band(const Cat& cat, const std::string& band, uint_t& bid) {
        vec1s notes = impl::astro_impl::get_band_get_notes_(
            cat, impl::astro_impl::get_band_has_notes_<Cat>{}
        );

        vec1u id = where(regex_match(cat.bands, band));
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
        vec1u id = where(regex_match(cat.bands, band) && regex_match(cat.notes, note_));
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
    void pick_sources(const vec<2,Type>& img, const vec1d& x, const vec1d& y,
        int_t hsize, vec<3,meta::rtype_t<Type>>& cube, vec1u& ids) {

        uint_t width = img.dims[0];
        uint_t height = img.dims[1];
        cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*x.size());
        ids.reserve(ids.size() + x.size());

        // Loop over all sources
        for (uint_t i = 0; i < x.size(); ++i) {
            // Discard any source that falls out of the boundaries of the image
            if (x[i]-hsize < 1 || x[i]+hsize >= width || y[i] - hsize < 1 || y[i] + hsize >= height) {
                continue;
            }

            uint_t ix = round(x[i]);
            uint_t iy = round(y[i]);

            // TODO: (optimization) optimize this using subsrc?
            // TODO: fix this (will probably not compile because x and y are double)
            auto cut = img((ix-hsize)-_-(ix+hsize), (iy-hsize)-_-(iy+hsize)).concretise();

            // Discard any source that contains a bad pixel (either infinite or NaN)
            if (count(!is_finite(cut)) != 0) {
                continue;
            }

            ids.push_back(i);
            cube.push_back(cut);
        }
    }

    struct filter_t {
        double rlam = dnan; // central wavelength (micron, integrated)
        vec1d lam;          // array of wavelength positions (micron)
        vec1d res;          // array of response for each wavelength (normalized to unit area)
    };

    // Note: filters are defined so that the measured flux is always obtained by:
    //                  flux = integrate(lam, res*sed)
    // In other words, it must be converted to an "energy counter" filter (also called RSR), and
    // normalized to unit integral (i.e. integrate(lam, res) == 1).

    template<typename TypeFL, typename TypeFR>
    double sed2flux_interpolate(const vec<1,TypeFL>& l, const vec<1,TypeFR>& r, double x, uint_t i) {
        return interpolate(r.safe[i], r.safe[i+1], l.safe[i], l.safe[i+1], x);
    }

    template<typename TypeFL, typename TypeFR, typename TypeL, typename TypeS>
    double sed2flux(const vec<1,TypeFL>& flam, const vec<1,TypeFR>& fres,
        const vec<1,TypeL>& lam, const vec<1,TypeS>& sed) {

        uint_t ised = lower_bound(lam, flam.safe[0]);
        if (ised == npos) return dnan;
        uint_t ifil = 0;

        double plam = flam.safe[0];
        double pval = fres.safe[0]*sed2flux_interpolate(lam, sed, plam, ised);

        double value = 0.0;
        const uint_t nfil = flam.size();
        const uint_t nsed = lam.size();
        while (ifil < nfil-1 && ised < nsed-1) {
            double nlam, nval;

            if (flam.safe[ifil+1] < lam.safe[ised+1]) {
                // Next point is from filter
                ++ifil;
                nlam = flam.safe[ifil];
                nval = fres.safe[ifil]*sed2flux_interpolate(lam, sed, nlam, ised);
            } else {
                // Next point is from SED
                ++ised;
                nlam = lam.safe[ised];
                nval = sed.safe[ised]*sed2flux_interpolate(flam, fres, nlam, ifil);
            }

            value += 0.5*(nval + pval)*(nlam - plam);

            plam = nlam;
            pval = nval;
        }

        if (ifil != nfil - 1) return dnan;

        return value;
    }

    template<typename TypeL, typename TypeS>
    double sed2flux(const filter_t& filter, const vec<1,TypeL>& lam, const vec<1,TypeS>& sed) {
        return sed2flux(filter.lam, filter.res, lam, sed);
    }

    template<typename TypeL, typename TypeS>
    vec1d sed2flux(const filter_t& filter, const vec<2,TypeL>& lam, const vec<2,TypeS>& sed) {
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
        const vec<1,TypeL>& lam, const vec<1,TypeS>& sed) {

        vec1d rflam = lam*(1.0 + z);
        vec1d rfsed = lsun2uJy(z, d, lam, sed);
        double flx = sed2flux(from, rflam, rfsed);
        double ref = sed2flux(to, lam, sed);

        return ref/flx;
    }

    // Compute the luminosity of a rest-frame SED by integrating a given wavelength region.
    // The return value is in the same units as the input SED, i.e. most likely Lsun.
    template<typename TL, typename TS>
    double sed_luminosity(const vec<1,TL>& lam, const vec<1,TS>& sed, double l0, double l1) {
        uint_t s = lower_bound(lam, l0);
        uint_t e = upper_bound(lam, l1);
        if (s == npos || e == npos) return dnan;

        return integrate(lam[s-_-e], sed[s-_-e]/lam[s-_-e]);
    }

    // Compute the LIR luminosity of a rest-frame SED (8um to 1000um).
    template<typename TL, typename TS>
    double lir_8_1000(const TL& lam, const TS& sed) {
        return sed_luminosity(lam, sed, 8.0, 1000.0);
    }

    template <typename T>
    bool make_psf(const std::array<uint_t,2>& dims, double x0, double y0,
        const std::string& psf_model, vec<2,T>& psf) {

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
        } else if (params[0] == "file") {
            #ifndef NO_CFITSIO

            uint_t narg = 1;
            if (params.size() < narg+1) {
                error("'file' PSF model requires one argument (file name), but "
                    "none was provided");
                return false;
            } else if (params.size() > narg+1) {
                warning("'file' PSF model requires one argument (file name) but ",
                    params.size()-1, " are provided");
            }

            vec2d tpsf;
            fits::read(params[1], tpsf);

            // Trim the tPSF, make sure that the peak is at the center, and that the dimensions
            // are odds
            vec1i idm = mult_ids(tpsf, max_id(tpsf));
            int_t hsize = 0;
            int_t imax = std::min(
                std::min(idm[0], int_t(tpsf.dims[0])-1-idm[0]),
                std::min(idm[1], int_t(tpsf.dims[1])-1-idm[1])
            );

            for (int_t i = 1; i <= imax; ++i) {
                if (tpsf(idm[0]-i,idm[1]) == 0.0 &&
                    tpsf(idm[0]+i,idm[1]) == 0.0 &&
                    tpsf(idm[0],idm[1]-i) == 0.0 &&
                    tpsf(idm[0],idm[1]+i) == 0.0) {
                    hsize = i;
                    break;
                }
            }

            if (hsize == 0) hsize = imax;
            tpsf = subregion(tpsf, {idm[0]-hsize, idm[1]-hsize, idm[0]+hsize, idm[1]+hsize});

            int_t ix0 = round(x0), iy0 = round(y0);
            double dx = x0 - ix0;
            double dy = y0 - iy0;

            tpsf = translate(tpsf, dx, dy);

            psf.clear();
            psf.resize(dims);
            vec1u idi, idp;
            subregion(psf, {ix0-hsize, iy0-hsize, ix0+hsize, iy0+hsize}, idi, idp);

            psf[idi] = tpsf[idp];

            #else

            error("'file' PSF model requires the CFITSIO library");
            return false;

            #endif
        } else {
            error("unknown PSF model '", params[0], "'");
            return false;
        }

        return true;
    }

    using filter_bank_t = vec<1, filter_t>;
    using filter_db_t = std::map<std::string, std::string>;

    inline filter_db_t read_filter_db(const std::string& filename) {
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
            line = trim(line);
            if (!line.empty() && line[0] == '#') continue;

            vec1s slice = split(line, "=");
            if (slice.size() != 2) continue;

            db.insert(std::make_pair(slice[0], db_dir+slice[1]));
        }

        return db;
    }

    inline bool get_filter(const filter_db_t& db, const std::string& str, filter_t& f) {
        vec1s spl = split(str, ":");

        auto read_filter_file = [](std::string filename, filter_t& tf) {
            if (ends_with(filename, ".fits")) {
                #ifndef NO_CFITSIO
                fits::input_table(filename).read_columns(fits::narrow, ftable(tf.lam, tf.res));
                #else
                error("get_filter: cannot read FITS filter without CFITSIO library");
                return false;
                #endif
            } else {
                ascii::read_table(filename, tf.lam, tf.res);
            }

            return true;
        };

        if (spl.size() == 1) {
            auto iter = db.find(str);
            if (iter == db.end()) {
                error("get_filter: unknown filter '", str,"'");
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

                return false;
            } else {
                if (!read_filter_file(iter->second, f)) {
                    return false;
                }

                f.rlam = integrate(f.lam, f.res*f.lam);
                return true;
            }
        } else {
            if (spl[1] == "range") {
                if (spl.size() != 4) {
                    error("get_filter: 'range' filter needs 4 arguments (name, 'range', "
                        "lambda_min, lambda_max)");
                    note("got: ", str);
                    return false;
                }

                float lmin, lmax;
                if (!from_string(spl[2], lmin)) {
                    error("get_filter: '"+str+"': cannot convert lambda_min to number");
                    return false;
                }

                if (!from_string(spl[3], lmax)) {
                    error("get_filter: '"+str+"': cannot convert lambda_max to number");
                    return false;
                }

                f.lam = {lmin*0.5, lmin*0.99, lmin, lmax, lmax*1.01, lmax*2.0};
                f.res = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
                f.res /= integrate(f.lam, f.res);
                f.rlam = 0.5*(lmin + lmax);

                return true;
            } else if (spl[1] == "file") {
                if (spl.size() != 3) {
                    error("get_filter: 'file' filter needs 3 arguments (name, 'file', "
                        "path to file)");
                    note("got: ", str);
                    return false;
                }

                if (!read_filter_file(spl[2], f)) {
                    return false;
                }

                f.rlam = integrate(f.lam, f.res*f.lam);

                return true;
            } else {
                error("unknown filter type '", spl[1], "'");
                return false;
            }
        }
    }

    template<typename Type = std::string>
    bool get_filters(const filter_db_t& db, const vec<1,Type>& str, filter_bank_t& fils) {
        for (auto& s : str) {
            filter_t f;
            if (get_filter(db, s, f)) {
                fils.push_back(f);
            } else {
                return false;
            }
        }

        return true;
    }

    inline void print_filters(const filter_db_t& db) {
        for (auto& sf : db) {
            filter_t f;
            std::string filename = sf.second;
            if (ends_with(filename, ".fits")) {
                #ifndef NO_CFITSIO
                fits::input_table(filename).read_columns(fits::narrow, ftable(f.lam, f.res));
                #endif
            } else {
                ascii::read_table(filename, f.lam, f.res);
            }

            f.rlam = integrate(f.lam, f.res*f.lam);
            print(sf.first, ": ", f.rlam);
        }
    }

    using filter_map_t = std::map<std::string, uint_t>;

    inline filter_map_t read_filter_map(const std::string& filename) {
        std::ifstream file(filename);

        filter_map_t db;
        if (!file.is_open()) {
            error("read_filter_map: cannot find '"+filename+"'");
            return db;
        }

        while (!file.eof()) {
            std::string line;
            std::getline(file, line);
            line = trim(line);
            if (!line.empty() && line[0] == '#') continue;

            vec1s slice = split(line, "=");
            if (slice.size() != 2) continue;


            uint_t id;
            if (from_string(slice[1], id)) {
                db.insert(std::make_pair(slice[0], id));
            }
        }

        return db;
    }

    inline bool get_filter_id(const filter_map_t& map, const std::string& str, uint_t& id) {
        auto iter = map.find(str);
        if (iter != map.end()) {
            id = iter->second;
            return true;
        } else {
            error("get_filter_id: unknown filter '", str,"'");
            uint_t best_d = -1;
            vec1s candidates;
            for (auto& i : map) {
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
                note("get_filter_id: did you mean '", candidates[0], "'? ");
            } else if (!candidates.empty()) {
                note("get_filter_id: did you mean one of ", candidates, "?");
            }

            id = npos;
            return false;
        }
    }

    template<typename Type = std::string>
    bool get_filter_id(const filter_map_t& map, const vec<1,Type>& str, vec1u& id) {
        id.resize(str.size());
        bool good = true;
        for (uint_t i : range(str.size())) {
            if (!get_filter_id(map, str[i], id[i])) {
                good = false;
            }
        }

        return good;
    }
}
}

#endif
