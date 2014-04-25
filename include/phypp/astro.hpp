#ifndef ASTRO_HPP
#define ASTRO_HPP

#include "phypp/vec.hpp"
#include "phypp/math.hpp"
#include "phypp/print.hpp"
#include "phypp/image.hpp"
#include "phypp/thread.hpp"
#include <map>

static std::string data_dir = system_var("PHYPP_DATA_DIR", "./");

struct psffit_result {
    double bg;
    double bg_err;
    double flux;
    double flux_err;
    double chi2;
};

template<typename TypeM, typename TypeE = TypeM, typename TypeP = TypeM>
psffit_result psffit(const vec_t<2,TypeM>& img, const vec_t<2,TypeE>& terr, const vec_t<2,TypeP>& psf, const vec1i& pos) {
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
        if (n_elements(z) < 2*npt) { \
            using rtype = rtype_t<Type>; \
            vec_t<1,rtype> r = z; \
            for (auto& t : r) { \
                t = name(t, args...); \
            } \
            return r; \
        } else { \
            auto idz = where(z > 0); \
            auto mi = min(z[idz]); \
            auto ma = max(z[idz]); \
            auto tz = merge(0.0, e10(rgen(log10(mi), log10(ma), npt))); \
            using rtype = rtype_t<Type>; \
            vec_t<1,rtype> td = arr<rtype>(npt+1); \
            for (uint_t i = 0; i < npt+1; ++i) { \
                td[i] = name(tz[i], args...); \
            } \
            vec_t<1,rtype> r = z; \
            for (auto& t : r) { \
                t = interpolate(td, tz, t); \
            } \
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
auto lsun2uJy(const T& z, const U& d, const V& lam, const W& lum) {
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
auto uJy2lsun(const T& z, const U& d, const V& lam, const W& flx) {
    const double Mpc = 3.0856e22; // [m/Mpc]
    const double Lsol = 3.839e26; // [W/Lsol]
    const double uJy = 1.0e32;    // [uJy/(W.m-2.Hz-1)]
    const double c = 2.9979e14;   // [um.s-1]
    const double factor = c*4.0*dpi*Mpc*Mpc/(uJy*Lsol);

    return factor*flx*d*d/lam;
}

// Flux in uJy to AB magnitude
template<typename T>
auto uJy2mag(const T& x, double zp = 23.9) {
    return -2.5*log10(x) + zp;
}

// AB magnitude to flux in uJy
template<typename T>
auto mag2uJy(const T& x, double zp = 23.9) {
    return e10(0.4*(zp - x));
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
        res[i] = angdist(tra1[i], tdec1[i], tra2[i], tdec2[i]);
    }

    return res;
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
    int_t rah, ram, dech, decm;
    double ras, decs;

    double signr = sign(ra);
    ra /= 15.0*signr;
    rah = ra;
    ram = (ra - rah)*60.0;
    ras = ((ra - rah)*60 - ram)*60;
    rah *= signr;

    double signd = sign(dec);
    dec *= signd;
    dech = dec;
    decm = (dec - dech)*60.0;
    decs = ((dec - dech)*60 - decm)*60;
    dech *= signd;

    auto format_sec = [](double sec) {
        std::string s = strn(sec);
        auto p = s.find_first_of('.');
        if (p == s.npos) {
            if (s.size() != 2) {
                return "0"+ s + ".0";
            } else {
                return s + ".0";
            }
        } else {
            if (p != 2u) {
                return "0"+ s;
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

struct qxmatch_res {
    vec2u id;
    vec2d d;
    vec1u rid;
    vec1d rd;

    // Reflection data
    MEMBERS1(id, d, rid, rd);
    MEMBERS2("qxmatch_res", MAKE_MEMBER(id), MAKE_MEMBER(d), MAKE_MEMBER(rid), MAKE_MEMBER(rd));
};

void qxmatch_save(const std::string& file, const qxmatch_res& r) {
    fits::write_table(file, ftable(r.id, r.d, r.rid, r.rd));
}

qxmatch_res qxmatch_restore(const std::string& file) {
    qxmatch_res r;
    fits::read_table(file, ftable(r.id, r.d, r.rid, r.rd));
    return r;
}

struct qxmatch_params {
    uint_t thread = 1u;
    uint_t nth = 1u;
    bool verbose = false;
    bool self = false;
};

template<typename C1, typename C2,
    typename enable = typename std::enable_if<!is_vec<C1>::value>::type>
qxmatch_res qxmatch(const C1& cat1, const C2& cat2, qxmatch_params params = qxmatch_params{}) {
    return qxmatch(cat1.ra, cat1.dec, cat2.ra, cat2.dec, params);
}

template<typename C1, typename enable = typename std::enable_if<!is_vec<C1>::value>::type>
qxmatch_res qxmatch(const C1& cat1, qxmatch_params params = qxmatch_params{}) {
    return qxmatch(cat1.ra, cat1.dec, params);
}

template<typename TypeR1, typename TypeD1>
qxmatch_res qxmatch(const vec_t<1,TypeR1>& ra1, const vec_t<1,TypeD1>& dec1,
    qxmatch_params params = qxmatch_params{}) {
    params.self = true;
    return qxmatch(ra1, dec1, ra1, dec1, params);
}

template<typename TypeR1, typename TypeD1, typename TypeR2, typename TypeD2>
qxmatch_res qxmatch(const vec_t<1,TypeR1>& ra1, const vec_t<1,TypeD1>& dec1,
    const vec_t<1,TypeR2>& ra2, const vec_t<1,TypeD2>& dec2,
    qxmatch_params params = qxmatch_params{}) {

    phypp_check(ra1.dims == dec1.dims, "first RA and Dec dimensions do not match (",
        ra1.dims, " vs ", dec1.dims, ")");
    phypp_check(ra2.dims == dec2.dims, "second RA and Dec dimensions do not match (",
        ra2.dims, " vs ", dec2.dims, ")");

    const double d2r = 3.14159265359/180.0;
    auto dra1  = ra1*d2r;
    auto ddec1 = dec1*d2r;
    auto dcdec1 = cos(ddec1);
    auto dra2  = ra2*d2r;
    auto ddec2 = dec2*d2r;
    auto dcdec2 = cos(ddec2);

    const uint_t n1 = n_elements(ra1);
    const uint_t n2 = n_elements(ra2);

    const uint_t nth = clamp(params.nth, 1u, npos);

    qxmatch_res res;
    res.id = replicate(npos, nth, n1);
    res.d  = dblarr(nth, n1)+dinf;
    res.rid = replicate(npos, n2);
    res.rd  = dblarr(n2)+dinf;

    auto work = [&, nth] (uint_t i, uint_t j, qxmatch_res& tres) {
        // For each pair of source, compute a distance indicator.
        // Note that this is not the 'true' distance in arseconds, but this is sufficient
        // to find the nearest neighbors (the true distance is obtained by computing
        // 2*asin(sqrt(sd))), but all these functions are monotonous, hence not applying
        // them does not change the relative distances).
        double sra = sin(0.5*(dra2[j] - dra1[i]));
        double sde = sin(0.5*(ddec2[j] - ddec1[i]));
        double sd = sde*sde + sra*sra*dcdec2[j]*dcdec1[i];

        // We then compare this new distance to the largest one that is in the Nth nearest
        // neighbor list. If it is lower than that, we insert it in the list, removing the
        // old one, and sort the whole thing so that the largest distance goes as the end of
        // the list.
        if (sd < tres.d(nth-1,i)) {
            tres.id(nth-1,i) = j;
            tres.d(nth-1,i) = sd;
            uint_t k = nth-2;
            while (k != npos && tres.d(k,i) > tres.d(k+1,i)) {
                std::swap(tres.d(k,i), tres.d(k+1,i));
                std::swap(tres.id(k,i), tres.id(k+1,i));
                --k;
            }
        }
        if (sd < tres.rd[j]) {
            tres.rid[j] = i;
            tres.rd[j] = sd;
        }
    };

    if (params.thread <= 1) {
        // When using a single thread, all the work is done in the main thread
        auto p = progress_start(n1);
        for (uint_t i = 0; i < n1; ++i) {
            for (uint_t j = 0; j < n2; ++j) {
                if (params.self && i == j) continue;
                work(i,j,res);
            }

            if (params.verbose) print_progress(p, i);
        }
    } else {
        // When using more than one thread, the work load is evenly separated between all the
        // available threads, such that they should more or less all end at the same time.
        std::atomic<uint_t> iter(0);

        // Create the thread pool and launch the threads
        std::vector<qxmatch_res> vres(params.thread);
        for (auto& r : vres) {
            r.id = replicate(npos, nth, n1);
            r.d  = dblarr(nth, n1)+dinf;
            r.rid = replicate(npos, n2);
            r.rd  = dblarr(n2)+dinf;
        }
        vec1u tbeg(params.thread);
        vec1u tend(params.thread);
        auto pool = thread::pool(params.thread);
        uint_t total = 0;
        uint_t assigned = floor(n1/float(params.thread));
        for (uint_t t = 0; t < params.thread; ++t) {
            if (t == params.thread-1) {
                assigned = n1 - total;
            }

            tbeg[t] = total;
            tend[t] = total+assigned;

            pool[t].start([&iter, &work, &vres, t, total, params, assigned, n2]() {
                for (uint_t i = total; i < total+assigned; ++i) {
                    for (uint_t j = 0; j < n2; ++j) {
                        if (params.self && i == j) continue;
                        work(i,j,vres[t]);
                    }

                    ++iter;
                }
            });

            total += assigned;
        }

        // Wait for the computation to finish
        // Here the main thread does nothing except sleeping, occasionally waking up every second to
        // update the progress bar if any.
        auto p = progress_start(n1);
        while (iter < n1) {
            thread::sleep_for(0.2);
            if (params.verbose) print_progress(p, iter);
        }

        // By now, all threads should have ended their tasks.
        // We must ask them to terminate nicely.
        for (auto& t : pool) {
            t.join();
        }

        // Merge back the results of each thread
        for (uint_t t = 0; t < params.thread; ++t) {
            auto ids = rgen(tend[t] - tbeg[t]) + tbeg[t];
            res.id(_,ids) = vres[t].id(_,ids);
            res.d(_,ids) = vres[t].d(_,ids);

            if (t == 0) {
                res.rid = vres[t].rid;
                res.rd = vres[t].rd;
            } else {
                for (uint_t j = 0; j < n2; ++j) {
                    if (res.rd[j] >= vres[t].rd[j]) {
                        res.rid[j] = vres[t].rid[j];
                        res.rd[j] = vres[t].rd[j];
                    }
                }
            }
        }
    }

    // Convert the distance estimator to a real distance
    res.d = 3600.0*(180.0/3.14159265359)*2*asin(sqrt(res.d));
    res.rd = 3600.0*(180.0/3.14159265359)*2*asin(sqrt(res.rd));

    return res;
}

struct id_pair {
    vec1u id1, id2;
    vec1u lost;
};

id_pair xmatch_clean_best(const qxmatch_res& r) {
    const uint_t ngal = dim(r.id)[1];

    id_pair c;
    c.id1.data.reserve(ngal);
    c.id2.data.reserve(ngal);
    c.lost.data.reserve(ngal/6);

    for (uint_t i = 0; i < ngal; ++i) {
        if (r.rid[r.id[i]] == i) {
            c.id1.data.push_back(i);
            c.id2.data.push_back(r.id[i]);
        } else {
            c.lost.data.push_back(i);
        }
    }

    c.id1.dims[0] = c.id1.data.size();
    c.id2.dims[0] = c.id2.data.size();
    c.lost.dims[0] = c.lost.data.size();

    return c;
}

void xmatch_check_lost(const id_pair& p) {
    if (n_elements(p.lost) != 0) {
        warning(n_elements(p.lost), " sources failed to cross match");
    }
}

void xmatch_save_lost(const id_pair& p, const std::string& save) {
    if (n_elements(p.lost) != 0) {
        warning(n_elements(p.lost), " sources failed to cross match");
        fits::write_table(save, ftable(p.lost));
    }
}

template<typename TypeR1, typename TypeD1, typename TypeR2, typename TypeD2>
vec1d qxcor(const vec_t<1,TypeR1>& ra1, const vec_t<1,TypeD1>& dec1,
    const vec_t<1,TypeR2>& ra2, const vec_t<1,TypeD2>& dec2) {
    phypp_check(ra1.dims == dec1.dims, "first RA and Dec dimensions do not match (",
        ra1.dims, " vs ", dec1.dims, ")");
    phypp_check(ra2.dims == dec2.dims, "second RA and Dec dimensions do not match (",
        ra2.dims, " vs ", dec2.dims, ")");

    const uint_t n1 = n_elements(ra1);
    const uint_t n2 = n_elements(ra2);

    vec1d res;
    res.reserve(n1*n2);

    const double d2r = 3.14159265359/180.0;
    auto dra1  = ra1*d2r;
    auto ddec1 = dec1*d2r;
    auto dcdec1 = cos(ddec1);
    auto dra2  = ra2*d2r;
    auto ddec2 = dec2*d2r;
    auto dcdec2 = cos(ddec2);

    for (uint_t i1 : range(n1))
    for (uint_t i2 : range(n2)) {
        double sra = sin(0.5*(dra2[i2] - dra1[i1]));
        double sde = sin(0.5*(ddec2[i2] - ddec1[i1]));
        double sd = sde*sde + sra*sra*dcdec2[i2]*dcdec1[i1];
        res.push_back(3600.0*(180.0/3.14159265359)*2*asin(sqrt(sd)));
    }

    return res;
}

template<typename TypeR, typename TypeD>
vec1d qxcor_self(const vec_t<1,TypeR>& ra, const vec_t<1,TypeD>& dec) {
    phypp_check(ra.dims == dec.dims, "RA and Dec dimensions do not match (",
        ra.dims, " vs ", dec.dims, ")");

    const uint_t n = n_elements(ra);

    vec1d res;
    if (n < 2) return res;

    res.reserve(n*(n-1));

    const double d2r = 3.14159265359/180.0;
    auto dra  = ra*d2r;
    auto ddec = dec*d2r;
    auto dcdec = cos(ddec);

    for (uint_t i1 : range(n))
    for (uint_t i2 : range(i1+1, n)) {
        double sra = sin(0.5*(dra[i2] - dra[i1]));
        double sde = sin(0.5*(ddec[i2] - ddec[i1]));
        double sd = sde*sde + sra*sra*dcdec[i2]*dcdec[i1];
        res.push_back(3600.0*(180.0/3.14159265359)*2*asin(sqrt(sd)));
    }

    return res;
}

struct comment_pool_t {
    std::map<std::string, std::string> var_pool;
    std::vector<std::string> com_list;
};

comment_pool_t create_comment_pool() {
    return comment_pool_t();
}

void add_comment(comment_pool_t& pool, const std::string& com) {
    pool.com_list.push_back(com);
}

void add_comment(comment_pool_t& pool, const std::string& var, const std::string& com) {
    pool.var_pool.insert(std::make_pair(var, com));
}

std::string build_comments(const comment_pool_t& pool, uint_t width = 80) {
    std::string str;

    auto add_cmt = [&str, width] (const std::string& c, const std::string& tab = "") {
        vec1s spl = wrap(c, width, tab);
        for (uint_t i = 0; i < spl.size(); ++i) {
            str += spl[i] + '\n';
        }
    };

    for (auto& s : pool.com_list) {
        add_cmt(s);
    }

    vec1s tree;
    for (auto& p : pool.var_pool) {
        vec1s ltree = trim(split(p.first, "."));
        uint_t nmax = std::min(tree.size(), ltree.size());
        uint_t icom;
        std::string tab = "";
        for (icom = 0; icom < nmax; ++icom) {
            if (tree[icom] != ltree[icom]) break;
            if (icom % 2 == 0) {
                tab += "| ";
            } else {
                tab += ": ";
            }
        }

        uint_t i;
        for (i = icom; i < ltree.size()-1; ++i) {
            str += tab+ltree[i]+'\n';
            if (i % 2 == 0) {
                tab += "| ";
            } else {
                tab += ": ";
            }
        }

        std::string header = ltree[i]+": ";
        add_cmt(tab+header+p.second, tab+std::string(header.size(), ' '));

        tree = ltree;
    }

    return str;
}

struct catalog_pool;

struct catalog_t {
    catalog_pool& pool;

    vec1u       idi;
    vec1u       idm;
    vec2d       d;

    std::string name;
    vec1s       sources;
    vec1s       files;
    std::string comment;
    std::string ref;

    template<std::size_t Dim, typename T, typename U, typename V,
        typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& in, const vec_t<Dim,U>& out, const V& def);

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& in, const U& out, const V& def);

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& in, const vec_t<1,U>& out, const V& def);

    template<std::size_t Dim, typename T, typename U, typename V,
        typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& in, const vec_t<Dim,U>& out, const V& def, const std::string& com);

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& in, const U& out, const V& def, const std::string& com);

    template<typename T, typename U, typename V>
    void merge(T& in, const U& out, const V& def);

    template<typename T, typename U, typename V>
    void merge(T& in, const U& out, const V& def, const std::string& com);

    template<typename T, typename U>
    void assign(T& in, const U& out, const std::string& com);

    void merge_flux(const vec2f& flux, const vec2f& err, const vec1s& bands, const vec1s& notes);
};

struct catalog_pool {
    comment_pool_t coms;
    uint_t ngal;
    std::list<catalog_t> pool;

    vec1u origin;
    vec1d ra, dec;
    vec2f flux, flux_err;
    vec1s bands, notes;

    bool xmatch;
    std::string xmatch_file;

    catalog_pool(bool xmatch_, const std::string& xmatch_file_) :
        xmatch(xmatch_), xmatch_file(xmatch_file_) {
        file::mkdir(file::get_directory(xmatch_file));
    }

    catalog_t& add_catalog(const std::string& name, const vec1s& sources, const vec1s& files,
        const std::string& comment = "") {

        std::string ref = "["+strn(pool.size()+1)+"]";
        pool.push_back({*this, uindgen(ngal), uindgen(ngal), {}, name, sources, files, comment, ref});
        return pool.back();
    }

private :
    catalog_t& add_catalog_(const vec1d& cra, const vec1d& cdec, bool no_new, vec1u sel,
        const std::string& name, const vec1s& sources, const vec1s& files,
        const std::string& comment = "") {

        phypp_check(cra.size() == cdec.size(), "need ra.size() == dec.size()");

        if (sel.empty()) {
            sel = uindgen(cra.size());
        }

        if (pool.empty()) {
            ngal = sel.size();
            origin = replicate(1, ngal);
            ra = cra[sel];
            dec = cdec[sel];
            vec1u idm = uindgen(ngal);
            std::string ref = "[1]";
            pool.push_back({*this, sel, idm, {}, name, sources, files, comment, ref});
            return pool.back();
        } else {
            print("cross-matching "+name+"...");

            std::string file = xmatch_file+"_"+strn(std::hash<std::string>()(
                name+strn(cra.size())+strn(sel.size())+strn(ngal)+strn(sources)+strn(files)
            ))+".fits";

            bool rematch = false;
            if (!xmatch && file::exists(file)) {
                struct {
                    vec1u sel;
                    uint_t ngal;
                } tmp;

                fits::read_table(file, ftable(tmp.sel, tmp.ngal));
                if (tmp.sel.size() != sel.size() || tmp.ngal != ngal) {
                    warning("incompatible cross-match data, re-matching");
                    rematch = true;
                }
            } else {
                rematch = true;
            }

            qxmatch_res xm;

            if (rematch) {
                qxmatch_params p; p.nth = 2; p.thread = 4; p.verbose = true;
                xm = qxmatch(cra[sel], cdec[sel], ra, dec, p);

                fits::write_table(file, ftable(
                    sel, ngal, xm.id, xm.d, xm.rid, xm.rd
                ));
            } else {
                print("reading data from "+file);
                fits::read_table(file, xm);
            }

            const uint_t n = sel.size();
            vec1u idm; idm.reserve(n);
            vec1u idn;
            vec1u tsel;
            if (no_new) {
                tsel.reserve(n);
            } else {
                tsel = uindgen(sel.size());
            }

            for (uint_t i = 0; i < n; ++i) {
                if (xm.rid[xm.id(0,i)] == i) {
                    idm.push_back(xm.id(0,i));
                    if (no_new) {
                        tsel.push_back(i);
                    }
                } else {
                    if (!no_new) {
                        idn.push_back(sel[i]);
                        idm.push_back(ra.size());
                        ra.push_back(cra[sel[i]]);
                        dec.push_back(cdec[sel[i]]);
                        origin.push_back(pool.size()+1);
                    }
                }
            }

            vec2d td = replicate(dnan, xm.d.dims[0], cra.size());
            td(_,sel[tsel]) = xm.d(_,tsel);

            sel = sel[tsel];

            if (no_new) {
                print("> ", idm.size(), " matched sources, ", n - idm.size(), " lost");
            } else {
                print("> ", n - idn.size(), " matched sources, ", idn.size(), " new sources");
                ngal += idn.size();
            }

            std::string ref = "["+strn(pool.size()+1)+"]";
            pool.push_back({*this, sel, idm, std::move(td), name, sources, files, comment, ref});
            return pool.back();
        }
    }

public :
    catalog_t& add_catalog(const vec1d& cra, const vec1d& cdec, vec1u sel,
        const std::string& name, const vec1s& sources, const vec1s& files,
        const std::string& comment = "") {
        return add_catalog_(cra, cdec, false, sel, name, sources, files, comment);
    }

    catalog_t& match_catalog(const vec1d& cra, const vec1d& cdec, vec1u sel,
        const std::string& name, const vec1s& sources, const vec1s& files,
        const std::string& comment = "") {
        return add_catalog_(cra, cdec, true, sel, name, sources, files, comment);
    }

    std::string build_catalog_list(uint_t width = 80) {
        std::string str;

        for (auto& c : pool) {
            str += c.ref+": "+c.name+"\n";
            if (!c.sources.empty()) {
                if (c.sources.size() == 1) {
                    str += "  source: "+c.sources[0]+"\n";
                } else {
                    std::string header = "  sources: ";
                    str += header+c.sources[0]+"\n";
                    for (uint_t i = 1; i < c.sources.size(); ++i) {
                        str += std::string(header.size(), ' ')+c.sources[i]+"\n";
                    }
                }
            }
            if (!c.files.empty()) {
                if (c.files.size() == 1) {
                    str += "  file: "+c.files[0]+"\n";
                } else {
                    std::string header = "  files: ";
                    str += header+c.files[0]+"\n";
                    for (uint_t i = 1; i < c.files.size(); ++i) {
                        str += std::string(header.size(), ' ')+c.files[i]+"\n";
                    }
                }
            }
            if (!c.comment.empty()) {
                vec1s spl = wrap("  comment: "+c.comment, width, std::string(4, ' '));
                for (uint_t i = 0; i < spl.size(); ++i) {
                    str += spl[i] + '\n';
                }
            }
            str += "\n";
        }

        return str;
    }
};

template<std::size_t Dim, typename T, typename U, typename V, typename enable>
void catalog_t::merge(vec_t<Dim,T>& in, const vec_t<Dim,U>& out, const V& def) {
    phypp_check(in.empty(), name+": merging twice into the same variable");
    in.dims = out.dims;
    in.dims[0] = pool.ngal;
    in.resize();
    in[_] = def;
    in(idm, repeat<Dim-1>(_)) = out(idi, repeat<Dim-1>(_));
}

template<typename T, typename U, typename V>
void catalog_t::merge(vec_t<1,T>& in, const U& out, const V& def) {
    phypp_check(in.empty(), name+": merging twice into the same variable");
    in = replicate(def, pool.ngal);
    in[idm] = out;
}

template<typename T, typename U, typename V>
void catalog_t::merge(vec_t<1,T>& in, const vec_t<1,U>& out, const V& def) {
    phypp_check(in.empty(), name+": merging twice into the same variable");
    in = replicate(def, pool.ngal);
    in[idm] = out[idi];
}

template<std::size_t Dim, typename T, typename U, typename V, typename enable>
void catalog_t::merge(vec_t<Dim,T>& in, const vec_t<Dim,U>& out, const V& def,
    const std::string& com) {
    merge(in, out, def);
    add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
}

template<typename T, typename U, typename V>
void catalog_t::merge(vec_t<1,T>& in, const U& out, const V& def, const std::string& com) {
    merge(in, out, def);
    add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
}

namespace astro_impl {
    template<typename M, typename V>
    struct do_catalog_merge_run {
        const reflex::member_t& m;
        const M& v;
        const V& def;
        catalog_t& cat;

        void operator () (reflex::member_t& n, M& p) {
            if (this->m.name == n.name) {
                this->cat.merge(p, this->v, this->def);
            }
        }

        template<typename P, typename enable2 =
            typename std::enable_if<reflex::is_struct<M>::value>::type>
        void operator () (reflex::member_t& n, reflex::struct_t<P> p) {
            this->cat.merge(p, this->v, this->def);
        }

        template<typename P>
        void operator () (reflex::member_t& n, P&& p) {
            phypp_check(this->m.name != n.name, "incompatible types in merging '",
                this->m.full_name(), "' into '", n.full_name(), "'"
            );
        }
    };

    template<typename T, typename V>
    struct do_catalog_merge {
        reflex::struct_t<T> t;
        const V& def;
        catalog_t& cat;

        template<typename M>
        void operator () (const reflex::member_t& m, const M& v) {
            do_catalog_merge_run<M,V> run{m, v, this->def, this->cat};
            reflex::foreach_member(this->t, run);
        }
    };
}

template<typename T, typename U, typename V>
void catalog_t::merge(T& in, const U& out, const V& def) {
    astro_impl::do_catalog_merge<T,V> run{reflex::wrap(in), def, *this};
    reflex::foreach_member(reflex::wrap(out), run);
}

template<typename T, typename U, typename V>
void catalog_t::merge(T& in, const U& out, const V& def, const std::string& com) {
    merge(in, out, def);
    add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
}

template<typename T, typename U>
void catalog_t::assign(T& in, const U& out, const std::string& com) {
    in = out;
    add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
}

void catalog_t::merge_flux(const vec2f& flux, const vec2f& err, const vec1s& bands,
    const vec1s& notes) {

    phypp_check(flux.dims[1] == err.dims[1], name+": flux and error do not match");
    phypp_check(flux.dims[1] == bands.size(), name+": flux and band do not match");

    vec1s tnotes = notes;
    if (tnotes.empty()) {
        tnotes = strarr(bands.size());
    }

    vec1u idne = where(!empty(tnotes));
    tnotes[idne] = tnotes[idne]+" ";
    tnotes += ref;

    vec2f tflux, terr;
    merge(tflux, flux, fnan);
    merge(terr, err, fnan);

    append(pool.bands, bands);
    append(pool.notes, tnotes);
    append<1>(pool.flux, tflux);
    append<1>(pool.flux_err, terr);
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

        auto cut = img(x[i]+r, y[i]+r).concretise();

        // Discard any source that contains a bad pixel (either infinite or NaN)
        if (total(!finite(cut)) != 0) {
            continue;
        }

        ids.push_back(i);
        cube.push_back(cut);
    }
}

struct qstack_params {
    bool keepnan = false;
};

template<typename Type>
void qstack(const vec1d& ra, const vec1d& dec, const std::string& filename, uint_t hsize,
    vec_t<3,Type>& cube, vec1u& ids, qstack_params params = qstack_params()) {

    phypp_check(file::exists(filename), "cannot stack on inexistant file '"+filename+"'");
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    // Open the FITS file
    fitsfile* fptr;
    int status = 0;

    fits_open_image(&fptr, filename.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

    // Read the header as a string and read the WCS data
    char* hstr = nullptr;
    int nkeys  = 0;
    fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
    fits::wcs astro(hstr);
    free(hstr);

    // Get the dimensions of the image
    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    phypp_check(naxis == 2, "cannot stack on image cubes (image dimensions: "+strn(naxis)+")");
    long naxes[2];
    fits_get_img_size(fptr, naxis, naxes, &status);
    long width = naxes[0], height = naxes[1];

    // Convert ra/dec to x/y
    vec1d x, y;
    fits::ad2xy(astro, ra, dec, x, y);

    // Allocate memory to hold all the cutouts
    if (cube.empty()) {
    	cube.dims[1] = cube.dims[2] = 2*hsize+1;
    }

    cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    ids.reserve(ids.size() + ra.size());

    // Loop over all sources
    for (uint_t i = 0; i < ra.size(); ++i) {
        long p0[2] = {long(round(x[i]-hsize)), long(round(y[i]-hsize))};
        long p1[2] = {long(round(x[i]+hsize)), long(round(y[i]+hsize))};

        // Discard any source that falls out of the boundaries of the image
        if (p0[0] < 1 || p1[0] >= width || p0[1] < 1 || p1[1] >= height) {
            continue;
        }

        vec_t<2,Type> cut(2*hsize+1, 2*hsize+1);

        Type null = fnan;
        int anynul = 0;
        long inc[2] = {1, 1};
        fits_read_subset(fptr, fits::traits<Type>::ttype, p0, p1, inc, &null,
            cut.data.data(), &anynul, &status);

        // Discard any source that contains a bad pixel (either infinite or NaN)
        if (!params.keepnan && total(!finite(cut)) != 0) {
            continue;
        }

        ids.push_back(i);
        cube.push_back(cut);
    }

    fits_close_file(fptr, &status);
}

template<typename Type>
void qstack(const vec1d& ra, const vec1d& dec, const std::string& ffile, const std::string& wfile,
	uint_t hsize, vec_t<3,Type>& cube, vec_t<3,Type>& wcube, vec1u& ids,
    qstack_params params = qstack_params()) {

    phypp_check(file::exists(ffile), "cannot stack on inexistant file '"+ffile+"'");
    phypp_check(file::exists(wfile), "cannot stack on inexistant file '"+wfile+"'");
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    // Open the FITS file
    fitsfile* fptr;
    fitsfile* wfptr;
    int status = 0;

    fits_open_image(&fptr, ffile.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+ffile+"'");
    fits_open_image(&wfptr, wfile.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+wfile+"'");

    // Read the header as a string and read the WCS data
    char* hstr = nullptr;
    int nkeys  = 0;
    fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
    fits::wcs astro(hstr);
    free(hstr);

    // Get the dimensions of the image
    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    phypp_check(naxis == 2, "cannot stack on image cubes (image dimensions: "+strn(naxis)+")");
    long naxes[2];
    fits_get_img_size(fptr, naxis, naxes, &status);
    long width = naxes[0], height = naxes[1];
    long wnaxes[2];
    fits_get_img_size(wfptr, 2, wnaxes, &status);
    phypp_check(naxes[0] == wnaxes[0] && naxes[1] == wnaxes[1], "image and weight map do not match");

    // Convert ra/dec to x/y
    vec1d x, y;
    fits::ad2xy(astro, ra, dec, x, y);

    // Allocate memory to hold all the cutouts
    if (cube.empty()) {
    	cube.dims[1] = cube.dims[2] = 2*hsize+1;
    }
    if (wcube.empty()) {
    	wcube.dims[1] = wcube.dims[2] = 2*hsize+1;
    }

    cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    wcube.reserve(wcube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    ids.reserve(ids.size() + ra.size());

    // Loop over all sources
    for (uint_t i = 0; i < ra.size(); ++i) {
        long p0[2] = {long(round(x[i]-hsize)), long(round(y[i]-hsize))};
        long p1[2] = {long(round(x[i]+hsize)), long(round(y[i]+hsize))};

        // Discard any source that falls out of the boundaries of the image
        if (p0[0] < 1 || p1[0] >= width || p0[1] < 1 || p1[1] >= height) {
            continue;
        }

        vec_t<2,Type> cut(2*hsize+1, 2*hsize+1);
        vec_t<2,Type> wcut(2*hsize+1, 2*hsize+1);

        Type null = fnan;
        int anynul = 0;
        long inc[2] = {1, 1};

        fits_read_subset(fptr,  fits::traits<Type>::ttype, p0, p1, inc, &null,
            cut.data.data(),  &anynul, &status);
        fits_read_subset(wfptr, fits::traits<Type>::ttype, p0, p1, inc, &null,
            wcut.data.data(), &anynul, &status);

        // Discard any source that contains a bad pixel (either infinite or NaN)
        if (!params.keepnan && total(!finite(cut) || !finite(wcut)) != 0) {
            continue;
        }

        ids.push_back(i);
        cube.push_back(cut);
        wcube.push_back(wcut);
    }

    fits_close_file(fptr, &status);
    fits_close_file(wfptr, &status);
}

template<typename Type>
vec_t<2,rtype_t<Type>> qstack_mean(const vec_t<3,Type>& fcube) {
    return mean(fcube, 0);
}

template<typename TypeF, typename TypeW>
vec_t<2,rtype_t<TypeF>> qstack_mean(const vec_t<3,TypeF>& fcube, const vec_t<3,TypeW>& wcube) {
    return total(fcube*wcube, 0)/total(wcube, 0);
}

template<typename Type>
vec_t<2,rtype_t<Type>> qstack_median(const vec_t<3,Type>& fcube) {
    return median(fcube, 0);
}

template<typename Type, typename TypeS, typename F>
void qstack_bootstrap(const vec_t<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed, F&& func) {

    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        auto tfcube = fcube(ids,_,_).concretise();
        func(tfcube);
    }
}

template<typename TypeF, typename TypeW, typename TypeS, typename F>
void qstack_bootstrap(const vec_t<3,TypeF>& fcube, const vec_t<3,TypeW>& wcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed, F&& func) {

    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        auto tfcube = fcube(ids,_,_).concretise();
        auto twcube = wcube(ids,_,_).concretise();
        func(tfcube, twcube);
    }
}

template<typename Type>
vec_t<3,rtype_t<Type>> qstack_bootstrap_apply_id_(const vec1u& ids, const vec_t<3,Type>& cube) {
    return cube(ids,_,_).concretise();
}

template<typename Type, typename ... Args>
uint_t qstack_bootstrap_get_size_(const vec_t<3,Type>& cube, const Args& ... cubes) {
    return cube.dims[0];
}

template<typename TypeS, typename F, typename ... Args>
void qstack_bootstrap(uint_t nbstrap, uint_t nsel, TypeS& seed, F&& func, const Args& ... cubes) {
    const uint_t nsrc = qstack_bootstrap_get_size_(cubes...);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, nsrc-1, nsel);
        func(qstack_bootstrap_apply_id_(ids, cubes)...);
    }
}

template<typename Type, typename TypeS>
vec_t<3,rtype_t<Type>> qstack_mean_bootstrap(const vec_t<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed) {

    vec_t<3,rtype_t<Type>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_mean(fcube(ids,_,_)));
    }

    return bs;
}

template<typename TypeF, typename TypeW, typename TypeS>
vec_t<3,rtype_t<TypeF>> qstack_mean_bootstrap(const vec_t<3,TypeF>& fcube,
    const vec_t<3,TypeW>& wcube, uint_t nbstrap, uint_t nsel, TypeS& seed) {

    vec_t<3,rtype_t<TypeF>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_mean(fcube(ids,_,_), wcube(ids,_,_)));
    }

    return bs;
}

template<typename Type, typename TypeS>
vec_t<3,rtype_t<Type>> qstack_median_bootstrap(const vec_t<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed) {

    vec_t<3,rtype_t<Type>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_median(fcube(ids,_,_)));
    }

    return bs;
}

// Compute the area covered by a field given a set of source coordinates [deg^2] and pre-computed
// convex hull (as obtained from convex_hull() with the same coordinates).
// Coordinates are assumed to be given in degrees.
// NOTE: this algorithm assumes that the field is well characterized by a *convex* hull, meaning
// that is has no hole (masked stars, ...) and the borders are not concave (no "zigzag" shape, ...).
// If these hypotheses do not hold, use field_area_h2d instead.
template<typename TX, typename TY, typename TH>
auto field_area_hull(const TX& ra, const TY& dec, const TH& hull) {
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
auto field_area_hull(const TX& ra, const TY& dec) {
    return field_area_hull(ra, dec, convex_hull(ra, dec));
}

// Compute the area covered by a field given a set of source coordinates [deg^2] by iteratively
// building a 2d histogram of point sources.
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
auto field_area_h2d(const TX& ra, const TY& dec) {
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    // Get rough extents of the field
    double min_ra  = min(ra),  max_ra  = max(ra);
    double min_dec = min(dec), max_dec = max(dec);
    double rough_area = (angdist(min_ra, min_dec, max_ra, min_dec)/3600.0)*
                        (angdist(min_ra, min_dec, min_ra, max_dec)/3600.0);

    // Optimal bin size
    // If field is perfectly uniform, then all cells of the field will contain 'optinum' objects
    uint_t optinum = 10;
    double dr = sqrt(optinum)*sqrt(rough_area/ra.size());

    // Since the field is possibly not uniform, we have to introduce a tolerence threshold.
    uint_t threshold = 0.5*optinum;

    // Now define the grid
    uint_t nra  = (max_ra  - min_ra)/dr;
    uint_t ndec = (max_dec - min_dec)/dr;
    vec2d rb = make_bins(min_ra,  max_ra,  nra);
    vec2d db = make_bins(min_dec, max_dec, ndec);

    // Compute the area occupied by a cell
    double cell_area = (rb(0,1) - rb(0,0))*(db(0,1) - db(0,0));

    // Build the 2d histogram
    vec2u counts = histogram2d(ra, dec, rb, db);

    // Sum up the areas of all cells with counts above the threshold
    return cell_area*total(counts > threshold);
}

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
auto field_area(const TX& ra, const TY& dec) {
    return field_area_h2d(ra, dec);
}

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename T>
auto field_area(const T& t) {
    return field_area(t.ra, t.dec);
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
auto sed2flux(const filter_t& filter, const vec_t<1,TypeL>& lam, const vec_t<1,TypeS>& sed) {
    if (!(filter.lam[0] > lam[0] && filter.lam[-1] < lam[-1])) {
        return dnan;
    }

    uint_t sp = lower_bound(filter.lam[0], lam);
    uint_t ep = upper_bound(filter.lam[-1], lam);

    vec1d nlam; nlam.reserve(filter.lam.size() + ep - sp + 1);
    vec1d nrs; nrs.reserve(filter.lam.size() + ep - sp + 1);

    uint_t j = sp;
    for (uint_t i : range(filter.lam)) {
        nlam.push_back(filter.lam[i]);
        nrs.push_back(filter.res[i]*interpolate(
            sed[j-1], sed[j], lam[j-1], lam[j], filter.lam[i]
        ));

        if (i != filter.lam.size() - 1) {
            while (lam[j] < filter.lam[i+1]) {
                nlam.push_back(lam[j]);
                nrs.push_back(sed[j]*interpolate(
                    filter.res[i], filter.res[i+1], filter.lam[i], filter.lam[i+1], lam[j]
                ));
                ++j;
            }
        }
    }

    return integrate(nlam, nrs);
}

template<typename TypeL, typename TypeS>
auto sed2flux(const filter_t& filter, const vec_t<2,TypeL>& lam, const vec_t<2,TypeS>& sed) {
    using rtype = decltype(sed[0]*filter.res[0]);
    const uint_t nsed = sed.dims[0];
    vec_t<1,rtype> r; r.reserve(nsed);

    for (uint_t s = 0; s < nsed; ++s) {
        r.push_back(sed2flux(filter, lam(s,_).concretise(), sed(s,_).concretise()));
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

template<typename TypeLib, typename TypeFi>
vec2d template_observed(const TypeLib& lib, double z, double d,
    const vec_t<1,TypeFi>& filters) {

    // Move each SED to the observed frame
    vec2d rflam = lib.lam*(1.0 + z);
    vec2d rfsed = lsun2uJy(z, d, lib.lam, lib.sed);

    // Convolve each SED with the response curve of the filters
    const uint_t nsed = rfsed.dims[0];
    const uint_t nfilter = filters.size();
    vec2d flux = dblarr(nsed, nfilter);
    for (uint_t f = 0; f < nfilter; ++f) {
        flux(_,f) = sed2flux(filters[f], rflam, rfsed);
    }

    return flux;
}

template<typename TypeLib, typename TypeFi, typename TypeZ, typename TypeD>
vec2d template_observed(const TypeLib& lib, const vec_t<1,TypeZ>& z, const vec_t<1,TypeD>& d,
    const vec_t<1,TypeFi>& filters) {

    phypp_check(z.size() == d.size(),
        "incompatible redshift and distance variables ("+strn(z.dims)+" vs "+strn(d.dims)+")");

    vec2d rflam = lib.lam*(1.0 + mean(z));
    vec2d rfsed = dblarr(lib.sed.dims);

    // Combine the SEDs of each source to build an effective "redshift convolved" SED in the
    // observed frame
    for (uint_t i = 0; i < z.size(); ++i) {
        auto tlam = lib.lam*(1.0 + z[i]);
        auto tsed = lsun2uJy(z[i], d[i], lib.lam, lib.sed);
        rfsed += interpolate(tsed, tlam, rflam);
    }

    rfsed /= z.size();

    // Convolve each SED with the response curve of the filters
    const uint_t nsed = rfsed.dims[0];
    const uint_t nfilter = filters.size();
    vec2d flux = dblarr(nsed, nfilter);
    for (uint_t f = 0; f < nfilter; ++f) {
        flux(_,f) = sed2flux(filters[f], rflam, rfsed);
    }

    return flux;
}

template<typename TypeLib, typename TypeFi, typename TypeZ, typename TypeD>
vec2d template_observed(const TypeLib& lib, const vec_t<2,TypeZ>& z, const vec_t<1,TypeD>& d,
    const vec_t<1,TypeFi>& filters) {

    phypp_check(z.dims[1] == d.dims[0],
        "incompatible redshift and distance variables ("+strn(z.dims)+" vs "+strn(d.dims)+")");

    vec2d rflam = lib.lam*(1.0 + mean(z));
    vec2d rfsed = dblarr(lib.sed.dims);

    // Combine the SEDs of each source to build an effective "redshift convolved" SED in the
    // observed frame
    double weight = 0.0;
    for (uint_t i = 0; i < z.dims[1]; ++i) {
        auto tlam = lib.lam*(1.0 + z[i]);
        auto tsed = z(1,i)*lsun2uJy(z(0,i), d[i], lib.lam, lib.sed);
        for (uint_t j = 0; j < lib.sed.dims[0]; ++j) {
            rfsed(j,_) += interpolate(tsed(j,_), tlam(j,_), rflam(j,_));
        }

        weight += z(1,i);
    }

    rfsed /= weight;

    // Convolve each SED with the response curve of the filters
    const uint_t nsed = rfsed.dims[0];
    const uint_t nfilter = filters.size();
    vec2d flux = dblarr(nsed, nfilter);
    for (uint_t f = 0; f < nfilter; ++f) {
        flux(_,f) = sed2flux(filters[f], rflam, rfsed);
    }

    return flux;
}

struct template_fit_res_t {
    uint_t bfit; // index of the best fit template in the library
    vec1d chi2;  // chi^2 of each template
    vec1d amp;   // renormalization amplitude of each template
    vec2d flux;  // best-fit flux for each template

    vec1u sed_sim;       // index of each error realization's best fit template
    vec1d amp_bfit_sim;  // renormalization amplitude of the best fit for each error realization
    vec1d amp_sim;       // renormalization amplitude for each error realization's best fit
};

struct template_fit_params {
    uint_t nsim = 200;   // number of random realizations to perform to estimate errors
    bool renorm = false; // if true, allow templates to be renormalized when fitted
};

template<typename TypeLib, typename TypeF, typename TypeE, typename TypeFi, typename TypeSeed,
    typename TypeZ, typename TypeD>
template_fit_res_t template_fit(const TypeLib& lib, TypeSeed& seed, const TypeZ& z,
   const TypeD& d, const vec_t<1,TypeF>& flux, const vec_t<1,TypeE>& err,
   const vec_t<1,TypeFi>& filters, template_fit_params params = template_fit_params()) {

    template_fit_res_t res;
    res.flux = template_observed(lib, z, d, filters);

    // Compute chi2 & renormalization factor
    auto weight = invsqr(err);
    using ttype = decltype(weight[0]*flux[0]*res.flux[0]);

    const uint_t nsed = res.flux.dims[0];
    const uint_t nfilter = filters.size();
    vec_t<1,ttype> tmp1(nsed), tmp2(nsed);
    for (uint_t i = 0; i < nsed; ++i) {
        auto tmp = res.flux(i,_);
        tmp1[i] = total(weight*flux*tmp);
        tmp2[i] = total(weight*tmp*tmp);
    }

    res.amp = tmp1/tmp2;
    if (params.renorm) {
        tmp1 *= res.amp;
        res.chi2 = total(weight*flux*flux) - tmp1;
    } else {
        res.chi2 = total(weight*flux*flux) - 2.0*tmp1 + tmp2;
    }

    // Find the best chi2 among all the SEDs
    res.bfit = min_id(res.chi2);

    // Compute the errors on the fit by adding a random offset to the measured photometry according
    // to the provided error. The fit is performed on each of these random realizations and the
    // error on the parameters are computed as the standard deviation of the fit results over all
    // the realizations.
    res.amp_bfit_sim.resize(params.nsim);
    res.amp_sim.resize(params.nsim);
    res.sed_sim.resize(params.nsim);

    const uint_t nflux = flux.size();
    for (uint_t i = 0; i < params.nsim; ++i) {
        auto fsim = flux + randomn(seed, nflux)*err;
        for (uint_t t = 0; t < nsed; ++t) {
            tmp1[t] = total(weight*fsim*res.flux(t,_));
        }

        auto amp = tmp1/tmp2;
        vec_t<1,ttype> chi2;
        if (params.renorm) {
            tmp1 *= amp;
            chi2 = total(weight*fsim*fsim) - tmp1;
        } else {
            chi2 = total(weight*fsim*fsim) - 2.0*tmp1 + tmp2;
        }

        auto ised = min_id(chi2);

        res.amp_bfit_sim[i] = amp[res.bfit];
        res.amp_sim[i] = amp[ised];
        res.sed_sim[i] = ised;
    }

    for (uint_t f = 0; f < nfilter; ++f) {
        res.flux(_,f) *= res.amp;
    }

    return res;
}

// Compute the LIR luminosity of a rest-frame SED (8um to 1000um).
template<typename TL, typename TS>
double lir_8_1000(const vec_t<1,TL>& lam, const vec_t<1,TS>& sed) {
    uint_t s = lower_bound(8.0, lam);
    uint_t e = upper_bound(1000.0, lam);
    if (s == npos || e == npos) return dnan;

    vec1u id = rgen(s, e);
    return integrate(lam[id], sed[id]/lam[id]);
}

bool make_psf(vec1u dims, double x0, double y0, const std::string& psf_model, vec2d& psf) {
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
