#ifndef ASTRO_HPP
#define ASTRO_HPP

#include <vec.hpp>
#include <math.hpp>
#include <print.hpp>
#include <image.hpp>
#include <thread.hpp>
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

// Convert a set of sexagesimal coordinates ('hh:mm:ss.ms') into degrees
void sex2deg(const std::string& sra, const std::string& sdec, double& ra, double& dec) {
    vec1s starr1 = split(sra, ":");
    vec1s starr2 = split(sdec, ":");

    if (starr1.size() != 3 || starr2.size() != 3) {
        ra = dnan;
        dec = dnan;
        return;
    }

    int_t rah, ram, dech, decm;
    double ras, decs;

    from_string(starr1[0], rah);
    from_string(starr1[1], ram);
    from_string(starr1[2], ras);
    from_string(starr2[0], dech);
    from_string(starr2[1], decm);
    from_string(starr2[2], decs);

    double signr = sign(rah);
    double signd = sign(dech);

    ra = (rah + ram*signr/60.0 + ras*signr/3600.0)*15.0;
    dec = dech + decm*signd/60.0 + decs*signd/3600.0;
}

template<std::size_t Dim, typename TSR, typename TSD, typename TR, typename TD>
void sex2deg(const vec_t<Dim,TSR>& sra, const vec_t<Dim,TSD>& sdec, vec_t<Dim,TR>& ra, vec_t<Dim,TD>& dec) {
    phypp_check(sra.size() == sdec.size(), "sex2deg: RA and Dec dimensions do not match (",
        sra.dims, " vs ", sdec.dims, ")");

    ra.resize(sra.dims);
    dec.resize(sra.dims);
    for (uint_t i : range(sra)) {
        sex2deg(sra[i], sdec[i], ra[i], dec[i]);
    }
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
    phypp_check(ra.size() == dec.size(), "deg2sex: RA and Dec dimensions do not match (",
        ra.dims, " vs ", dec.dims, ")");

    sra.resize(ra.dims);
    sdec.resize(ra.dims);
    for (uint_t i : range(ra)) {
        deg2sex(ra[i], dec[i], sra[i], sdec[i]);
    }
}

struct qxmatch_res {
    vec2u id;
    vec2d d;
    vec2u rid;
    vec2d rd;

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

template<typename C1, typename C2, typename ... Args,
    typename enable = typename std::enable_if<!is_vec<C1>::value>::type>
qxmatch_res qxmatch(const C1& cat1, const C2& cat2, Args&& ... args) {
    return qxmatch(cat1.ra, cat1.dec, cat2.ra, cat2.dec, std::forward<Args>(args)...);
}

template<typename TypeR1, typename TypeD1, typename TypeR2, typename TypeD2>
qxmatch_res qxmatch(const vec_t<1,TypeR1>& ra1, const vec_t<1,TypeD1>& dec1,
    const vec_t<1,TypeR2>& ra2, const vec_t<1,TypeD2>& dec2,
    declare_keywords(_thread(1u), _nth(1u), _verbose(false), _self(false))) {

    phypp_check(n_elements(ra1) == n_elements(dec1), "qxmatch: not as many RA as there are Dec");
    phypp_check(n_elements(ra2) == n_elements(dec2), "qxmatch: not as many RA as there are Dec");

    const double d2r = 3.14159265359/180.0;
    auto dra1  = ra1*d2r;
    auto ddec1 = dec1*d2r;
    auto dcdec1 = cos(ddec1);
    auto dra2  = ra2*d2r;
    auto ddec2 = dec2*d2r;
    auto dcdec2 = cos(ddec2);

    const uint_t n1 = n_elements(ra1);
    const uint_t n2 = n_elements(ra2);

    const uint_t nth = clamp(get_keyword(_nth), 1u, npos);
    const bool self = get_keyword(_self);

    qxmatch_res res;
    res.id = replicate(npos, nth, n1);
    res.d  = dblarr(nth, n1)+dinf;
    res.rid = replicate(npos, nth, n2);
    res.rd  = dblarr(nth, n2)+dinf;

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
        if (sd < tres.rd(nth-1,j)) {
            tres.rid(nth-1,j) = i;
            tres.rd(nth-1,j) = sd;
            uint_t k = nth-2;
            while (k != npos && tres.rd(k,j) > tres.rd(k+1,j)) {
                std::swap(tres.rd(k,j), tres.rd(k+1,j));
                std::swap(tres.rid(k,j), tres.rid(k+1,j));
                --k;
            }
        }
    };

    const uint_t nthread = get_keyword(_thread);
    if (nthread <= 1) {
        // When using a single thread, all the work is done in the main thread
        auto p = progress_start(n1);
        for (uint_t i = 0; i < n1; ++i) {
            for (uint_t j = 0; j < n2; ++j) {
                if (self && i == j) continue;
                work(i,j,res);
            }

            if (get_keyword(_verbose)) print_progress(p, i);
        }
    } else {
        // When using more than one thread, the work load is evenly separated between all the
        // available threads, such that they should more or less all end at the same time.
        std::atomic<uint_t> iter(0);

        // Create the thread pool and launch the threads
        std::vector<qxmatch_res> vres(nthread);
        for (auto& r : vres) {
            r.id = replicate(npos, nth, n1);
            r.d  = dblarr(nth, n1)+dinf;
            r.rid = replicate(npos, nth, n2);
            r.rd  = dblarr(nth, n2)+dinf;
        }
        vec1u tbeg(nthread);
        vec1u tend(nthread);
        auto pool = thread::pool(nthread);
        uint_t total = 0;
        uint_t assigned = floor(n1/float(nthread));
        for (uint_t t = 0; t < nthread; ++t) {
            if (t == nthread-1) {
                assigned = n1 - total;
            }

            tbeg[t] = total;
            tend[t] = total+assigned;

            pool[t].start([&iter, &work, &vres, t, total, self, assigned, n2]() {
                for (uint_t i = total; i < total+assigned; ++i) {
                    for (uint_t j = 0; j < n2; ++j) {
                        if (self && i == j) continue;
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
            if (get_keyword(_verbose)) print_progress(p, iter);
        }

        // By now, all threads should have ended their tasks.
        // We must ask them to terminate nicely.
        for (auto& t : pool) {
            t.join();
        }

        // Merge back the results of each thread
        for (uint_t t = 0; t < nthread; ++t) {
            auto ids = rgen(tend[t] - tbeg[t]) + tbeg[t];
            res.id(_,ids) = vres[t].id(_,ids);
            res.d(_,ids) = vres[t].d(_,ids);

            if (t == 0) {
                res.rid = vres[t].rid;
                res.rd = vres[t].rd;
            } else {
                for (uint_t j = 0; j < n2; ++j) {
                    for (uint_t n = 0; n < nth; ++n) {
                        if (res.rd(nth-1,j) < vres[t].rd(n,j)) break;

                        res.rid(nth-1,j) = vres[t].rid(n,j);
                        res.rd(nth-1,j) = vres[t].rd(n,j);
                        uint_t k = nth-2;
                        while (k != npos && res.rd(k,j) > res.rd(k+1,j)) {
                            std::swap(res.rd(k,j), res.rd(k+1,j));
                            std::swap(res.rid(k,j), res.rid(k+1,j));
                            --k;
                        }
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

void xmatch_check_lost(const id_pair& p, declare_keywords(_save(""))) {
    if (n_elements(p.lost) != 0) {
        warning(n_elements(p.lost), " sources failed to cross match");
        if (!get_keyword(_save).empty()) {
            fits::write_table(get_keyword(_save), ftable(p.lost));
        }
    }
}

template<typename TypeR1, typename TypeD1, typename TypeR2, typename TypeD2>
vec1d qxcor(const vec_t<1,TypeR1>& ra1, const vec_t<1,TypeD1>& dec1,
    const vec_t<1,TypeR2>& ra2, const vec_t<1,TypeD2>& dec2) {

    phypp_check(n_elements(ra1) == n_elements(dec1), "qxcor: not as many RA as there are Dec");
    phypp_check(n_elements(ra2) == n_elements(dec2), "qxcor: not as many RA as there are Dec");

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
    phypp_check(n_elements(ra) == n_elements(dec), "qxcor_self: not as many RA as there are Dec");

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

    template<std::size_t Dim, typename T, typename U, typename V,
        typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& in, const vec_t<Dim,U>& out, const V& def, const std::string& com);

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& in, const U& out, const V& def, const std::string& com);

    template<typename T, typename U, typename V>
    void merge(T& in, const U& out, const V& def);

    template<typename T, typename U, typename V>
    void merge(T& in, const U& out, const V& def, const std::string& com);

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
        pool.push_back({*this, uindgen(ngal), {}, name, sources, files, comment, ref});
        return pool.back();
    }

    catalog_t& add_catalog(const vec1d& cra, const vec1d& cdec, const std::string& name,
        const vec1s& sources, const vec1s& files, const std::string& comment = "") {

        phypp_check(cra.size() == cdec.size(), "add_catalog: need ra.size() == dec.size()");

        if (pool.empty()) {
            ngal = cra.size();
            origin = replicate(1, ngal);
            ra = cra;
            dec = cdec;
            vec1u idm = uindgen(ngal);
            std::string ref = "[1]";
            pool.push_back({*this, idm, {}, name, sources, files, comment, ref});
            return pool.back();
        } else {
            print("cross-matching "+name+"...");

            std::string file = xmatch_file+"_"+strn(std::hash<std::string>()(
                name+strn(sources)+strn(files)+comment
            ))+".fits";

            bool rematch = false;
            if (!xmatch && file::exists(file)) {
                struct {
                    vec1u idm;
                    uint_t ngal;
                } tmp;

                fits::read_table(file, ftable(tmp.idm, tmp.ngal));
                if (tmp.idm.size() != cra.size() || tmp.ngal != ngal) {
                    warning("incompatible cross-match data, re-matching");
                    rematch = true;
                }
            } else {
                rematch = true;
            }

            if (rematch) {
                qxmatch_res xm = qxmatch(cra, cdec, ra, dec,
                    keywords(_nth(2), _thread(4), _verbose(true))
                );

                const uint_t n = cra.size();
                vec1u idm(n);
                vec1u idn;

                for (uint_t i = 0; i < n; ++i) {
                    if (xm.rid[xm.id[i]] == i) {
                        idm[i] = xm.id[i];
                    } else {
                        idn.push_back(i);
                        idm[i] = ra.size();
                        ra.push_back(cra[i]);
                        dec.push_back(cdec[i]);
                        origin.push_back(pool.size()+1);
                    }
                }

                fits::write_table(file, ftable(xm.d, idm, idn, ngal));

                print("> ", n - idn.size(), " matched sources, ", idn.size(), " new sources");
                ngal += idn.size();

                std::string ref = "["+strn(pool.size()+1)+"]";
                pool.push_back({*this, idm, std::move(xm.d), name, sources, files, comment, ref});
                return pool.back();
            } else {
                print("reading data from "+file);

                vec1u idm, idn;
                vec2d d;

                fits::read_table(file, ftable(d, idm, idn));

                for (uint_t i = 0; i < idn.size(); ++i) {
                    ra.push_back(cra[idn[i]]);
                    dec.push_back(cdec[idn[i]]);
                    origin.push_back(pool.size()+1);
                }

                print("> ", cra.size() - idn.size(), " matched sources, ", idn.size(),
                    " new sources");

                ngal += idn.size();

                std::string ref = "["+strn(pool.size()+1)+"]";
                pool.push_back({*this, idm, std::move(d), name, sources, files, comment, ref});
                return pool.back();
            }
        }
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
    in[ix(idm,rep<Dim-1>(_))] = out;
}

template<typename T, typename U, typename V>
void catalog_t::merge(vec_t<1,T>& in, const U& out, const V& def) {
    phypp_check(in.empty(), name+": merging twice into the same variable");
    in = replicate(def, pool.ngal);
    in[idm] = out;
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

template<typename T, typename U, typename V>
void catalog_t::merge(T& in, const U& out, const V& def) {
    struct {
        reflex::struct_t<T> t;
        const V& def;
        catalog_t& cat;

        template<typename M>
        void operator () (const reflex::member_t& m, const M& v) {
            struct {
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
            } do_run{m, v, this->def, this->cat};

            reflex::foreach_member(this->t, do_run);
        }
    } do_merge{reflex::wrap(in), def, *this};

    reflex::foreach_member(reflex::wrap(out), do_merge);
}


template<typename T, typename U, typename V>
void catalog_t::merge(T& in, const U& out, const V& def, const std::string& com) {
    merge(in, out, def);
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
    uint_t width = naxes[0], height = naxes[1];

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
        // Discard any source that falls out of the boundaries of the image
        if (x[i]-hsize < 1 || x[i]+hsize >= width || y[i] - hsize < 1 || y[i] + hsize >= height) {
            continue;
        }

        vec_t<2,Type> cut(2*hsize+1, 2*hsize+1);

        Type null = fnan;
        int anynul = 0;
        long p0[2] = {long(round(x[i]-hsize)), long(round(y[i]-hsize))};
        long p1[2] = {long(round(x[i]+hsize)), long(round(y[i]+hsize))};
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
    uint_t width = naxes[0], height = naxes[1];
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
        // Discard any source that falls out of the boundaries of the image
        if (x[i]-hsize < 1 || x[i]+hsize >= width || y[i] - hsize < 1 || y[i] + hsize >= height) {
            continue;
        }

        vec_t<2,Type> cut(2*hsize+1, 2*hsize+1);
        vec_t<2,Type> wcut(2*hsize+1, 2*hsize+1);

        Type null = fnan;
        int anynul = 0;
        long p0[2] = {long(round(x[i]-hsize)), long(round(y[i]-hsize))};
        long p1[2] = {long(round(x[i]+hsize)), long(round(y[i]+hsize))};
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

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
auto field_area(const TX& ra, const TY& dec) {
    vec1u hull = convex_hull(ra, dec);

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
    // NOTE: this should be more correct, and will not miss any feature of both the
    // SED and the filter, but is much slower

    // vec1d xgrid = merge(filter.lam, lam);
    // xgrid = xgrid[sort(xgrid)];
    // vec1d ifil = interpolate(filter.res, filter.lam, xgrid);
    // double mx = min(filter.lam);
    // ifil[where(xgrid < mx)] = 0;
    // mx = max(filter.lam);
    // ifil[where(xgrid > mx)] = 0;
    // return integrate(xgrid, ifil*interpolate(sed, lam, xgrid))/integrate(xgrid, ifil);

    // NOTE: this is a faster approximation, where the SED is interpolated on the
    // response curve of the filter
    return integrate(filter.lam, filter.res*interpolate(sed, lam, filter.lam));
}

template<typename TypeL, typename TypeS>
auto sed2flux(const filter_t& filter, const vec_t<2,TypeL>& lam, const vec_t<2,TypeS>& sed) {
    using rtype = decltype(sed[0]*filter.res[0]);
    const uint_t nsed = sed.dims[0];
    vec_t<1,rtype> r; r.reserve(nsed);

    for (uint_t s = 0; s < nsed; ++s) {
        r.push_back(integrate(filter.lam, filter.res*interpolate(sed(s,_), lam(s,_), filter.lam)));
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

struct template_fit_res_t {
    uint_t bfit; // index of the best fit template in the library
    vec1d chi2;  // chi^2 of each template
    vec1d amp;   // renormalization amplitude of each template
    vec2d flux;  // best-fit flux for each template

    vec1u sed_sim;       // index of each error realization's best fit template
    vec1d amp_bfit_sim;  // renormalization amplitude of the best fit for each error realization
    vec1d amp_sim;       // renormalization amplitude for each error realization's best fit
};

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

template<typename TypeLib, typename TypeF, typename TypeE, typename TypeFi, typename TypeSeed,
    typename TypeZ, typename TypeD>
template_fit_res_t template_fit_renorm(const TypeLib& lib, TypeSeed& seed, const TypeZ& z,
   const TypeD& d, const vec_t<1,TypeF>& flux, const vec_t<1,TypeE>& err,
   const vec_t<1,TypeFi>& filters) {

    template_fit_res_t res;
    res.flux = template_observed(lib, z, d, filters);

    // Compute chi2 & renormalization factor
    auto weight = invsqr(err);
    using ttype = decltype(weight[0]*flux[0]*res.flux[0]);

    const uint_t nsed = res.flux.dims[0];
    const uint_t nfilter = filters.size();
    auto tmp1 = arr<ttype>(nsed);
    auto tmp2 = arr<ttype>(nsed);
    for (uint_t i = 0; i < nsed; ++i) {
        auto tmp = res.flux(i,_);
        tmp1[i] = total(weight*flux*tmp);
        tmp2[i] = total(weight*tmp*tmp);
    }

    res.amp = tmp1/tmp2;
    tmp1 *= res.amp;

    res.chi2 = total(weight*flux*flux) - tmp1;

    // Find the best chi2 among all the SEDs
    res.bfit = min_id(res.chi2);

    // Compute the errors on the fit by adding a random offset to the measured photometry according
    // to the provided error. The fit is performed on each of these random realizations and the
    // error on the parameters are computed as the standard deviation of the fit results over all
    // the realizations.
    const uint_t nsim = 200;
    res.amp_bfit_sim = fltarr(nsim);
    res.amp_sim = fltarr(nsim);
    res.sed_sim = fltarr(nsim);

    const uint_t nflux = flux.size();
    for (uint_t i = 0; i < nsim; ++i) {
        auto fsim = flux + randomn(seed, nflux)*err;
        for (uint_t t = 0; t < nsed; ++t) {
            tmp1[t] = total(weight*fsim*res.flux(t,_));
        }

        auto amp = tmp1/tmp2;
        tmp1 *= amp;

        auto chi2 = total(weight*fsim*fsim) - tmp1;
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

template<typename TypeLib, typename TypeF, typename TypeE, typename TypeFi, typename TypeSeed,
    typename TypeZ, typename TypeD>
template_fit_res_t template_fit(const TypeLib& lib, TypeSeed& seed, const TypeZ& z, const TypeD& d,
    const vec_t<1,TypeF>& flux, const vec_t<1,TypeE>& err, const vec_t<1,TypeFi>& filters) {

    template_fit_res_t res;
    res.flux = template_observed(lib, z, d, filters);

    // Compute chi2 & renormalization factor
    auto weight = invsqr(err);
    using ttype = decltype(weight[0]*flux[0]*res.flux[0]);

    uint_t nsed = res.flux.dims[0];
    const uint_t nfilter = filters.size();
    auto tmp1 = arr<ttype>(nsed);
    auto tmp2 = arr<ttype>(nsed);
    for (uint_t i = 0; i < nsed; ++i) {
        auto tmp = res.flux(i,_);
        tmp1[i] = total(weight*flux*tmp);
        tmp2[i] = total(weight*tmp*tmp);
    }

    res.amp = tmp1/tmp2;

    res.chi2 = total(weight*flux*flux) - 2.0*tmp1 + tmp2;

    // Find the best chi2 among all the SEDs
    res.bfit = min_id(res.chi2);

    // Compute the errors on the fit by adding a random offset to the measured photometry according
    // to the provided error. The fit is performed on each of these random realizations and the
    // error on the parameters are computed as the standard deviation of the fit results over all
    // the realizations.
    const uint_t nsim = 1000;
    res.amp_bfit_sim = fltarr(nsim);
    res.amp_sim = fltarr(nsim);
    res.sed_sim = fltarr(nsim);

    const uint_t nflux = flux.size();
    for (uint_t i = 0; i < nsim; ++i) {
        auto fsim = flux + randomn(seed, nflux)*err;
        for (uint_t t = 0; t < nsed; ++t) {
            tmp1[t] = total(weight*fsim*res.flux(t,_));
        }

        auto amp = tmp1/tmp2;
        auto chi2 = total(weight*fsim*fsim) - 2.0*tmp1 + tmp2;
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

#endif

