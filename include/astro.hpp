#ifndef ASTRO_HPP
#define ASTRO_HPP

#include <vec.hpp>
#include <math.hpp>
#include <print.hpp>
#include <image.hpp>
#include <thread.hpp>
#include <map>

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
    
    int_t cnt;
    vec1i idnan = where(!finite(icut) || !finite(ecut), cnt);
    if (cnt != 0) {
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
    return (1.0+z)*(2.99792458e5/cosmo.H0)*integrate([&](T t) { return pow(pow((1.0+t),3)*cosmo.wm + cosmo.wL, -0.5); }, 0, z);
}

// Lookback time [Gyr] as a function of redshift 'z'.
// There is no analytic form for this function, hence it must be numerically integrated.
// Note that this function is costly. If you need to compute it for many redshift values, you are
// advised to use the interpolated version below.
template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
T lookback_time(const T& z, const cosmo_t& cosmo) {
    if (z <= 0) return 0.0;
    return (3.09/(cosmo.H0*3.155e-3))*integrate([&](T t) { return pow(pow((1+t),3)*cosmo.wm + cosmo.wL + pow(1+t,2)*cosmo.wk, -0.5)/(1+t); }, 0, z);
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
            for (int_t i = 0; i < npt+1; ++i) { \
                td[i] = name(tz[i], args...); \
            } \
            vec_t<1,rtype> r = z; \
            for (auto& t : r) { \
                t = interpol_fast(td, tz, t); \
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
double angdist(double ra1, double dec1, double ra2, double dec2) {
    double sra = sin(0.5*(ra2 - ra1));
    double sde = sin(0.5*(dec2 - dec1));
    return 2.0*asin(sqrt(sde*sde + sra*sra*cos(dec2)*cos(dec1)));
}

struct qxmatch_res {
    vec2i id;
    vec2d d;
    vec2i rid;
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
    declare_keywords(_thread(1), _nth(1), _verbose(false), _self(false))) {

    phypp_check(n_elements(ra1) == n_elements(dec1), "qxmatch: not as many RA as there are Dec");
    phypp_check(n_elements(ra2) == n_elements(dec2), "qxmatch: not as many RA as there are Dec");

    const double d2r = 3.14159265359/180.0;
    auto dra1  = ra1*d2r;
    auto ddec1 = dec1*d2r;
    auto dcdec1 = cos(ddec1);
    auto dra2  = ra2*d2r;
    auto ddec2 = dec2*d2r;
    auto dcdec2 = cos(ddec2);

    const int_t n1 = n_elements(ra1); 
    const int_t n2 = n_elements(ra2);
    
    const int_t nth = get_keyword(_nth);
    const bool self = get_keyword(_self);

    qxmatch_res res;
    res.id = intarr(nth, n1)-1;
    res.d  = dblarr(nth, n1)+dinf;
    res.rid = intarr(nth, n2)-1;
    res.rd  = dblarr(nth, n2)+dinf;

    auto work = [&] (int_t i, int_t j) {
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
        if (sd < res.d(nth-1,i)) {
            res.id(nth-1,i) = j;
            res.d(nth-1,i) = sd;
            int_t k = nth-2;
            while (k >= 0 && res.d(k,i) > res.d(k+1,i)) {
                std::swap(res.d(k,i), res.d(k+1,i));
                std::swap(res.id(k,i), res.id(k+1,i));
                --k;
            }
        }
        if (sd < res.rd(nth-1,j)) {
            res.rid(nth-1,j) = i;
            res.rd(nth-1,j) = sd;
            int_t k = nth-2;
            while (k >= 0 && res.rd(k,j) > res.rd(k+1,j)) {
                std::swap(res.rd(k,j), res.rd(k+1,j));
                std::swap(res.rid(k,j), res.rid(k+1,j));
                --k;
            }
        }
    };
    
    int_t nthread = get_keyword(_thread);
    if (nthread <= 1) {
        // When using a single thread, all the work is done in the main thread
        auto p = progress_start(n1);
        for (int_t i = 0; i < n1; ++i) {
            for (int_t j = 0; j < n2; ++j) {
                if (self && i == j) continue;
                work(i,j);
            }
            
            if (get_keyword(_verbose)) print_progress(p, i);
        }
    } else {
        // When using more than one thread, the work load is evenly separated between all the
        // available threads, such that they should more or less all end at the same time.
        std::atomic<int_t> iter(0);
    
        // Create the thread pool and launch the threads
        auto pool = thread::pool(nthread);
        int_t total = 0;
        int_t total2 = 0;
        for (int_t t = 0; t < nthread; ++t) {
            int_t assigned = floor(n1/float(nthread));
            if (t == nthread-1) {
                assigned = n1 - total;
            }
            
            pool[t].start([&iter, &work, total, self, assigned, n2]() {
                for (int_t i = total; i < total+assigned; ++i) {
                    for (int_t j = 0; j < n2; ++j) {
                        if (self && i == j) continue;
                        work(i,j);
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
    }

    // Convert the distance estimator to a real distance
    res.d = 3600.0*(180.0/3.14159265359)*2*asin(sqrt(res.d));

    return res;
}

struct id_pair {
    vec1i id1, id2;
    vec1i lost;
};

id_pair xmatch_clean_best(const qxmatch_res& r) {
    const int_t ngal = dim(r.id)[1];
    
    id_pair c;
    c.id1.data.reserve(ngal);
    c.id2.data.reserve(ngal);
    c.lost.data.reserve(ngal/6);

    for (int_t i = 0; i < ngal; ++i) {
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

struct xmatch_merge_from {
    std::size_t ngal;
    id_pair id;

    xmatch_merge_from(std::size_t ngal_, const id_pair& id_) : ngal(ngal_), id(id_) {}

    template<std::size_t Dim, typename T, typename U, typename V, typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& t, const vec_t<Dim,U>& u, const V& def) {
        t.dims = u.dims;
        t.dims[Dim-1] = ngal;
        t.resize();
        t[_] = def;
        t[ix(rep<Dim-1>(_),id.id2)] = u[ix(rep<Dim-1>(_),id.id1)];
    }

    template<std::size_t Dim, typename T, typename U, typename V, typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& t, const vec_t<Dim,U>& u) {
        t.dims = u.dims;
        t.dims[Dim-1] = ngal;
        t.resize();
        t[ix(rep<Dim-1>(_),id.id2)] = u[ix(rep<Dim-1>(_),id.id1)];
    }

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& t, const U& u, const V& def) {
        t = replicate(def, ngal);
        t[id.id2] = u[id.id1];
    }

    template<typename T, typename U>
    void merge(vec_t<1,T>& t, const U& u) {
        t = arr<T>(ngal);
        t[id.id2] = u[id.id1];
    }
};

struct xmatch_merge_to {
    std::size_t ngal;
    id_pair id;

    xmatch_merge_to(std::size_t ngal_, const id_pair& id_) : ngal(ngal_), id(id_) {}

    template<std::size_t Dim, typename T, typename U, typename V, typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& t, const vec_t<Dim,U>& u, const V& def) {
        t.dims = u.dims;
        t.dims[Dim-1] = ngal;
        t.resize();
        t[_] = def;
        t[id.id1] = u[ix(rep<Dim-1>(_),id.id2)];
    }

    template<std::size_t Dim, typename T, typename U, typename V, typename enable = typename std::enable_if<Dim != 1>::type>
    void merge(vec_t<Dim,T>& t, const vec_t<Dim,U>& u) {
        t.dims = u.dims;
        t.dims[Dim-1] = ngal;
        t.resize();
        t[id.id1] = u[ix(rep<Dim-1>(_),id.id2)];
    }

    template<typename T, typename U, typename V>
    void merge(vec_t<1,T>& t, const U& u, const V& def) {
        t = arr<T>(ngal);
        t[_] = def;
        t[id.id1] = u[id.id2];
    }

    template<typename T, typename U>
    void merge(vec_t<1,T>& t, const U& u) {
        t = arr<T>(ngal);
        t[id.id1] = u[id.id2];
    }
};

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

struct catalog_t {
    std::string name;
    vec1s       sources;
    vec1s       files;
    std::string comment;
    std::string ref;
};

using catalog_pool_t = std::list<catalog_t>;

catalog_pool_t create_catalog_pool() {
    return catalog_pool_t();
}

catalog_t& add_catalog(catalog_pool_t& pool, const std::string& name, const vec1s& sources, 
    const vec1s& files, const std::string& comment = "") {

    std::string ref = "["+strn(pool.size()+1)+"]";
    pool.push_back({name, sources, files, comment, ref});
    return pool.back();
}

std::string build_catalog_list(const catalog_pool_t& pool, uint_t width = 80) {
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

// Compute the area covered by a field given a set of source coordinates [deg^2].
// Coordinates are assumed to be given in degrees.
template<typename TX, typename TY>
auto field_area(const TX& ra, const TY& dec) {
    vec1i hull = convex_hull(ra, dec);

    decltype(1.0*ra[0]*dec[0]) area = 0;

    auto d2r = dpi/180.0;
    typename TX::effective_type hx = ra[hull]*d2r;
    typename TY::effective_type hy = dec[hull]*d2r;
    uint_t nh = hull.size();
    for (uint_t i = 2; i < nh; ++i) {
        double e1 = angdist(hx[0],   hy[0],   hx[i-1], hy[i-1]);
        double e2 = angdist(hx[i-1], hy[i-1], hx[i],   hy[i]);
        double e3 = angdist(hx[i],   hy[i],   hx[0],   hy[0]);
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

#endif

