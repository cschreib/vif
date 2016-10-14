#ifndef PHYPP_IO_FITS_WCS_HPP
#define PHYPP_IO_FITS_WCS_HPP

#ifndef NO_WCSLIB
#include <wcslib/wcshdr.h>
#include <wcslib/wcserr.h>
#endif

#include "phypp/core/error.hpp"
#include "phypp/io/fits/fits.hpp"

namespace phypp {
namespace fits {
    struct make_wcs_header_params {
        // The pixel size in arcsec
        double pixel_scale = dnan;
        // The reference position
        double sky_ref_ra = dnan, sky_ref_dec = dnan;
        // The pixel corresponding to the reference position
        double pixel_ref_x = dnan, pixel_ref_y = dnan;
        // The number of pixels in X and Y axis
        uint_t dims_x = npos, dims_y = npos;
    };

    // Add WCS data to a FITS header, computed from a set of simple parameters.
    inline bool make_wcs_header(const make_wcs_header_params& params, fits::header& hdr) {
        if (hdr.empty()) {
            hdr = "END" + std::string(77, ' ');
        }

        if (is_finite(params.pixel_scale)) {
            if (!setkey(hdr, "CDELT1", -params.pixel_scale/3600.0)) {
                error("make_wcs_header: could not set keyword 'CDELT1' to '",
                    -params.pixel_scale, "'");
                return false;
            }
            if (!setkey(hdr, "CDELT2", params.pixel_scale/3600.0)) {
                error("make_wcs_header: could not set keyword 'CDELT2' to '",
                    params.pixel_scale, "'");
                return false;
            }
            if (!setkey(hdr, "CTYPE1", "'RA---TAN'")) {
                error("make_wcs_header: could not set keyword 'CTYPE1' to 'RA---TAN'");
                return false;
            }
            if (!setkey(hdr, "CTYPE2", "'DEC--TAN'")) {
                error("make_wcs_header: could not set keyword 'CTYPE2' to 'DEC--TAN'");
                return false;
            }
            if (!setkey(hdr, "EQUINOX", 2000.0)) {
                error("make_wcs_header: could not set keyword 'EQUINOX' to '",
                    2000.0, "'");
                return false;
            }
        }

        if (is_finite(params.pixel_ref_x) && is_finite(params.pixel_ref_y)) {
            if (!setkey(hdr, "CRPIX1", params.pixel_ref_x)) {
                error("make_wcs_header: could not set keyword 'CRPIX1' to '",
                    params.pixel_ref_x, "'");
                return false;
            }
            if (!setkey(hdr, "CRPIX2", params.pixel_ref_y)) {
                error("make_wcs_header: could not set keyword 'CRPIX2' to '",
                    params.pixel_ref_y, "'");
                return false;
            }
        }

        if (is_finite(params.sky_ref_ra) && is_finite(params.sky_ref_dec)) {
            if (!setkey(hdr, "CRVAL1", params.sky_ref_ra)) {
                error("make_wcs_header: could not set keyword 'CRVAL1' to '",
                    params.sky_ref_ra, "'");
                return false;
            }
            if (!setkey(hdr, "CRVAL2", params.sky_ref_dec)) {
                error("make_wcs_header: could not set keyword 'CRVAL2' to '",
                    params.sky_ref_dec, "'");
                return false;
            }
        }

        if (params.dims_x != npos && params.dims_y != npos) {
            if (!setkey(hdr, "META_0", 2u)) {
                error("make_wcs_header: could not set keyword 'META_0' to '", 2u, "'");
                return false;
            }
            if (!setkey(hdr, "META_1", params.dims_x)) {
                error("make_wcs_header: could not set keyword 'META_1' to '",
                    params.dims_x, "'");
                return false;
            }
            if (!setkey(hdr, "META_2", params.dims_y)) {
                error("make_wcs_header: could not set keyword 'META_2' to '",
                    params.dims_y, "'");
                return false;
            }
        }

        return true;
    }

    // Add WCS data to a FITS header, computed from a set of simple parameters.
    // Format: {"pixel_scale:0.06", "sky_ref:-3.56985,52.6456", ...}
    // Parameters:
    //  - pixel_scale [float]: the pixel size in arcsec
    //  - sky_ref [float,float]: the reference position
    //  - pixel_ref [float,float]: the pixel corresponding to the reference position
    //  - dims [uint,uint]: number of pixels in X and Y axis
    inline bool make_wcs_header(const vec1s& string_params, fits::header& hdr) {
        make_wcs_header_params params;

        for (auto& p : string_params) {
            vec1s spl = split(p, ":");

            if (spl.size() != 2) {
                error("make_wcs_header: parameter '", p, "' is ill formed");
                return false;
            }

            spl[0] = trim(tolower(spl[0]));

            if (spl[0] == "pixel_scale") {
                if (!from_string(spl[1], params.pixel_scale)) {
                    error("make_wcs_header: could not read pixel scale '",
                        spl[1], "' as double");
                    return false;
                }
            } else if (spl[0] == "pixel_ref") {
                vec1s tspl = split(spl[1], ",");
                if (tspl.size() != 2) {
                    error("make_wcs_header: ill formed 'pixel_ref' parameter: '", spl[0], "'");
                    note("make_wcs_header: expecting two comma separated coordinates of "
                        "reference pixel");
                    return false;
                }
                if (!from_string(tspl[0], params.pixel_ref_x)) {
                    error("make_wcs_header: could not read X pixel reference '",
                        tspl[0], "' as double");
                    return false;
                }
                if (!from_string(tspl[1], params.pixel_ref_y)) {
                    error("make_wcs_header: could not read Y pixel reference '",
                        tspl[1], "' as double");
                    return false;
                }
            } else if (spl[0] == "sky_ref") {
                vec1s tspl = split(spl[1], ",");
                if (tspl.size() != 2) {
                    error("make_wcs_header: ill formed 'sky_ref' parameter: '", spl[0], "'");
                    note("make_wcs_header: expecting two comma separated coordinates of "
                        "reference sky position");
                    return false;
                }
                if (!from_string(tspl[0], params.sky_ref_ra)) {
                    error("make_wcs_header: could not read RA sky position reference '",
                        tspl[0], "' as double");
                    return false;
                }
                if (!from_string(tspl[1], params.sky_ref_dec)) {
                    error("make_wcs_header: could not read Dec sky position reference '",
                        tspl[1], "' as double");
                    return false;
                }
            } else if (spl[0] == "dims") {
                vec1s tspl = split(spl[1], ",");
                if (tspl.size() != 2) {
                    error("make_wcs_header: ill formed 'dims' parameter: '", spl[0], "'");
                    note("make_wcs_header: expecting two comma separated number of pixels");
                    return false;
                }
                if (!from_string(tspl[0], params.dims_x)) {
                    error("make_wcs_header: could not read number of pixels in first axis '",
                        tspl[0], "' as unsigned integer");
                    return false;
                }
                if (!from_string(tspl[1], params.dims_y)) {
                    error("make_wcs_header: could not read number of pixels in second axis '",
                        tspl[1], "' as unsigned integer");
                    return false;
                }
            } else {
                error("make_wcs_header: unknown parameter '", spl[0], "'");
                return false;
            }
        }

        return make_wcs_header(params, hdr);
    }

    inline fits::header filter_wcs(const fits::header& hdr) {
        // List of keywords taken from 'cphead' (WCSTools)
        vec1s keywords = {"RA", "DEC", "EPOCH", "EQUINOX", "RADECSYS", "SECPIX", "IMWCS",
            "CD1_1", "CD1_2", "CD2_1", "CD2_2", "PC1_1", "PC1_2", "PC2_1", "PC2_2",
            "PC001001", "PC001002", "PC002001", "PC002002", "LATPOLE", "LONPOLE",
            "SECPIX", "CTYPE", "CRVAL", "CDELT", "CRPIX", "CROTA",
            "CUNIT", "CO1_", "CO2_", "PROJP", "PV1_", "PV2_", "END"};

        vec1s okeys = cut(hdr, 80);
        vec1s nkeys;
        for (auto& k : okeys) {
            for (auto& wk : keywords) {
                if (start_with(k, wk)) {
                    nkeys.push_back(k);
                }
            }
        }

        return collapse(nkeys);
    }

#ifndef NO_WCSLIB
    // Extract astrometry from a FITS image header
    struct wcs {
        wcsprm* w = nullptr;
        int nwcs  = 0;

        explicit wcs(uint_t naxis = 2) : w(new wcsprm), nwcs(1) {
            w->flag = -1;
            wcsini(true, naxis, w);
        }

        explicit wcs(const fits::header& hdr) {
            // Enable error reporting
            wcserr_enable(1);

            // Feed the header to WCSLib to extract the astrometric parameters
            int nreject = 0;
            int success = wcspih(
                const_cast<char*>(hdr.c_str()), hdr.size()/80 + 1,
                WCSHDR_all, 0, &nreject, &nwcs, &w
            );

            if ((success != 0 || nwcs == 0) && w) {
                wcsvfree(&nwcs, &w);
                w = nullptr;
            }

            // Try a dummy coordinate conversion to see if everything is recognized
            if (w) {
                double map[2] = {0.0, 0.0};
                double itmp[2];
                double phi, theta;
                double world[2];
                int status = 0;

                int ret = wcsp2s(w, 1, 2, map, itmp, &phi, &theta, world, &status);
                if (ret != 0) {
                    wcserr_prt(w->err, "error: ");
                    wcsvfree(&nwcs, &w);
                    w = nullptr;
                }
            }
        }

        wcs(const wcs&) = delete;
        wcs& operator = (const wcs&) = delete;
        wcs(wcs&& tw) {
            std::swap(w, tw.w);
            std::swap(nwcs, tw.nwcs);
        }
        wcs& operator = (wcs&& tw) {
            if (w) {
                wcsvfree(&nwcs, &w);
            }
            w = tw.w; tw.w = nullptr;
            nwcs = tw.nwcs; tw.nwcs = 0;
            return *this;
        }

        bool is_valid() const {
            return w != nullptr;
        }

        ~wcs() {
            if (w) {
                wcsvfree(&nwcs, &w);
                w = nullptr;
            }
        }
    };
#else
    struct wcs {
        template <typename T, typename ... Args>
        wcs(Args&&...) {
            static_assert(!std::is_same<T,T>::value, "WCS support is is disabled, "
                "please enable the WCSLib library to use this function");
        }
    };
#endif

    template<typename Dummy = void>
    fits::wcs extast(const fits::header& hdr) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        return fits::wcs(hdr);
#endif
    }

    template<std::size_t D = 1, typename T = double, typename U = double, typename V, typename W>
    void ad2xy(const fits::wcs& w, const vec<D,T>& ra, const vec<D,U>& dec,
        vec<D,V>& x, vec<D,W>& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(ra.size() == dec.size(), "RA and Dec arrays do not match sizes ("+
            strn(ra.size())+" vs "+strn(dec.size())+")");

        uint_t ngal = ra.size();
        if (ngal == 0) {
            x.clear(); y.clear();
            return;
        }

        uint_t naxis = w.w->naxis;
        vec2d world(ngal, naxis);
        for (uint_t i : range(ngal)) {
            world.safe(i,0) = ra.safe[i];
            world.safe(i,1) = dec.safe[i];
        }

        vec2d pos(ngal, naxis);
        std::vector<double> phi(ngal), theta(ngal);
        std::vector<double> itmp(naxis*ngal);
        std::vector<int>    stat(ngal);

        int status = wcss2p(w.w, ngal, naxis, world.data.data(), phi.data(), theta.data(),
            itmp.data(), pos.data.data(), stat.data());

        if (status != 0) {
            wcserr_prt(w.w->err, "error: ");
        }

        phypp_check(status == 0, "error in WCS conversion");

        x.dims = ra.dims; x.resize();
        y.dims = dec.dims; y.resize();

        for (uint_t i : range(ngal)) {
            x.safe[i] = pos.safe(i,0);
            y.safe[i] = pos.safe(i,1);
        }
#endif
    }

    template<std::size_t D = 1, typename T = double, typename U = double, typename V, typename W>
    void xy2ad(const fits::wcs& w, const vec<D,T>& x, const vec<D,U>& y,
        vec<D,V>& ra, vec<D,W>& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(x.size() == y.size(), "x and y arrays do not match sizes ("+
            strn(x.size())+" vs "+strn(y.size())+")");

        uint_t ngal = x.size();
        if (ngal == 0) {
            ra.clear(); dec.clear();
            return;
        }

        uint_t naxis = w.w->naxis;
        vec2d map(ngal, naxis);
        for (uint_t i : range(ngal)) {
            map.safe(i,0) = x.safe[i];
            map.safe(i,1) = y.safe[i];
        }

        vec2d world(ngal, naxis);
        std::vector<double> phi(ngal), theta(ngal);
        std::vector<double> itmp(naxis*ngal);
        std::vector<int>    stat(ngal);

        int status = wcsp2s(w.w, ngal, naxis, map.data.data(), itmp.data(),
            phi.data(), theta.data(), world.data.data(), stat.data());

        if (status != 0) {
            wcserr_prt(w.w->err, "error: ");
        }

        phypp_check(status == 0, "error in WCS conversion");

        ra.dims = x.dims; ra.resize();
        dec.dims = y.dims; dec.resize();

        for (uint_t i : range(ngal)) {
            ra.safe[i] = world.safe(i,0);
            dec.safe[i] = world.safe(i,1);
        }
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value &&
            !meta::is_vec<U>::value && !meta::is_vec<V>::value && !meta::is_vec<W>::value>::type>
    void ad2xy(const fits::wcs& w, const T& ra, const U& dec, V& x, W& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        phypp_check(w.is_valid(), "invalid WCS data");

        uint_t naxis = w.w->naxis;
        vec1d world(naxis);
        world.safe[0] = ra;
        world.safe[1] = dec;

        vec1d pos(naxis);

        double phi, theta;
        std::vector<double> itmp(naxis);
        int stat;

        int status = wcss2p(w.w, 1, naxis, world.data.data(), &phi, &theta, itmp.data(),
            pos.data.data(), &stat);

        if (status != 0) {
            wcserr_prt(w.w->err, "error: ");
        }

        phypp_check(status == 0, "error in WCS conversion");

        x = pos.safe[0];
        y = pos.safe[1];
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value &&
            !meta::is_vec<U>::value && !meta::is_vec<V>::value && !meta::is_vec<W>::value>::type>
    void xy2ad(const fits::wcs& w, const T& x, const U& y, V& ra, W& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        phypp_check(w.is_valid(), "invalid WCS data");

        uint_t naxis = w.w->naxis;
        vec1d map(naxis);
        map.safe[0] = x;
        map.safe[1] = y;

        vec1d world(naxis);

        double phi, theta;
        std::vector<double> itmp(naxis);
        int stat;

        int status = wcsp2s(w.w, 1, naxis, map.data.data(), itmp.data(), &phi, &theta,
            world.data.data(), &stat);

        if (status != 0) {
            wcserr_prt(w.w->err, "error: ");
        }

        phypp_check(status == 0, "error in WCS conversion");

        ra = world.safe[0];
        dec = world.safe[1];
#endif
    }

    // Obtain the pixel size of a given image in arsec/pixel.
    // Will fail (return false) if no WCS information is present in the image.
    template<typename Dummy = void>
    bool get_pixel_size(const fits::wcs& wcs, double& aspix) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        if (!wcs.is_valid()) {
            return false;
        }

        // Convert radius to number of pixels
        vec1d r, d;
        fits::xy2ad(wcs, {0, 1}, {0, 0}, r, d);
        aspix = angdist(r.safe[0], d.safe[0], r.safe[1], d.safe[1]);

        return true;
#endif
    }

    // Obtain the pixel size of a given image in arsec/pixel.
    // Will fail (return false) if no WCS information is present in the image.
    template<typename Dummy = void>
    bool get_pixel_size(const std::string& file, double& aspix) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        if (end_with(file, ".sectfits")) {
            vec1s sects = fits::read_sectfits(file);
            return get_pixel_size(sects[0], aspix);
        } else {
            fits::header hdr = fits::read_header(file);
            auto wcs = fits::wcs(hdr);
            if (!get_pixel_size(wcs, aspix)) {
                warning("could not extract WCS information");
                note("parsing '", file, "'");
                return false;
            }

            return true;
        }
#endif
    }
}
}

#endif

