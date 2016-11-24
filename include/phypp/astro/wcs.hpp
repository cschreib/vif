#ifndef PHYPP_IO_ASTRO_WCS_HPP
#define PHYPP_IO_ASTRO_WCS_HPP

#ifndef NO_WCSLIB
#include <wcslib/wcshdr.h>
#include <wcslib/wcserr.h>
#endif

#include "phypp/core/error.hpp"
#include "phypp/io/fits.hpp"
#include "phypp/astro/astro.hpp"

namespace phypp {
namespace astro {
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
            if (!fits::setkey(hdr, "CDELT1", -params.pixel_scale/3600.0)) {
                error("make_wcs_header: could not set keyword 'CDELT1' to '",
                    -params.pixel_scale, "'");
                return false;
            }
            if (!fits::setkey(hdr, "CDELT2", params.pixel_scale/3600.0)) {
                error("make_wcs_header: could not set keyword 'CDELT2' to '",
                    params.pixel_scale, "'");
                return false;
            }
            if (!fits::setkey(hdr, "CTYPE1", "'RA---TAN'")) {
                error("make_wcs_header: could not set keyword 'CTYPE1' to 'RA---TAN'");
                return false;
            }
            if (!fits::setkey(hdr, "CTYPE2", "'DEC--TAN'")) {
                error("make_wcs_header: could not set keyword 'CTYPE2' to 'DEC--TAN'");
                return false;
            }
            if (!fits::setkey(hdr, "EQUINOX", 2000.0)) {
                error("make_wcs_header: could not set keyword 'EQUINOX' to '",
                    2000.0, "'");
                return false;
            }
        }

        if (is_finite(params.pixel_ref_x) && is_finite(params.pixel_ref_y)) {
            if (!fits::setkey(hdr, "CRPIX1", params.pixel_ref_x)) {
                error("make_wcs_header: could not set keyword 'CRPIX1' to '",
                    params.pixel_ref_x, "'");
                return false;
            }
            if (!fits::setkey(hdr, "CRPIX2", params.pixel_ref_y)) {
                error("make_wcs_header: could not set keyword 'CRPIX2' to '",
                    params.pixel_ref_y, "'");
                return false;
            }
        }

        if (is_finite(params.sky_ref_ra) && is_finite(params.sky_ref_dec)) {
            if (!fits::setkey(hdr, "CRVAL1", params.sky_ref_ra)) {
                error("make_wcs_header: could not set keyword 'CRVAL1' to '",
                    params.sky_ref_ra, "'");
                return false;
            }
            if (!fits::setkey(hdr, "CRVAL2", params.sky_ref_dec)) {
                error("make_wcs_header: could not set keyword 'CRVAL2' to '",
                    params.sky_ref_dec, "'");
                return false;
            }
        }

        if (params.dims_x != npos && params.dims_y != npos) {
            if (!fits::setkey(hdr, "META_0", 2u)) {
                error("make_wcs_header: could not set keyword 'META_0' to '", 2u, "'");
                return false;
            }
            if (!fits::setkey(hdr, "META_1", params.dims_x)) {
                error("make_wcs_header: could not set keyword 'META_1' to '",
                    params.dims_x, "'");
                return false;
            }
            if (!fits::setkey(hdr, "META_2", params.dims_y)) {
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

        vec1u dims;

        uint_t ra_axis = 0, dec_axis = 1;
        uint_t x_axis = 0, y_axis = 1;

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

            if (w) {
                // Get dimensions from the FITS header
                uint_t naxis = axis_count();
                dims.resize(naxis);
                for (uint_t i : range(dims)) {
                    uint_t dim = npos;
                    if (fits::getkey(hdr, "NAXIS"+strn(i+1), dim)) {
                        dims[dims.size()-1-i] = dim;
                    }
                }

                // Identify RA and Dec axis
                ra_axis = x_axis = find_axis("RA");
                dec_axis = y_axis = find_axis("DEC");

                // Y is by definition the first axis, so swap them if
                // the input file has DEC/RA instead of RA/DEC
                if (x_axis < y_axis) std::swap(x_axis, y_axis);

                // Try a dummy coordinate conversion to see if everything is recognized
                vec1d map = replicate(0.0, w->naxis);
                vec1d world(w->naxis);
                vec1d itmp(w->naxis);
                double phi, theta;
                int status = 0;

                int ret = wcsp2s(w, 1, w->naxis, map.data.data(), itmp.data.data(), &phi, &theta,
                    world.data.data(), &status);

                if (ret != 0) {
                    wcserr_prt(w->err, "error: ");
                    wcsvfree(&nwcs, &w);
                    w = nullptr;
                }
            }
        }

        uint_t axis_count() const {
            return (w ? w->naxis : 0);
        }

        uint_t find_axis(const std::string& name) const {
            for (uint_t i : range(axis_count())) {
                std::string ctype = split(w->ctype[i], "-")[0];
                if (ctype == toupper(name)) {
                    return w->naxis-1 - i;
                }
            }

            return npos;
        }

        wcs(const wcs&) = delete;
        wcs& operator = (const wcs&) = delete;

        wcs(wcs&& tw) {
            std::swap(w, tw.w);
            std::swap(nwcs, tw.nwcs);
            std::swap(dims, tw.dims);
            std::swap(ra_axis, tw.ra_axis);
            std::swap(dec_axis, tw.dec_axis);
            std::swap(x_axis, tw.x_axis);
            std::swap(y_axis, tw.y_axis);
        }

        wcs& operator = (wcs&& tw) {
            if (w) {
                wcsvfree(&nwcs, &w);
            }

            w = tw.w; tw.w = nullptr;
            nwcs = tw.nwcs; tw.nwcs = 0;
            dims = tw.dims; tw.dims.clear();
            ra_axis = tw.ra_axis;
            dec_axis = tw.dec_axis;
            x_axis = tw.x_axis;
            y_axis = tw.y_axis;

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
    astro::wcs extast(const fits::header& hdr) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        return astro::wcs(hdr);
#endif
    }

}

namespace impl {
    namespace wcs_impl {
        vec2d world2pix(const astro::wcs& w, const vec2d& world) {
            vec2d pix(world.dims);

            uint_t npt = world.dims[0];
            uint_t naxis = world.dims[1];

            std::vector<double> phi(npt), theta(npt);
            std::vector<double> itmp(naxis*npt);
            std::vector<int>    stat(npt);

            int status = wcss2p(w.w, npt, naxis, world.data.data(), phi.data(), theta.data(),
                itmp.data(), pix.data.data(), stat.data());

            if (status != 0) {
                wcserr_prt(w.w->err, "error: ");
            }

            phypp_check(status == 0, "error in WCS conversion");

            return pix;
        }

        vec2d pix2world(const astro::wcs& w, vec2d pix) {
            vec2d world(pix.dims);

            uint_t npt = pix.dims[0];
            uint_t naxis = pix.dims[1];

            std::vector<double> phi(npt), theta(npt);
            std::vector<double> itmp(naxis*npt);
            std::vector<int>    stat(npt);

            int status = wcsp2s(w.w, npt, naxis, pix.data.data(), itmp.data(),
                phi.data(), theta.data(), world.data.data(), stat.data());

            if (status != 0) {
                wcserr_prt(w.w->err, "error: ");
            }

            phypp_check(status == 0, "error in WCS conversion");

            return world;
        }
    }
}

namespace astro {

    template<std::size_t D = 1, typename T = double, typename U = double, typename V, typename W>
    void ad2xy(const astro::wcs& w, const vec<D,T>& ra, const vec<D,U>& dec,
        vec<D,V>& x, vec<D,W>& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(ra.dims == dec.dims, "RA and Dec arrays do not match sizes ("+
            strn(ra.dims)+" vs "+strn(dec.dims)+")");

        uint_t npt = ra.size();
        if (npt == 0) {
            x.clear(); y.clear();
            return;
        }

        uint_t naxis = w.axis_count();
        vec2d world(npt, naxis);
        for (uint_t i : range(npt)) {
            world.safe(i,naxis-1-w.ra_axis) = ra.safe[i];
            world.safe(i,naxis-1-w.dec_axis) = dec.safe[i];
        }

        vec2d pix = impl::wcs_impl::world2pix(w, world);

        x.resize(ra.dims);
        y.resize(ra.dims);

        for (uint_t i : range(npt)) {
            x.safe[i] = pix.safe(i,naxis-1-w.x_axis);
            y.safe[i] = pix.safe(i,naxis-1-w.y_axis);
        }
#endif
    }

    template<std::size_t D = 1, typename T = double, typename U = double, typename V, typename W>
    void xy2ad(const astro::wcs& w, const vec<D,T>& x, const vec<D,U>& y,
        vec<D,V>& ra, vec<D,W>& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(x.dims == y.dims, "x and y arrays do not match sizes ("+
            strn(x.dims)+" vs "+strn(y.dims)+")");

        uint_t npt = x.size();
        if (npt == 0) {
            ra.clear(); dec.clear();
            return;
        }

        uint_t naxis = w.axis_count();
        vec2d pix(npt, naxis);
        for (uint_t i : range(npt)) {
            pix.safe(i,naxis-1-w.x_axis) = x.safe[i];
            pix.safe(i,naxis-1-w.y_axis) = y.safe[i];
        }

        vec2d world = impl::wcs_impl::pix2world(w, pix);

        ra.resize(x.dims);
        dec.resize(x.dims);

        for (uint_t i : range(npt)) {
            ra.safe[i] = world.safe(i,naxis-1-w.ra_axis);
            dec.safe[i] = world.safe(i,naxis-1-w.dec_axis);
        }
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value &&
            !meta::is_vec<U>::value && !meta::is_vec<V>::value && !meta::is_vec<W>::value>::type>
    void ad2xy(const astro::wcs& w, const T& ra, const U& dec, V& x, W& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        vec<1,T> tra  = replicate(ra,  1);
        vec<1,T> tdec = replicate(dec, 1);
        vec<1,V> tx;
        vec<1,W> ty;

        ad2xy(w, tra, tdec, tx, ty);

        x = tx.safe[0];
        y = ty.safe[0];
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!meta::is_vec<T>::value &&
            !meta::is_vec<U>::value && !meta::is_vec<V>::value && !meta::is_vec<W>::value>::type>
    void xy2ad(const astro::wcs& w, const T& x, const U& y, V& ra, W& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        vec<1,T> tx = replicate(x, 1);
        vec<1,T> ty = replicate(y, 1);
        vec<1,V> tra;
        vec<1,W> tdec;

        xy2ad(w, tx, ty, tra, tdec);

        ra = tra.safe[0];
        dec = tdec.safe[0];
#endif
    }

    // Obtain the pixel size of a given image in arsec/pixel.
    // Will fail (return false) if no WCS information is present in the image.
    template<typename Dummy = void>
    bool get_pixel_size(const astro::wcs& wcs, double& aspix) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        if (!wcs.is_valid()) {
            return false;
        }

        // Convert radius to number of pixels
        vec1d r, d;
        astro::xy2ad(wcs, {0, 1}, {0, 0}, r, d);
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
            auto wcs = astro::wcs(hdr);
            if (!get_pixel_size(wcs, aspix)) {
                warning("could not extract WCS information");
                note("parsing '", file, "'");
                return false;
            }

            return true;
        }
#endif
    }

    enum class axis_unit {
        wcslib_default,

        wave_m,
        wave_cm,
        wave_mm,
        wave_um,
        wave_nm,
        wave_Angstrom,

        freq_Hz,
        freq_kHz,
        freq_MHz,
        freq_GHz,
        freq_THz,

        sky_deg,
        sky_rad
    };

}

namespace impl {
    namespace wcs_impl {
        double conv_si2unit(astro::axis_unit unit) {
            switch (unit) {
                // WCSlib returns data in SI units or degrees
                case astro::axis_unit::wcslib_default: return 1.0;

                case astro::axis_unit::wave_m:         return 1.0;
                case astro::axis_unit::wave_cm:        return 1e2;
                case astro::axis_unit::wave_mm:        return 1e3;
                case astro::axis_unit::wave_um:        return 1e6;
                case astro::axis_unit::wave_nm:        return 1e9;
                case astro::axis_unit::wave_Angstrom:  return 1e10;

                case astro::axis_unit::freq_Hz:        return 1.0;
                case astro::axis_unit::freq_kHz:       return 1e-3;
                case astro::axis_unit::freq_MHz:       return 1e-6;
                case astro::axis_unit::freq_GHz:       return 1e-9;
                case astro::axis_unit::freq_THz:       return 1e-12;

                case astro::axis_unit::sky_deg:
                case astro::axis_unit::sky_rad:        return dpi/180.0;
            }
        }

        template <std::size_t D, typename T>
        void si2unit(vec<D,T>& data, astro::axis_unit unit) {
            double conv = conv_si2unit(unit);
            if (conv != 1.0) {
                data *= conv;
            }
        }

        template <std::size_t D, typename T>
        void unit2si(vec<D,T>& data, astro::axis_unit unit) {
            double conv = conv_si2unit(unit);
            if (conv != 1.0) {
                data /= conv;
            }
        }
    }
}

namespace astro {

    template<std::size_t D = 1, typename T = double, typename U>
    void x2w(const astro::wcs& wcs, uint_t axis, const vec<D,T>& x, vec<D,U>& w,
        axis_unit unit = axis_unit::wcslib_default) {

#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        uint_t naxis = wcs.axis_count();

        phypp_check(wcs.is_valid(), "invalid WCS data");
        phypp_check(axis < naxis, "trying to use an axis that does not exist (",
            axis, " vs ", naxis, ")");

        uint_t npix = x.size();

        vec2d pix(npix, naxis);
        pix.safe(_,naxis-1-axis) = x;
        for (uint_t i : range(naxis)) {
            if (i == axis) continue;
            pix(_,naxis-1-i) = wcs.dims[i]/2.0 + 1.0;
        }

        vec2d world = impl::wcs_impl::pix2world(wcs, pix);
        w = world.safe(_,naxis-1-axis);
        impl::wcs_impl::si2unit(w, unit);
#endif
    }

    template<std::size_t D = 1, typename T = double, typename U>
    void w2x(const astro::wcs& wcs, uint_t axis, const vec<D,T>& w, vec<D,U>& x,
        axis_unit unit = axis_unit::wcslib_default) {

#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        uint_t naxis = wcs.axis_count();

        phypp_check(wcs.is_valid(), "invalid WCS data");
        phypp_check(axis < naxis, "trying to use an axis that does not exist (",
            axis, " vs ", naxis, ")");

        uint_t npix = w.size();

        vec1d tw = w;
        impl::wcs_impl::unit2si(tw, unit);

        vec2d world(npix, naxis);
        world.safe(_,naxis-1-axis) = tw;

        for (uint_t i : range(naxis)) {
            if (i == axis) continue;
            world(_,naxis-1-i) = wcs.dims[i]/2.0 + 1.0;
        }

        vec2d pix = impl::wcs_impl::world2pix(wcs, world);
        x = pix.safe(_,naxis-1-axis);
#endif
    }

    template <typename Dummy = void>
    void x2w(const astro::wcs& wcs, uint_t axis, double x, double& w,
        axis_unit unit = axis_unit::wcslib_default) {

#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        vec1d tx = replicate(x, 1);
        vec1d tw;

        x2w(wcs, axis, tx, tw, unit);

        w = tw.safe[0];
#endif
    }

    template <typename Dummy = void>
    void w2x(const astro::wcs& wcs, uint_t axis, double w, double& x,
        axis_unit unit = axis_unit::wcslib_default) {

#ifdef NO_WCSLIB
        static_assert(!std::is_same<Dummy,Dummy>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        vec1d tw = replicate(w, 1);
        vec1d tx;

        w2x(wcs, axis, tw, tx, unit);

        x = tw.safe[0];
#endif
    }

    // Obtain the pixel size of a given image in arsec/pixel.
    // Will fail (return false) if no WCS information is present in the image.
    template<typename Type = double>
    vec<1,Type> build_axis(const astro::wcs& wcs, uint_t axis, axis_unit unit = axis_unit::wcslib_default) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<Type,Type>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        uint_t naxis = wcs.axis_count();

        phypp_check(wcs.is_valid(), "invalid WCS data");
        phypp_check(axis < naxis, "trying to use an axis that does not exist (",
            axis, " vs ", naxis, ")");

        uint_t npix = wcs.dims[axis];

        vec2d pix(npix, naxis);
        pix(_,naxis-1-axis) = dindgen(npix)+1;
        for (uint_t i : range(naxis)) {
            if (i == axis) continue;
            pix(_,naxis-1-i) = wcs.dims[i]/2.0 + 1.0;
        }

        vec2d world = impl::wcs_impl::pix2world(wcs, pix);

        vec<1,Type> ret = world(_,naxis-1-axis);
        impl::wcs_impl::si2unit(ret, unit);

        return ret;
#endif
    }
}
}

#endif

