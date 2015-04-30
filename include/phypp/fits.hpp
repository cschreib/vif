#ifndef FITS_HPP
#define FITS_HPP

#ifndef NO_CCFITS
#include <CCfits/CCfits>
#endif
#ifndef NO_WCSLIB
#include <wcslib/wcshdr.h>
#include <wcslib/wcserr.h>
#endif
#include <fitsio.h>
#include <string>
#include "phypp/vec.hpp"
#include "phypp/reflex.hpp"
#include "phypp/file.hpp"
#include "phypp/print.hpp"
#include "phypp/math.hpp"

namespace fits {
#ifndef NO_CCFITS
    using namespace CCfits;
#endif

    struct exception : std::exception {
        explicit exception(const std::string& m) : msg(m) {}
        std::string msg;

        const char* what() const noexcept override {
            return msg.c_str();
        }
    };

    template<typename T>
    struct traits;

    template<>
    struct traits<std::string> {
        using dtype = char;

        static const char tform = 'B';
        static const int ttype = TBYTE;

        static std::string def() {
            return "";
        }

        static bool is_convertible(int type) {
            if (type == TSTRING) return true;
            if (type == TBYTE) return true;
            return false;
        }
    };

    template<>
    struct traits<bool> {
        using dtype = char;

        static const char tform = 'B';
        static const int ttype = TBYTE;
        static const int image_type = BYTE_IMG;

        static bool def() {
            return false;
        }

        static bool is_convertible(int type) {
            if (type == TLOGICAL) return true;
            if (type == TBIT) return true;
            if (type == TBYTE) return true;
            return false;
        }
    };

    template<>
    struct traits<char> {
        using dtype = char;

        static const char tform = 'S';
        static const int ttype = TSBYTE;
        static const int image_type = SBYTE_IMG;

        static bool def() {
            return '\0';
        }

        static bool is_convertible(int type) {
            if (type == TBYTE) return true;
            return false;
        }
    };

    template<>
    struct traits<uint_t> {
        using dtype = uint_t;

        static const char tform = 'J';
        static const int ttype = TLONG;
        static const int image_type = LONG_IMG;

        static uint_t def() {
            return 0;
        }

        static bool is_convertible(int type) {
            if (type == TSHORT) return true;
            if (type == TLONG) return true;
            if (type == TBIT) return true;
            if (type == TBYTE) return true;
            if (type == TINT32BIT) return true;
            return false;
        }
    };

    template<>
    struct traits<int_t> {
        using dtype = int_t;

        static const char tform = 'J';
        static const int ttype = TLONG;
        static const int image_type = LONG_IMG;

        static int_t def() {
            return 0;
        }

        static bool is_convertible(int type) {
            if (type == TSHORT) return true;
            if (type == TLONG) return true;
            if (type == TBIT) return true;
            if (type == TBYTE) return true;
            if (type == TINT32BIT) return true;
            return false;
        }
    };

    template<>
    struct traits<float> {
        using dtype = float;

        static const char tform = 'E';
        static const int ttype = TFLOAT;
        static const int image_type = FLOAT_IMG;

        static float def() {
            return fnan;
        }

        static bool is_convertible(int type) {
            if (type == TSHORT) return true;
            if (type == TLONG) return true;
            if (type == TFLOAT) return true;
            if (type == TBIT) return true;
            if (type == TBYTE) return true;
            if (type == TINT32BIT) return true;
            return false;
        }
    };

    template<>
    struct traits<double> {
        using dtype = double;

        static const char tform = 'D';
        static const int ttype = TDOUBLE;
        static const int image_type = DOUBLE_IMG;

        static float def() {
            return dnan;
        }

        static bool is_convertible(int type) {
            if (type == TSHORT) return true;
            if (type == TLONG) return true;
            if (type == TFLOAT) return true;
            if (type == TDOUBLE) return true;
            if (type == TBIT) return true;
            if (type == TBYTE) return true;
            if (type == TINT32BIT) return true;
            return false;
        }
    };

    #define phypp_check_fits(assertion, msg) \
        do { if (!(assertion)) throw fits::exception(msg); } while(0)

    void phypp_check_cfitsio(int status, const std::string& msg) {
        if (status != 0) {
            char txt[FLEN_STATUS];
            fits_get_errstatus(status, txt);
            phypp_check_fits(status == 0, "error: cfitsio: "+std::string(txt)+"\nerror: "+msg);
        }
    }

    using header = std::string;

    // Return the number of dimensions of a FITS file
    template<typename Dummy = void>
    uint_t file_axes(const std::string& name) {
#ifdef NO_CCFITS
        static_assert(!std::is_same<Dummy,Dummy>::value, "CCfits is disabled, "
            "please enable the CCfits library to use this function");
#else
        try {
            return fits::FITS(name, fits::Read).pHDU().axes();
        } catch (fits::FitsException& e) {
            error("reading: "+name);
            error("FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
#endif
    }

    // Return the dimensions of a FITS image
    // Note: will return an empty array for FITS tables
    vec1u file_dimensions(const std::string& filename) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_image(&fptr, filename.c_str(), READONLY, &status);
            fits::phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            vec1u dims;
            int naxis = 0;
            fits_get_img_dim(fptr, &naxis, &status);
            if (naxis != 0) {
                std::vector<long> naxes(naxis);
                fits_get_img_size(fptr, naxis, naxes.data(), &status);
                for (auto& l : naxes) {
                    dims.push_back(l);
                }
            }

            fits_close_file(fptr, &status);
            return dims;
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename Dummy = void>
    bool is_cube(const std::string& name) {
        return file_axes<Dummy>(name) == 3;
    }

    template<typename Dummy = void>
    bool is_image(const std::string& name) {
        return file_axes<Dummy>(name) == 2;
    }

#ifndef NO_CCFITS
    template<typename THDU, std::size_t Dim, typename Type>
    void read_impl__(THDU& hdu, vec<Dim, Type>& v, std::string& hdr) {
        phypp_check_fits(hdu.axes() == Dim, "FITS file does not match array dimensions ("+
            strn(hdu.axes())+" vs "+strn(Dim)+")");

        std::size_t n = 1;
        for (uint_t j = 0; j < Dim; ++j) {
            v.dims[j] = hdu.axis(Dim-1-j);
            n *= v.dims[j];
        }

        std::valarray<Type> tv(n);
        hdu.read(tv, 1, n);

        v.data.assign(std::begin(tv), std::end(tv));

        // Read the header as a string
        char* hstr = nullptr;
        int nkeys  = 0;
        int status = 0;
        fits_hdr2str(hdu.fitsPointer(), 0, nullptr, 0, &hstr, &nkeys, &status);
        hdr = hstr;
        free(hstr);
    }
#endif

    // Load the content of a FITS file into an array.
    template<std::size_t Dim, typename Type>
    void read(const std::string& name, vec<Dim, Type>& v, fits::header& hdr) {
#ifdef NO_CCFITS
        static_assert(!std::is_same<Type,Type>::value, "CCfits is disabled, "
            "please enable the CCfits library to use this function");
#else
        try {
            fits::FITS file(name, fits::Read);
            fits::PHDU& phdu = file.pHDU();
            if (phdu.axes() == 0) {
                fits::ExtHDU* hdu = nullptr;
                uint_t i = 1;
                while (!hdu && i <= file.extension().size()) {
                    hdu = &file.extension(i);
                    if (hdu->axes() == 0) hdu = nullptr;
                }

                phypp_check_fits(hdu, "FITS file does not contain any data");

                read_impl__(*hdu, v, hdr);
            } else {
                read_impl__(phdu, v, hdr);
            }
        } catch (fits::FitsException& e) {
            error("reading: "+name);
            error("FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
#endif
    }

    template<std::size_t Dim, typename Type>
    void read(const std::string& name, vec<Dim, Type>& v) {
        std::string hdr;
        read(name, v, hdr);
    }

    // Load the content of a FITS file into an array.
    template<std::size_t Dim, typename Type>
    void read_hdu(const std::string& name, vec<Dim, Type>& v, uint_t hdu, fits::header& hdr) {
#ifdef NO_CCFITS
        static_assert(!std::is_same<Type,Type>::value, "CCfits is disabled, "
            "please enable the CCfits library to use this function");
#else
        try {
            fits::FITS file(name, fits::Read);
            fits::ExtHDU& thdu = file.extension(hdu+1);
            phypp_check_fits(thdu.axes() != 0, "FITS HDU does not contain any data");
            read_impl__(thdu, v, hdr);
        } catch (fits::FitsException& e) {
            error("reading: "+name);
            error("FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
#endif
    }

    template<std::size_t Dim, typename Type>
    void read_hdu(const std::string& name, vec<Dim, Type>& v, uint_t hdu) {
        std::string hdr;
        read_hdu(name, v, hdu, hdr);
    }

    template<std::size_t Dim = 2, typename Type = double>
    vec<Dim, Type> read(const std::string& name) {
        vec<Dim, Type> v;
        read(name, v);
        return v;
    }

    vec1s read_sectfits(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw fits::exception("could not find file '"+filename+"'");
        }

        std::string dir = file::get_directory(filename);

        vec1s files;
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '#') {
                files.push_back(dir+line);
            }
        }

        if (files.empty()) {
            throw fits::exception("empty sectfits '"+filename+"'");
        }

        return files;
    }

    fits::header read_header(const std::string& filename) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_image(&fptr, filename.c_str(), READONLY, &status);
            fits::phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            int naxis = 0;
            fits_get_img_dim(fptr, &naxis, &status);
            if (naxis == 0) {
                fits_close_file(fptr, &status);
                fits_open_table(&fptr, filename.c_str(), READONLY, &status);
                fits::phypp_check_cfitsio(status, "cannot open file '"+filename+"'");
            }

            // Read the header as a string
            char* hstr = nullptr;
            int nkeys  = 0;
            fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
            fits::header hdr = hstr;
            free(hstr);

            fits_close_file(fptr, &status);
            return hdr;
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    fits::header read_header_hdu(const std::string& filename, uint_t hdu) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_file(&fptr, filename.c_str(), READONLY, &status);
            fits::phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            int nhdu = 0;
            fits_get_num_hdus(fptr, &nhdu, &status);
            phypp_check(hdu < uint_t(nhdu), "requested HDU does not exists in this FITS file "
                "(", hdu, " vs. ", nhdu, ")");

            fits_movabs_hdu(fptr, hdu+1, nullptr, &status);

            // Read the header as a string
            char* hstr = nullptr;
            int nkeys  = 0;
            fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
            fits::header hdr = hstr;
            free(hstr);

            fits_close_file(fptr, &status);
            return hdr;
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    fits::header read_header_sectfits(const std::string& filename, uint_t sect) {
        if (end_with(filename, ".sectfits")) {
            vec1s sects = read_sectfits(filename);
            phypp_check(sect < sects.size(), "no section ", sect, " in '", filename,
                "' (only ", sects.size()," available)");

            return read_header(sects[sect]);
        } else {
            return read_header(filename);
        }
    }

    template<typename T>
    bool getkey(const fits::header& hdr, const std::string& key, T& v) {
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i = 0; i < nentry; ++i) {
            std::string entry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
            std::size_t eqpos = entry.find_first_of("=");
            if (eqpos == entry.npos) continue;
            std::size_t cpos = entry.find_first_of("/", eqpos);

            std::string nam = trim(entry.substr(0, eqpos));

            if (nam == key) {
                std::string value = trim(entry.substr(eqpos+1, cpos-eqpos-1));
                return from_string(value, v);
            }
        }

        return false;
    }

    bool getkey(const fits::header& hdr, const std::string& key, std::string& v) {
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i = 0; i < nentry; ++i) {
            std::string entry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
            std::size_t eqpos = entry.find_first_of("=");
            if (eqpos == entry.npos) continue;
            std::size_t cpos = entry.find_first_of("/", eqpos);

            std::string nam = trim(entry.substr(0, eqpos));

            if (nam == key) {
                std::string value = trim(entry.substr(eqpos+1, cpos-eqpos-1));
                if (value.find_first_of("'") == 0) {
                    value = trim(value, "'");
                }
                v = value;
                return true;
            }
        }

        return false;
    }

    template<typename T>
    bool setkey(fits::header& hdr, const std::string& key, const T& v,
        const std::string& comment = "") {

        if (key.size() > 8) return false;

        // Build new entry
        std::string entry = key+std::string(8-key.size(), ' ')+"= ";
        std::string value = strn(v);
        if (!comment.empty()) {
            value += " / "+comment;
        }

        entry += value;
        if (entry.size() > 80) return false;

        // Insertion point
        std::size_t ipos = hdr.npos;
        if (hdr.size() >= 80) {
            ipos = hdr.substr(hdr.size() - 80, 80).find("END");
            if (ipos != hdr.npos) ipos += hdr.size() - 80;
        }

        // Look for existing entry in header
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i = 0; i < nentry; ++i) {
            // Extract entry 'i'
            std::string tentry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
            std::size_t eqpos = tentry.find_first_of("=");
            if (eqpos == tentry.npos) continue;

            // Extract entry's name
            std::string nam = trim(tentry.substr(0, eqpos));
            // Look for comments
            std::size_t cpos = tentry.find_first_of("/", eqpos);

            if (nam == key) {
                // The entry already exists
                ipos = i*80;

                // Copy comments of the previous entry
                if (comment.empty() && cpos != tentry.npos) {
                    std::string cmt = tentry.substr(cpos+1);
                    entry += " /"+cmt;
                    if (entry.size() > 80) {
                        if (entry.back() == '&') {
                            entry.resize(80);
                            entry.back() = '&';
                        } else {
                            entry.resize(80);
                        }
                    }
                }

                hdr.erase(ipos, 80);

                break;
            }
        }

        entry += std::string(80-entry.size(), ' ');
        if (ipos == hdr.npos) {
            hdr += entry;
        } else {
            hdr.insert(ipos, entry);
        }

        return true;
    }

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
    bool make_wcs_header(const make_wcs_header_params& params, fits::header& hdr) {
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
    bool make_wcs_header(const vec1s& string_params, fits::header& hdr) {
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

    fits::header filter_wcs(const fits::header& hdr) {
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

        wcs() = default;

        explicit wcs(const fits::header& hdr) {
            // Feed it to WCSLib to extract the astrometric parameters
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

                wcserr_enable(1);
                int ret = wcsp2s(w, 1, 2, map, itmp, &phi, &theta, world, &status);
                if (ret != 0) {
                    wcserr_prt(w->err, "error: ");
                    wcsvfree(&nwcs, &w);
                    w = nullptr;
                }
                wcserr_enable(0);
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
    struct wcs {};
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

    template<typename T = double, typename U = double, typename V, typename W>
    void ad2xy(const fits::wcs& w, const vec<1,T>& ra, const vec<1,U>& dec,
        vec<1,V>& x, vec<1,W>& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(ra.size() == dec.size(), "RA and Dec arrays do not match sizes ("+
            strn(ra.size())+" vs "+strn(dec.size())+")");

        std::size_t ngal = ra.size();

        vec1d world(2*ngal);
        vec1u ids1 = 2*uindgen(ngal);
        vec1u ids2 = ids1+1;
        world.safe[ids1] = ra;
        world.safe[ids2] = dec;

        vec1d pos(2*ngal);

        std::vector<double> phi(ngal), theta(ngal);
        std::vector<double> itmp(2*ngal);
        std::vector<int>    stat(ngal);

        wcss2p(w.w, ngal, 2, world.data.data(), phi.data(), theta.data(),
            itmp.data(), pos.data.data(), stat.data());

        x = pos.safe[ids1];
        y = pos.safe[ids2];
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W>
    void xy2ad(const fits::wcs& w, const vec<1,T>& x, const vec<1,U>& y,
        vec<1,V>& ra, vec<1,W>& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else

        phypp_check(w.is_valid(), "invalid WCS data");
        phypp_check(x.size() == y.size(), "x and y arrays do not match sizes ("+
            strn(x.size())+" vs "+strn(y.size())+")");

        std::size_t ngal = x.size();

        vec1d map(2*ngal);
        vec1u ids1 = 2*uindgen(ngal);
        vec1u ids2 = ids1+1;
        map.safe[ids1] = x;
        map.safe[ids2] = y;

        vec1d world(2*ngal);

        std::vector<double> phi(ngal), theta(ngal);
        std::vector<double> itmp(2*ngal);
        std::vector<int>    stat(ngal);

        wcsp2s(w.w, ngal, 2, map.data.data(), itmp.data(), phi.data(), theta.data(),
            world.data.data(), stat.data());

        ra = world.safe[ids1];
        dec = world.safe[ids2];
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!is_vec<T>::value &&
            !is_vec<U>::value && !is_vec<V>::value && !is_vec<W>::value>::type>
    void ad2xy(const fits::wcs& w, const T& ra, const U& dec, V& x, W& y) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        phypp_check(w.is_valid(), "invalid WCS data");

        vec1d world(2);
        world.safe[0] = ra;
        world.safe[1] = dec;

        vec1d pos(2);

        double phi, theta;
        std::vector<double> itmp(2);
        int stat;

        wcss2p(w.w, 1, 2, world.data.data(), &phi, &theta, itmp.data(),
            pos.data.data(), &stat);

        x = pos.safe[0];
        y = pos.safe[1];
#endif
    }

    template<typename T = double, typename U = double, typename V, typename W,
        typename enable = typename std::enable_if<!is_vec<T>::value &&
            !is_vec<U>::value && !is_vec<V>::value && !is_vec<W>::value>::type>
    void xy2ad(const fits::wcs& w, const T& x, const U& y, V& ra, W& dec) {
#ifdef NO_WCSLIB
        static_assert(!std::is_same<T,T>::value, "WCS support is disabled, "
            "please enable the WCSLib library to use this function");
#else
        phypp_check(w.is_valid(), "invalid WCS data");

        vec1d map(2);
        map.safe[0] = x;
        map.safe[1] = y;

        vec1d world(2);

        double phi, theta;
        std::vector<double> itmp(2);
        int stat;

        wcsp2s(w.w, 1, 2, map.data.data(), itmp.data(), &phi, &theta,
            world.data.data(), &stat);

        ra = world.safe[0];
        dec = world.safe[1];
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
            if (!wcs.is_valid()) {
                warning("could not extract WCS information");
                note("parsing '", file, "'");
                return false;
            }

            // Convert radius to number of pixels
            vec1d r, d;
            fits::xy2ad(wcs, {0, 1}, {0, 0}, r, d);
            aspix = angdist(r.safe[0], d.safe[0], r.safe[1], d.safe[1]);

            return true;
        }
#endif
    }

    void header2fits_(fitsfile* fptr, const std::string& hdr) {
        int status = 0;
        std::size_t nentry = hdr.size()/80 + 1;
        fits_set_hdrsize(fptr, nentry, &status);

        for (uint_t i = 0; i < nentry; ++i) {
            std::string entry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
            if (start_with(entry, "END ")) continue;

            // Skip if it is an internal FITS keyword
            std::size_t eqpos = entry.find_first_of("=");
            if (eqpos != entry.npos) {
                std::string nam = trim(entry.substr(0, eqpos));
                if (nam == "SIMPLE" || nam == "BITPIX" || start_with(nam, "NAXIS") || nam == "EXTEND" ||
                    nam == "XTENSION" || nam == "EXTNAME" || nam == "PCOUNT" || nam == "GCOUNT") {
                    continue;
                }
            }

            fits_write_record(fptr, entry.c_str(), &status);
        }
    }

    // Write an image in a FITS file
    template<std::size_t Dim, typename Type>
    void write(const std::string& name, const vec<Dim,Type>& v, const fits::header& hdr) {
#ifdef NO_CCFITS
        static_assert(!std::is_same<Type,Type>::value, "CCfits is disabled, "
            "please enable the CCfits library to use this function");
#else
        std::array<long,Dim> isize;
        for (uint_t i = 0; i < Dim; ++i) {
            isize[i] = v.dims[Dim-1-i];
        }

        std::valarray<rtype_t<Type>> tv(v.size());
        for (uint_t i = 0; i < v.size(); ++i) {
            tv[i] = v.safe[i];
        }

        try {
            fits::FITS f("!"+name, traits<rtype_t<Type>>::image_type, Dim, isize.data());
            f.pHDU().write(1, tv.size(), tv);
            if (!hdr.empty()) {
                header2fits_(f.pHDU().fitsPointer(), hdr);
            }
            f.flush();
        } catch (fits::FitsException& e) {
            print("error: FITS: "+e.message());
            throw;
        }
#endif
    }

    template<std::size_t Dim, typename Type>
    void write(const std::string& name, const vec<Dim,Type>& v) {
        write(name, v, "");
    }


    struct column_info {
        std::string name;
        enum type_t {
            string, boolean, byte, integer, float_simple, float_double
        } type;
        vec1u dims;
        uint_t length = 1;
    };

    vec<1,column_info> read_table_columns(const std::string& filename) {
        vec<1,column_info> cols;

        fitsfile* fptr;
        int status = 0;

        fits::header hdr = read_header(filename);

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file");

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);

            for (uint_t c : range(ncol)) {
                column_info ci;

                char name[80];
                char type = 0;
                long repeat = 0;
                fits_get_bcolparms(fptr, c+1, name, nullptr, &type, &repeat,
                    nullptr, nullptr, nullptr, nullptr, &status);

                const uint_t max_dim = 20;
                long axes[max_dim];
                int naxis = 0;
                fits_read_tdim(fptr, c+1, max_dim, &naxis, axes, &status);

                if (!getkey(hdr, "TTYPE"+strn(c+1), ci.name)) {
                    continue;
                }

                ci.name = trim(ci.name);

                switch (type) {
                case 'B' : {
                    vec<1,char> v(repeat);
                    char def = traits<char>::def();
                    int null;
                    fits_read_col(fptr, traits<char>::ttype, c+1, 1, 1, repeat,
                        &def, v.data.data(), &null, &status);

                    if (min(v) == 0 && max(v) <= 1) {
                        ci.type = column_info::boolean;
                    } else {
                        ci.type = column_info::string;
                    }
                    break;
                }
                case 'S' : {
                    ci.type = column_info::byte;
                    break;
                }
                case 'J' :
                case 'I' : {
                    ci.type = column_info::integer;
                    break;
                }
                case 'E' : {
                    ci.type = column_info::float_simple;
                    break;
                }
                case 'D' : {
                    ci.type = column_info::float_double;
                    break;
                }
                default : {
                    warning("unhandled column type '", type, "' for ", ci.name);
                    continue;
                }
                }

                for (uint_t i : range(naxis)) {
                    // First dimension for string is the string length
                    if (i == 0 && ci.type == column_info::string) {
                        ci.length = axes[i];
                        if (naxis == 1) {
                            break;
                        }

                        continue;
                    }

                    ci.dims.push_back(axes[i]);
                }

                ci.dims = reverse(ci.dims);

                cols.push_back(ci);
            }

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }

        return cols;
    }

    std::string bake_macroed_name(const std::string& s) {
        return toupper(file::bake_macroed_name(s));
    }

    std::string type_to_string_(int type) {
        if (type == TSTRING) return "string";
        if (type == TSHORT) return "short";
        if (type == TLONG) return "long";
        if (type == TFLOAT) return "float";
        if (type == TDOUBLE) return "double";
        if (type == TLOGICAL) return "bool";
        if (type == TBIT) return "bool";
        if (type == TBYTE) return "char";
        if (type == TINT32BIT) return "int32";
        if (type == TCOMPLEX) return "complex float";
        if (type == TDBLCOMPLEX) return "complex double";
        return "unknown";
    }

    // Load the content of a FITS file into a set of arrays.
    template<std::size_t Dim, typename Type,
        typename enable = typename std::enable_if<!std::is_same<Type,std::string>::value>::type>
    void read_table_impl_(fitsfile* fptr, const std::string& tcolname,
        vec<Dim,Type>& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASEINSEN, const_cast<char*>(colname.c_str()), &id, &status);
        if (loose) {
            if (status != 0) return;
        } else {
            phypp_check_cfitsio(status, "cannot find collumn '"+colname+"' in FITS file");
        }

        int type, naxis;
        long repeat, width;
        std::array<long,Dim> naxes;
        fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
        fits_read_tdim(fptr, id, Dim, &naxis, naxes.data(), &status);
        phypp_check_fits(naxis == Dim, "wrong dimension for column '"+colname+"' "
            "(expected "+strn(Dim)+", got "+strn(naxis)+")");
        phypp_check_fits(traits<Type>::is_convertible(type), "wrong type for column '"+colname+"' "
            "(expected "+pretty_type_t(Type)+", got "+type_to_string_(type)+")");
        status = 0;

        for (uint_t i = 0; i < Dim; ++i) {
            v.dims[i] = naxes[Dim-1-i];
        }

        v.resize();

        Type def = traits<Type>::def();
        int null;
        fits_read_col(
            fptr, traits<Type>::ttype, id, 1, 1, repeat, &def, v.data.data(), &null, &status
        );
    }

    template<typename Type>
    void read_table_impl_(fitsfile* fptr, const std::string& tcolname, Type& v, bool loose = false) {
        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASEINSEN, const_cast<char*>(colname.c_str()), &id, &status);
        if (loose) {
            if (status != 0) return;
        } else {
            phypp_check_cfitsio(status, "cannot find collumn '"+colname+"' in FITS file");
        }

        int type, naxis;
        long repeat, width;
        std::array<long,1> naxes;
        fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
        fits_read_tdim(fptr, id, 1, &naxis, naxes.data(), &status);
        phypp_check_fits(naxis == 1, "wrong dimension for column '"+colname+"' "
            "(expected 1, got "+strn(naxis)+")");
        phypp_check_fits(traits<Type>::is_convertible(type), "wrong type for column '"+colname+"' "
            "(expected "+pretty_type_t(Type)+", got "+type_to_string_(type)+")");
        status = 0;

        Type def = traits<Type>::def();
        int null;
        fits_read_col(
            fptr, traits<Type>::ttype, id, 1, 1, repeat, &def,
            reinterpret_cast<typename traits<Type>::dtype*>(&v), &null, &status
        );
    }

    template<std::size_t Dim>
    void read_table_impl_(fitsfile* fptr, const std::string& tcolname,
        vec<Dim,std::string>& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASEINSEN, const_cast<char*>(colname.c_str()), &id, &status);
        if (loose) {
            if (status != 0) return;
        } else {
            phypp_check_cfitsio(status, "cannot find collumn '"+colname+"' in FITS file");
        }

        int type, naxis;
        long repeat, width;
        std::array<long,Dim+1> naxes;
        fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
        fits_read_tdim(fptr, id, Dim+1, &naxis, naxes.data(), &status);
        phypp_check_fits(naxis == Dim+1, "wrong dimension for column '"+colname+"' "
            "(expected "+strn(Dim)+", got "+strn(naxis-1)+")");
        phypp_check_fits(traits<std::string>::is_convertible(type), "wrong type for column '"+colname+"' "
            "(expected "+pretty_type_t(std::string)+", got "+type_to_string_(type)+")");
        status = 0;

        for (uint_t i = 0; i < Dim; ++i) {
            v.dims[i] = naxes[Dim-i];
        }

        v.resize();

        char* buffer = new char[n_elements(v)*naxes[0]];
        char def = '\0';
        int null;
        fits_read_col(
            fptr, traits<std::string>::ttype, id, 1, 1, n_elements(v)*naxes[0], &def,
            buffer, &null, &status
        );

        for (uint_t i = 0; i < n_elements(v); ++i) {
            v[i].reserve(naxes[0]);
            for (uint_t j = 0; j < uint_t(naxes[0]); ++j) {
                char c = buffer[i*naxes[0] + j];
                if (c == '\0') break;
                v[i].push_back(c);
            }
        }

        delete[] buffer;
    }

    void read_table_impl_(fitsfile* fptr, const std::string& tcolname,
        std::string& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASEINSEN, const_cast<char*>(colname.c_str()), &id, &status);
        if (loose) {
            if (status != 0) return;
        } else {
            phypp_check_cfitsio(status, "cannot find collumn '"+colname+"' in FITS file");
        }

        int type, naxis;
        long repeat, width;
        std::array<long,1> naxes;
        fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
        fits_read_tdim(fptr, id, 1, &naxis, naxes.data(), &status);
        phypp_check_fits(naxis == 1, "wrong dimension for column '"+colname+"' "
            "(expected 1, got "+strn(naxis)+")");
        phypp_check_fits(traits<std::string>::is_convertible(type), "wrong type for column '"+colname+"' "
            "(expected "+pretty_type_t(std::string)+", got "+type_to_string_(type)+")");
        status = 0;

        char* buffer = new char[naxes[0]+1];
        char def = '\0';
        int null;
        fits_read_col(
            fptr, traits<std::string>::ttype, id, 1, 1, naxes[0], &def,
            buffer, &null, &status
        );

        buffer[naxes[0]] = '\0';
        v = buffer;

        delete[] buffer;
    }

    template<typename T>
    void read_table_impl_(fitsfile* fptr, const std::string& colname,
        reflex::struct_t<T> data, bool loose = false);

    namespace impl {
        struct do_read_member {
            fitsfile* fptr;
            std::string base;
            bool loose;

            template<typename P>
            void operator () (reflex::member_t& m, P&& v) {
                read_table_impl_(fptr, base+toupper(m.name), std::forward<P>(v), loose);
            }
        };
    }

    template<typename T>
    void read_table_impl_(fitsfile* fptr, const std::string& colname,
        reflex::struct_t<T> data, bool loose) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        impl::do_read_member run{fptr, toupper(colname)+".", loose};
        reflex::foreach_member(data, run);
    }

    void read_table_(fitsfile* fptr, bool loose) {}

    template<typename T, typename ... Args>
    void read_table_(fitsfile* fptr, bool loose, const std::string& name, T& v, Args&& ... args) {
        read_table_impl_(fptr, name, v, loose);
        read_table_(fptr, loose, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void read_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file");

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);
            phypp_check_fits((sizeof...(Args)+1)/2 <= std::size_t(ncol),
                "too few collumns in this FITS file ("+
                strn(ncol)+" vs "+strn((sizeof...(Args)+1)/2)+")");

            read_table_(fptr, false, name, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename ... Args>
    void read_table_loose(const std::string& filename, const std::string& name, Args&& ... args) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file");

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);
            phypp_check_fits((sizeof...(Args)+1)/2 <= std::size_t(ncol),
                "too few collumns in this FITS file ("+
                strn(ncol)+" vs "+strn((sizeof...(Args)+1)/2)+")");

            read_table_(fptr, true, name, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    void read_table_(file::macroed_t, fitsfile* fptr, const std::string& names) {}

    template<typename T, typename ... Args>
    void read_table_(file::macroed_t, fitsfile* fptr, const std::string& names, T& v, Args&& ... args) {
        std::size_t pos = names.find_first_of(',');
        read_table_impl_(fptr, bake_macroed_name(names.substr(0, pos)), v);

        if (pos != names.npos) {
            read_table_(file::macroed_t(), fptr, names.substr(pos+1), std::forward<Args>(args)...);
        }
    }

    template<typename ... Args>
    void read_table(const std::string& filename, file::macroed_t,
        const std::string& names, Args&& ... args) {

        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file");

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);
            phypp_check_fits(sizeof...(Args) <= std::size_t(ncol), "too few collumns in this FITS "
                "file  ("+strn(ncol)+" vs "+strn(sizeof...(Args))+")");

            read_table_(file::macroed_t(), fptr, names, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    namespace impl {
        struct do_read_struct {
            fitsfile* fptr;
            bool loose;

            template<typename P>
            void operator () (reflex::member_t& m, P&& v) {
                read_table_impl_(fptr, toupper(m.name), std::forward<P>(v), loose);
            }
        };
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void read_table(const std::string& filename, T& t) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            reflex::foreach_member(reflex::wrap(t), impl::do_read_struct{fptr,false});

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void read_table_loose(const std::string& filename, T& t) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file '"+filename+"'");


            reflex::foreach_member(reflex::wrap(t), impl::do_read_struct{fptr,true});

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    // Write several columns in a FITS file.
    void write_table_(fitsfile* fptr, int id) {}

    template<std::size_t Dim, typename Type,
        typename enable = typename std::enable_if<!std::is_same<Type,std::string>::value>::type>
    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname,
        const vec<Dim,Type>& v) {

        std::string tform = strn(n_elements(v))+traits<Type>::tform;
        int status = 0;
        fits_insert_col(
            fptr, id, const_cast<char*>(colname.c_str()), const_cast<char*>(tform.c_str()), &status
        );

        std::array<long,Dim> dims;
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i] = v.dims[Dim-1-i];
        }

        fits_write_tdim(fptr, id, Dim, dims.data(), &status);

        fits_write_col(fptr, traits<Type>::ttype, id, 1, 1, n_elements(v),
            const_cast<typename vec<Dim,Type>::dtype*>(v.data.data()), &status);

        ++id;
    }

    template<typename Type>
    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname, const Type& v) {
        std::string tform = std::string(1, traits<Type>::tform);
        int status = 0;
        fits_insert_col(
            fptr, id, const_cast<char*>(colname.c_str()), const_cast<char*>(tform.c_str()), &status
        );

        using dtype = typename traits<Type>::dtype;
        fits_write_col(fptr, traits<Type>::ttype, id, 1, 1, 1,
            const_cast<dtype*>(reinterpret_cast<const dtype*>(&v)), &status);

        ++id;
    }

    template<std::size_t Dim>
    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname,
        const vec<Dim,std::string>& v) {

        std::size_t nmax = 0;
        for (auto& s : v) {
            if (s.size() > nmax) {
                nmax = s.size();
            }
        }

        std::string tform = strn(nmax*n_elements(v))+traits<std::string>::tform;

        int status = 0;
        fits_insert_col(
            fptr, id, const_cast<char*>(colname.c_str()), const_cast<char*>(tform.c_str()), &status
        );

        std::array<long,Dim+1> dims;
        dims[0] = nmax;
        for (uint_t i = 0; i < Dim; ++i) {
            dims[i+1] = v.dims[Dim-1-i];
        }

        fits_write_tdim(fptr, id, Dim+1, dims.data(), &status);

        char* buffer = new char[nmax*n_elements(v)];
        uint_t p = 0;
        for (auto& s : v) {
            for (auto c : s) {
                buffer[p] = c;
                ++p;
            }
            for (uint_t i = s.size(); i < nmax; ++i) {
                buffer[p] = '\0';
                ++p;
            }
        }

        fits_write_col(
            fptr, traits<std::string>::ttype, id, 1, 1, nmax*n_elements(v), buffer, &status
        );

        delete[] buffer;

        ++id;
    }

    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname,
        const std::string& v) {

        std::string tform = strn(v.size())+traits<std::string>::tform;
        int status = 0;
        fits_insert_col(
            fptr, id, const_cast<char*>(colname.c_str()), const_cast<char*>(tform.c_str()), &status
        );

        long dims = v.size();
        fits_write_tdim(fptr, id, 1, &dims, &status);

        fits_write_col(fptr, traits<std::string>::ttype, id, 1, 1, dims,
            const_cast<char*>(v.c_str()), &status);

        ++id;
    }

    template<typename T>
    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname,
        reflex::struct_t<T> data);

    namespace impl {
        struct do_write_member {
            fitsfile* fptr;
            int& id;
            std::string base;

            template<typename P>
            void operator () (const reflex::member_t& m, const P& v) {
                write_table_impl_(fptr, id, base+toupper(m.name), v);
            }
        };
    }

    template<typename T>
    void write_table_impl_(fitsfile* fptr, int& id, const std::string& colname,
        reflex::struct_t<T> data) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        impl::do_write_member run{fptr, id, colname+"."};
        reflex::foreach_member(data, run);
    }

    template<typename T, typename ... Args>
    void write_table_(fitsfile* fptr, int id, const std::string& name, const T& v, Args&& ... args) {
        write_table_impl_(fptr, id, toupper(name), v);
        write_table_(fptr, id, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void write_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_create_file(&fptr, ("!"+filename).c_str(), &status);
            phypp_check_cfitsio(status, "cannot open file");

            fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

            write_table_(fptr, 1, name, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    void write_table_(file::macroed_t, fitsfile* fptr, int id, const std::string& names) {}


    template<typename T, typename ... Args>
    void write_table_(file::macroed_t, fitsfile* fptr, int id, const std::string& names,
        const named_t<T>& v, Args&& ... args);

    template<typename T, typename ... Args>
    void write_table_(file::macroed_t, fitsfile* fptr, int id, const std::string& names,
        const T& v, Args&& ... args) {

        std::size_t pos = names.find_first_of(',');
        write_table_impl_(fptr, id, bake_macroed_name(names.substr(0, pos)), v);

        if (pos != names.npos) {
            write_table_(file::macroed_t(), fptr, id, names.substr(pos+1), std::forward<Args>(args)...);
        }
    }

    template<typename T, typename ... Args>
    void write_table_(file::macroed_t, fitsfile* fptr, int id, const std::string& names,
        const named_t<T>& v, Args&& ... args) {

        write_table_impl_(fptr, id, v.name, v.obj);

        std::size_t pos = names.find_first_of(')');
        pos = names.find_first_of(',', pos);

        if (pos != names.npos) {
            write_table_(file::macroed_t(), fptr, id, names.substr(pos+1), std::forward<Args>(args)...);
        }
    }

    template<typename ... Args>
    void write_table(const std::string& filename, file::macroed_t,
        const std::string& names, Args&& ... args) {

        fitsfile* fptr;
        int status = 0;

        try {
            fits_create_file(&fptr, ("!"+filename).c_str(), &status);
            phypp_check_cfitsio(status, "cannot open file");

            fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

            write_table_(file::macroed_t(), fptr, 1, names, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    namespace impl {
        struct do_write_struct {
            do_write_struct(fitsfile* p, int i = 1) : fptr(p), id(i) {}

            fitsfile* fptr;
            int id = 1;

            template<typename P>
            void operator () (const reflex::member_t& m, const P& v) {
                write_table_impl_(fptr, id, toupper(m.name), v);
            }
        };
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void write_table(const std::string& filename, const T& t) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        fitsfile* fptr;
        int status = 0;

        try {
            fits_create_file(&fptr, ("!"+filename).c_str(), &status);
            phypp_check_cfitsio(status, "cannot open file");
            fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

            reflex::foreach_member(reflex::wrap(t), impl::do_write_struct(fptr, 1));

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    // Append an new column to an existing FITS table
    // Note: if the column already exsits in the file, it will be overwritten

    void remove_column_(fitsfile* fptr, const std::string& name) {
        int id = 0, status = 0;
        fits_get_colnum(fptr, CASEINSEN, const_cast<char*>(name.c_str()), &id, &status);
        if (status == 0) {
            fits_delete_col(fptr, id, &status);
        }
    }

    void remove_columns_(fitsfile* fptr) {}

    template<typename T, typename ... Args>
    void remove_columns_(fitsfile* fptr, const std::string& name, const T& t, const Args& ... args) {
        remove_column_(fptr, name);
        remove_columns_(fptr, args...);
    }

    template<typename ... Args>
    void remove_columns_(file::macroed_t, fitsfile* fptr, const std::string& names, const Args& ... args) {
        vec1s vnames = split(names.substr(names.find_first_of(')')+1), ",");
        for (auto& s : vnames) {
            s = bake_macroed_name(s);
            if (s.empty()) continue;
            remove_column_(fptr, s);
        }
    }

    template<typename ... Args>
    void update_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READWRITE, &status);
            phypp_check_cfitsio(status, "cannot open file");

            remove_columns_(fptr, name, args...);

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);
            write_table_(fptr, ncol+1, name, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename ... Args>
    void update_table(const std::string& filename, file::macroed_t,
        const std::string& names, Args&& ... args) {

        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READWRITE, &status);
            phypp_check_cfitsio(status, "cannot open file");

            remove_columns_(file::macroed_t(), fptr, names, std::forward<Args>(args)...);

            int ncol;
            fits_get_num_cols(fptr, &ncol, &status);
            write_table_(file::macroed_t(), fptr, ncol+1, names, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    void display(const std::string& name) {
        fork("ds9 "+name);
    }

    void display(const std::string& name1, const std::string& name2) {
        fork("ds9 -rgb -red "+name1+" -green "+name2);
    }

    void display(const std::string& name1, const std::string& name2, const std::string& name3) {
        fork("ds9 -rgb -red "+name1+" -green "+name2+" -blue "+name3);
    }
}

#endif

