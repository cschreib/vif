#ifndef FITS_HPP
#define FITS_HPP

#include <CCfits/CCfits>
#include <wcslib/wcshdr.h>
#include <string>
#include "vec.hpp"
#include "reflex.hpp"
#include "file.hpp"
#include "print.hpp"
#include "math.hpp"

namespace fits {
    using namespace CCfits;

    struct exception {
        explicit exception(const std::string& m) : msg(m) {}
        std::string msg;
    };

    template<typename T>
    struct traits;

    template<>
    struct traits<std::string> {
        using dtype = char;

        static const char tform = 'B';
        static const int ttype = TBYTE;
        static const fits::ValueType type = (fits::ValueType)TBYTE;

        static std::string def() {
            return "";
        }
    };

    template<>
    struct traits<bool> {
        using dtype = char;

        static const char tform = 'B';
        static const int ttype = TBYTE;
        static const fits::ValueType type = (fits::ValueType)TBYTE;
        static const int image_type = BYTE_IMG;

        static bool def() {
            return false;
        }
    };

    template<>
    struct traits<char> {
        using dtype = char;

        static const char tform = 'S';
        static const int ttype = TSBYTE;
        static const fits::ValueType type = (fits::ValueType)TSBYTE;
        static const int image_type = SBYTE_IMG;

        static bool def() {
            return '\0';
        }
    };

    template<>
    struct traits<uint_t> {
        using dtype = uint_t;

        static const char tform = 'J';
        static const int ttype = TLONG;
        static const fits::ValueType type = fits::Tlong;
        static const int image_type = LONG_IMG;
        static uint_t def() {
            return 0;
        }
    };

    template<>
    struct traits<int_t> {
        using dtype = int_t;

        static const char tform = 'J';
        static const int ttype = TLONG;
        static const fits::ValueType type = fits::Tlong;
        static const int image_type = LONG_IMG;
        static int_t def() {
            return 0;
        }
    };

    template<>
    struct traits<float> {
        using dtype = float;

        static const char tform = 'E';
        static const int ttype = TFLOAT;
        static const fits::ValueType type = fits::Tfloat;
        static const int image_type = FLOAT_IMG;
        static float def() {
            return fnan;
        }
    };

    template<>
    struct traits<double> {
        using dtype = double;

        static const char tform = 'D';
        static const int ttype = TDOUBLE;
        static const fits::ValueType type = fits::Tdouble;
        static const int image_type = DOUBLE_IMG;
        static float def() {
            return dnan;
        }
    };

    #define phypp_check_fits(assertion, msg) \
        do { if (!(assertion)) throw fits::exception(msg); } while(0)

    void phypp_check_cfitsio(int status, const std::string& msg) {
        if (status != 0) {
            char txt[FLEN_STATUS];
            fits_get_errstatus(status, txt);
            phypp_check_fits(status == 0, "error: cfitsio: "+std::string(txt)+"\n"+msg);
        }
    }

    using header = std::string;

    template<std::size_t Dim, typename Type>
    void read_impl_(fits::FITS& file, vec_t<Dim, Type>& v, std::string& hdr, bool gethdr) {
        fits::PHDU& phdu = file.pHDU();
        if (phdu.axes() == 0) {
            fits::ExtHDU* hdu = nullptr;
            uint_t i = 1;
            while (!hdu && i <= file.extension().size()) {
                hdu = &file.extension(i);
                if (hdu->axes() == 0) hdu = nullptr;
            }

            phypp_check_fits(hdu, "FITS file does not contain any data");

            fits::ExtHDU& thdu = *hdu;

            phypp_check_fits(thdu.axes() == Dim, "FITS file does not match array dimensions ("+
                strn(thdu.axes())+" vs "+strn(Dim)+")");

            std::size_t n = 1;
            for (uint_t j = 0; j < Dim; ++j) {
                v.dims[j] = thdu.axis(Dim-1-j);
                n *= v.dims[j];
            }

            std::valarray<Type> tv(n);
            thdu.read(tv, 1, n);

            v.data.assign(std::begin(tv), std::end(tv));

            if (gethdr) {
                // Read the header as a string
                char* hstr = nullptr;
                int nkeys  = 0;
                int status = 0;
                fits_hdr2str(thdu.fitsPointer(), 0, nullptr, 0, &hstr, &nkeys, &status);
                hdr = hstr;
                free(hstr);
            }
        } else {
            phypp_check_fits(phdu.axes() == Dim, "FITS file does not match array dimensions ("+
                strn(phdu.axes())+" vs "+strn(Dim)+")");

            std::size_t n = 1;
            for (uint_t i = 0; i < Dim; ++i) {
                v.dims[i] = phdu.axis(Dim-1-i);
                n *= v.dims[i];
            }

            std::valarray<Type> tv(n);
            phdu.read(tv, 1, n);

            v.data.assign(std::begin(tv), std::end(tv));

            if (gethdr) {
                // Read the header as a string
                char* hstr = nullptr;
                int nkeys  = 0;
                int status = 0;
                fits_hdr2str(phdu.fitsPointer(), 0, nullptr, 0, &hstr, &nkeys, &status);
                hdr = hstr;
                free(hstr);
            }
        }
    }

    // Load the content of a FITS file into an array.
    // Be careful that C++ is a row-major language, in the sense that to get the (x,y) pixel of a
    // given image, one has to get the (y,x) element in the array.
    template<std::size_t Dim, typename Type>
    vec_t<Dim, Type> read(const std::string& name) {
        try {
            fits::FITS file(name, fits::Read);
            vec_t<Dim, Type> v;
            std::string hdr;
            read_impl_(file, v, hdr, false);
            return v;
        } catch (fits::FitsException& e) {
            error("reading: "+name);
            error("FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
    }

    template<std::size_t Dim, typename Type>
    void read(const std::string& name, vec_t<Dim, Type>& v) {
        try {
            fits::FITS file(name, fits::Read);
            std::string hdr;
            read_impl_(file, v, hdr, false);
        } catch (fits::FitsException& e) {
            print("error: FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
    }

    template<std::size_t Dim, typename Type>
    void read(const std::string& name, vec_t<Dim, Type>& v, fits::header& hdr) {
        try {
            fits::FITS file(name, fits::Read);
            read_impl_(file, v, hdr, true);
        } catch (fits::FitsException& e) {
            print("error: FITS: "+e.message());
            throw;
        } catch (fits::exception& e) {
            error("reading: "+name);
            error(e.msg);
            throw;
        }
    }

    fits::header read_header(const std::string& filename) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_image(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

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
        std::string entry = key+std::string(8-key.size(), ' ')+"= ";
        std::string value = strn(v);
        if (!comment.empty()) {
            value += " / "+comment;
        }
        if (value.size() > 80) return false;
        entry += value;
        entry += std::string(80-entry.size(), ' ');

        std::size_t pos = hdr.find_last_of('E');
        hdr.insert(pos, entry);
        return true;
    }

    // Extract astrometry from a FITS image header
    struct wcs {
        wcsprm* w = nullptr;
        int nwcs  = 0;

        wcs() = default;

        explicit wcs(const fits::header& hdr) {
            // Feed it to WCSLib to extract the astrometric parameters
            int nreject = 0;
            wcspih(
                const_cast<char*>(hdr.c_str()), hdr.size()/80 + 1,
                WCSHDR_all, 0, &nreject, &nwcs, &w
            );
        }

        wcs(const wcs&) = delete;
        wcs& operator = (const wcs&) = delete;
        wcs(wcs&& tw) {
            std::swap(w, tw.w);
            std::swap(nwcs, tw.nwcs);
        }
        wcs& operator = (wcs&& tw) {
            w = tw.w; tw.w = nullptr;
            nwcs = tw.nwcs; nwcs = 0;
            return *this;
        }

        ~wcs() {
            if (w) {
                wcsvfree(&nwcs, &w);
            }
        }
    };

    fits::wcs extast(const fits::header& hdr) {
        return fits::wcs(hdr);
    }

    template<typename T, typename U, typename V, typename W>
    void ad2xy(const fits::wcs& w, const vec_t<1,T>& ra, const vec_t<1,U>& dec,
        vec_t<1,V>& x, vec_t<1,W>& y) {

        phypp_check(w.w != nullptr, "uninitialized WCS structure");
        phypp_check(ra.size() == dec.size(), "RA and Dec arrays do not match sizes ("+
            strn(ra.size())+" vs "+strn(dec.size())+")");

        std::size_t ngal = ra.size();

        vec1d world = dblarr(2*ngal);
        vec1u ids1 = 2*uindgen(ngal);
        vec1u ids2 = ids1+1;
        world[ids1] = ra;
        world[ids2] = dec;

        vec1d pos = dblarr(2*ngal);

        std::vector<double> phi(ngal), theta(ngal);
        std::vector<double> itmp(2*ngal);
        std::vector<int>    stat(ngal);

        wcss2p(w.w, ngal, 2, world.data.data(), phi.data(), theta.data(),
            itmp.data(), pos.data.data(), stat.data());

        x = pos[ids1];
        y = pos[ids2];
    }

    // Write an image in a FITS file
    template<std::size_t Dim, typename Type>
    void write(const std::string& name, const vec_t<Dim,Type>& v) {
        std::array<long,Dim> isize;
        for (uint_t i = 0; i < Dim; ++i) {
            isize[i] = v.dims[Dim-1-i];
        }

        std::valarray<rtype_t<Type>> tv(v.size());
        for (uint_t i = 0; i < v.size(); ++i) {
            tv[i] = dref(v.data[i]);
        }

        try {
            fits::FITS f("!"+name, traits<rtype_t<Type>>::image_type, Dim, isize.data());
            f.pHDU().write(1, tv.size(), tv);
            f.flush();
        } catch (fits::FitsException& e) {
            print("error: FITS: "+e.message());
            throw;
        }
    }

    void header2fits_(fits::HDU& hdu, const std::string& hdr) {
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i = 0; i < nentry; ++i) {
            std::string entry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
            std::size_t eqpos = entry.find_first_of("=");
            if (eqpos == entry.npos) continue;
            std::size_t cpos = entry.find_first_of("/", eqpos);

            std::string nam = trim(entry.substr(0, eqpos));
            if (nam == "SIMPLE" || nam == "BITPIX" || nam.find("NAXIS") == 0 || nam == "EXTEND") {
                continue;
            }

            std::string value = trim(entry.substr(eqpos+1, cpos-eqpos-1));
            std::string comment = trim(entry.substr(cpos+1));

            if (value.find_first_of("'") == 0) {
                value = trim(value, "'");
                hdu.addKey(nam, value, comment);
            } else if (value.find(".") != value.npos) {
                double dvalue; from_string(value, dvalue);
                hdu.addKey(nam, dvalue, comment);
            } else {
                std::ptrdiff_t ivalue; from_string(value, ivalue);
                hdu.addKey(nam, ivalue, comment);
            }
        }
    }

    template<std::size_t Dim, typename Type>
    void write(const std::string& name, const vec_t<Dim,Type>& v, const fits::header& hdr) {
        std::array<long,Dim> isize;
        for (uint_t i = 0; i < Dim; ++i) {
            isize[i] = v.dims[Dim-1-i];
        }

        std::valarray<rtype_t<Type>> tv(v.size());
        for (uint_t i = 0; i < v.size(); ++i) {
            tv[i] = dref(v.data[i]);
        }

        try {
            fits::FITS f("!"+name, traits<rtype_t<Type>>::image_type, Dim, isize.data());
            f.pHDU().write(1, tv.size(), tv);
            header2fits_(f.pHDU(), hdr);
            f.flush();
        } catch (fits::FitsException& e) {
            print("error: FITS: "+e.message());
            throw;
        }
    }

    struct macroed_t {};
    #define ftable(...) fits::macroed_t(), #__VA_ARGS__, __VA_ARGS__

    std::string bake_macroed_name(const std::string& s) {
        std::string t = toupper(trim(s));
        auto p = t.find_first_of('.');
        if (p != t.npos) {
            t = t.substr(p+1);
        }
        return t;
    }

    // Load the content of a FITS file into a set of arrays.
    void read_table_(fitsfile* fptr) {}

    template<std::size_t Dim, typename Type,
        typename enable = typename std::enable_if<!std::is_same<Type,std::string>::value>::type>
    void read_table_impl_(fitsfile* fptr, const std::string& tcolname,
        vec_t<Dim,Type>& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASESEN, const_cast<char*>(colname.c_str()), &id, &status);
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
        phypp_check_fits(naxis == Dim, "wrong dimension for column '"+colname+"'");
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
        fits_get_colnum(fptr, CASESEN, const_cast<char*>(colname.c_str()), &id, &status);
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
        phypp_check_fits(naxis == 1, "wrong dimension for column '"+colname+"'");
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
        vec_t<Dim,std::string>& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASESEN, const_cast<char*>(colname.c_str()), &id, &status);
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
        phypp_check_fits(naxis == Dim+1, "wrong dimension for column '"+colname+"'");
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

        delete buffer;
    }

    void read_table_impl_(fitsfile* fptr, const std::string& tcolname,
        std::string& v, bool loose = false) {

        int id, status = 0;
        std::string colname = toupper(tcolname);
        fits_get_colnum(fptr, CASESEN, const_cast<char*>(colname.c_str()), &id, &status);
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
        phypp_check_fits(naxis == 1, "wrong dimension for column '"+colname+"'");
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

        delete buffer;
    }

    template<typename T>
    void read_table_impl_(fitsfile* fptr, const std::string& colname,
        reflex::struct_t<T> data, bool loose = false) {

        struct {
            fitsfile* fptr;
            std::string base;
            bool loose;

            template<typename P>
            void operator () (reflex::member_t& m, P&& v) {
                read_table_impl_(this->fptr, this->base+toupper(m.name), std::forward<P>(v), this->loose);
            }
        } do_read{fptr, toupper(colname)+".", loose};

        reflex::foreach_member(data, do_read);
    }

    template<typename T, typename ... Args>
    void read_table_(fitsfile* fptr, const std::string& name, T& v, Args&& ... args) {
        read_table_impl_(fptr, name, v);
        read_table_(fptr, std::forward<Args>(args)...);
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

            read_table_(fptr, name, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+name);
            error(e.msg);
            throw;
        }
    }

    void read_table_(macroed_t, fitsfile* fptr, const std::string& names) {}

    template<typename T, typename ... Args>
    void read_table_(macroed_t, fitsfile* fptr, const std::string& names, T& v, Args&& ... args) {
        std::size_t pos = names.find_first_of(',');
        read_table_impl_(fptr, bake_macroed_name(names.substr(0, pos)), v);

        if (pos != names.npos) {
            read_table_(macroed_t(), fptr, names.substr(pos+1), std::forward<Args>(args)...);
        }
    }

    template<typename ... Args>
    void read_table(const std::string& filename, macroed_t,
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

            read_table_(macroed_t(), fptr, names, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("reading: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void read_table(const std::string& filename, T& t) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            struct {
                fitsfile* fptr;

                template<typename P>
                void operator () (reflex::member_t& m, P&& v) {
                    read_table_impl_(this->fptr, toupper(m.name), std::forward<P>(v), true);
                }
            } do_read{fptr};

            reflex::foreach_member(reflex::wrap(t), do_read);

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
        fitsfile* fptr;
        int status = 0;

        try {
            fits_open_table(&fptr, filename.c_str(), READONLY, &status);
            phypp_check_cfitsio(status, "cannot open file '"+filename+"'");

            struct {
                fitsfile* fptr;

                template<typename P>
                void operator () (reflex::member_t& m, P&& v) {
                    read_table_impl_(this->fptr, toupper(m.name), std::forward<P>(v), true);
                }
            } do_read{fptr};

            reflex::foreach_member(reflex::wrap(t), do_read);

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
        const vec_t<Dim,Type>& v) {

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
            const_cast<typename vec_t<Dim,Type>::dtype*>(v.data.data()), &status);

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
        const vec_t<Dim,std::string>& v) {

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

        delete buffer;

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
        reflex::struct_t<T> data) {

        struct {
            fitsfile* fptr;
            int id;
            std::string base;

            template<typename P>
            void operator () (const reflex::member_t& m, const P& v) {
                write_table_impl_(this->fptr, this->id, this->base+toupper(m.name), v);
            }
        } do_write{fptr, id, colname+"."};

        reflex::foreach_member(data, do_write);
    }

    template<typename T, typename ... Args>
    void write_table_(fitsfile* fptr, int id, const std::string& name, const T& v, Args&& ... args) {
        write_table_impl_(fptr, id, name, v);
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

    void write_table_(macroed_t, fitsfile* fptr, int id, const std::string& names) {}

    template<typename T, typename ... Args>
    void write_table_(macroed_t, fitsfile* fptr, int id, const std::string& names,
        const T& v, Args&& ... args) {

        std::size_t pos = names.find_first_of(',');
        write_table_impl_(fptr, id, bake_macroed_name(names.substr(0, pos)), v);

        if (pos != names.npos) {
            write_table_(macroed_t(), fptr, id, names.substr(pos+1), std::forward<Args>(args)...);
        }
    }

    template<typename ... Args>
    void write_table(const std::string& filename, macroed_t,
        const std::string& names, Args&& ... args) {

        fitsfile* fptr;
        int status = 0;

        try {
            fits_create_file(&fptr, ("!"+filename).c_str(), &status);
            phypp_check_cfitsio(status, "cannot open file");

            fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

            write_table_(macroed_t(), fptr, 1, names, std::forward<Args>(args)...);

            fits_close_file(fptr, &status);
        } catch (fits::exception& e) {
            fits_close_file(fptr, &status);
            error("writing: "+filename);
            error(e.msg);
            throw;
        }
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void write_table(const std::string& filename, const T& t) {
        fitsfile* fptr;
        int status = 0;

        try {
            fits_create_file(&fptr, ("!"+filename).c_str(), &status);
            phypp_check_cfitsio(status, "cannot open file");
            fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

            struct {
                fitsfile* fptr;
                int id = 1;

                template<typename P>
                void operator () (const reflex::member_t& m, const P& v) {
                    write_table_impl_(this->fptr, this->id, toupper(m.name), v);
                }
            } do_write{fptr};

            reflex::foreach_member(reflex::wrap(t), do_write);

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

