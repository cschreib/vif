#ifndef VIF_IO_FITS_BASE_HPP
#define VIF_IO_FITS_BASE_HPP

#include <string>
#include "vif/core/vec.hpp"
#include "vif/core/print.hpp"
#include "vif/core/error.hpp"
#include "vif/utility/generic.hpp"
#include "vif/io/filesystem.hpp"
#include "vif/io/ascii.hpp"

#ifndef NO_CFITSIO

#include <fitsio.h>

namespace vif {
namespace fits {
    struct exception : std::exception {
        explicit exception(const std::string& m) : msg(m) {}
        std::string msg;

        const char* what() const noexcept override {
            return msg.c_str();
        }
    };

    #define vif_check_fits(assertion, msg) \
        do { if (!(assertion)) throw vif::fits::exception(msg); } while(0)

    inline void vif_check_cfitsio(int status, const std::string& msg) {
        if (status != 0) {
            char txt[FLEN_STATUS];
            fits_get_errstatus(status, txt);
            char tdetails[FLEN_ERRMSG];
            std::string details;
            while (fits_read_errmsg(tdetails)) {
                if (!details.empty()) details += "\ncfitsio: ";
                details += std::string(tdetails);
            }

            vif_check_fits(status == 0, "error: cfitsio: "+std::string(txt)+"\ncfitsio: "+
                details+"\nerror: "+msg);
        }
    }

    using header = std::string;

    enum hdu_type {
        null_hdu = 0,
        empty_hdu,
        image_hdu,
        table_hdu
    };

    // Format of FITS table (row/column oriented)
    enum class table_format {
        column_oriented,
        row_oriented
    };
}

namespace impl {
    namespace fits_impl {
        template<typename T>
        struct traits;

        template<>
        struct traits<std::string> {
            using dtype = char;

            static const char tform = 'A';
            static const int ttype = TSTRING;

            static std::string def() {
                return "";
            }

            static bool is_convertible(int type) {
                if (type == TSTRING) return true;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSTRING) return true;
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

            static bool is_convertible_narrow(int type) {
                if (type == TLOGICAL) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TINT32BIT) return true;
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

            static bool is_convertible_narrow(int type) {
                if (type == TBYTE) return true;
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TINT32BIT) return true;
                return false;
            }
        };

        template<>
        struct traits<uint_t> {
            using dtype = uint_t;

            static const char tform = 'K';
            static const int ttype = TLONGLONG;
            static const int image_type = LONGLONG_IMG;

            static uint_t def() {
                return 0;
            }

            static bool is_convertible(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }
        };

        template<>
        struct traits<int_t> {
            using dtype = int_t;

            static const char tform = 'K';
            static const int ttype = TLONGLONG;
            static const int image_type = LONGLONG_IMG;

            static int_t def() {
                return 0;
            }

            static bool is_convertible(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }
        };

        template<>
        struct traits<short> {
            using dtype = short;

            static const char tform = 'I';
            static const int ttype = TSHORT;
            static const int image_type = SHORT_IMG;

            static short def() {
                return 0;
            }

            static bool is_convertible(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return false;
                if (type == TLONGLONG) return false;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return false;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
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
                if (type == TLONGLONG) return true;
                if (type == TFLOAT) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TFLOAT) return true;
                if (type == TDOUBLE) return true;
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
                if (type == TLONGLONG) return true;
                if (type == TFLOAT) return true;
                if (type == TDOUBLE) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }

            static bool is_convertible_narrow(int type) {
                if (type == TSHORT) return true;
                if (type == TLONG) return true;
                if (type == TLONGLONG) return true;
                if (type == TFLOAT) return true;
                if (type == TDOUBLE) return true;
                if (type == TBIT) return true;
                if (type == TBYTE) return true;
                if (type == TINT32BIT) return true;
                return false;
            }
        };

        inline int bitpix_to_type(int bitpix) {
            switch(bitpix) {
                case BYTE_IMG     : return TBYTE;
                case SHORT_IMG    : return TSHORT;
                case LONG_IMG     : return TLONG;
                case LONGLONG_IMG : return TLONGLONG;
                case FLOAT_IMG    : return TFLOAT;
                case DOUBLE_IMG   : return TDOUBLE;
                default : throw fits::exception("unknown image type '"+to_string(bitpix)+"'");
            }
        }

        inline std::string type_to_string_(int type) {
            if (type == TSTRING) return "string";
            if (type == TSHORT) return "short";
            if (type == TLONG) return "long";
            if (type == TLONGLONG) return "long long";
            if (type == TFLOAT) return "float";
            if (type == TDOUBLE) return "double";
            if (type == TLOGICAL) return "bool";
            if (type == TBIT) return "bool";
            if (type == TBYTE) return "char";
            if (type == TINT32BIT) return "int32";
            if (type == TCOMPLEX) return "complex float";
            if (type == TDBLCOMPLEX) return "complex double";
            return "unknown ("+to_string(type)+")";
        }

        enum file_type {
            generic_file,
            image_file,
            table_file
        };

        enum access_right  {
            write_only,
            read_only,
            read_write
        };

        template<typename T>
        struct is_string : std::false_type {};

        template<>
        struct is_string<std::string> : std::true_type {};
        template<>
        struct is_string<char*> : std::true_type {};
        template<>
        struct is_string<const char*> : std::true_type {};
        template<std::size_t N>
        struct is_string<char[N]> : std::true_type {};
        template<std::size_t N>
        struct is_string<const char[N]> : std::true_type {};

        class file_base {
        public :
            file_base(file_type type, access_right rights) :
                type_(type), rights_(rights) {}

            file_base(file_type type, const std::string& filename, access_right rights) :
                type_(type), rights_(rights), filename_(filename) {

                open(filename);
            }

            file_base(file_base&& in) noexcept : type_(in.type_), rights_(in.rights_),
                filename_(in.filename_), fptr_(in.fptr_), status_(in.status_) {
                in.fptr_ = nullptr;
            }

            virtual ~file_base() noexcept {
                close();
            }

            void open(const std::string& filename) {
                close();

                if (rights_ == write_only || (rights_ == read_write && !file::exists(filename))) {
                    fits_create_file(&fptr_, ("!"+filename).c_str(), &status_);
                    fits::vif_check_cfitsio(status_, "cannot create file '"+filename+"'");
                } else {
                    int trights = (rights_ == read_only ? READONLY : READWRITE);

                    switch (type_) {
                        case generic_file :
                            fits_open_file(&fptr_, filename.c_str(), trights, &status_);
                            break;
                        case image_file :
                            fits_open_image(&fptr_, filename.c_str(), trights, &status_);
                            break;
                        case table_file :
                            fits_open_table(&fptr_, filename.c_str(), trights, &status_);
                            break;
                    }

                    fits::vif_check_cfitsio(status_, "cannot open file '"+filename+"'");
                }

                filename_ = filename;

                update_internal_state();
            }

            void close() noexcept {
                if (!fptr_) return;

                status_ = 0;
                fits_close_file(fptr_, &status_);
                fptr_ = nullptr;
            }

            bool is_open() const noexcept {
                return fptr_ != nullptr;
            }

            const std::string& filename() const noexcept {
                return filename_;
            }

            int cfitsio_status() const noexcept {
                return status_;
            }

            fitsfile* cfitsio_ptr() noexcept {
                return fptr_;
            }

            const fitsfile* cfitsio_ptr() const noexcept {
                return fptr_;
            }

            virtual void update_internal_state() {
                get_table_format_();
            }

        protected :

            void check_is_open_() const {
                if (!is_open()) {
                    throw fits::exception("no FITS file opened, cannot proceed");
                }
            }

            void get_table_format_() {
                // Check if data exists
                uint_t nhdu = hdu_count();
                if (nhdu != 0) {
                    int ncols = 0;
                    fits_get_num_cols(fptr_, &ncols, &status_);
                    if (status_ == 0 && ncols != 0) {
                        // Data exists, see if row or column-oriented
                        uint_t nrow;
                        if (read_keyword("NAXIS2", nrow)) {
                            if (nrow > 1 || nrow == 0) {
                                format_ = fits::table_format::row_oriented;
                                return;
                            }
                        }
                    }

                    status_ = 0;
                }

                // Default
                format_ = fits::table_format::column_oriented;
            }

        public :

            fits::header read_header() const {
                check_is_open_();

                status_ = 0;
                fits::header hdr;
                char* hstr = nullptr;
                int nkeys  = 0;
                fits_hdr2str(fptr_, 0, nullptr, 0, &hstr, &nkeys, &status_);
                fits::vif_check_cfitsio(status_, "could not dump header to string");
                hdr = hstr;
                fits_free_memory(hstr, &status_);
                return hdr;
            }

            bool has_keyword(const std::string& name) const {
                check_is_open_();

                char comment[80];
                char content[80];
                int status = 0;
                fits_read_keyword(fptr_, name.c_str(), content, comment, &status);
                return status == 0;
            }

            bool read_keyword(const std::string& name, std::string& value) const {
                check_is_open_();

                char comment[80];
                char content[80];
                int status = 0;
                fits_read_keyword(fptr_, name.c_str(), content, comment, &status);
                value = content;
                value = trim(value, " '");
                return status == 0;
            }

            template <typename T>
            bool read_keyword(const std::string& name, T& value) const {
                check_is_open_();

                std::string content;
                if (!read_keyword(name, content)) return false;

                // Convert FITS boolean values to something we can read
                if (content == "T") {
                    content = "1";
                } else if (content == "F") {
                    content = "0";
                }

                return from_string(content, value);
            }

            uint_t hdu_count() const {
                check_is_open_();

                int nhdu = 0;
                fits_get_num_hdus(fptr_, &nhdu, &status_);
                fits::vif_check_cfitsio(status_, "could not get the number of HDUs");
                return nhdu;
            }

            uint_t current_hdu() const {
                check_is_open_();

                int hdu = 1;
                fits_get_hdu_num(fptr_, &hdu);
                return hdu-1;
            }

            fits::hdu_type hdu_type() const {
                check_is_open_();

                bool simple = false;
                if (read_keyword("SIMPLE", simple)) {
                    if (axis_count() == 0) {
                        return fits::empty_hdu;
                    } else {
                        return fits::image_hdu;
                    }
                } else {
                    std::string xtension;
                    if (read_keyword("XTENSION", xtension)) {
                        if (xtension == "IMAGE") {
                            if (axis_count() == 0) {
                                return fits::empty_hdu;
                            } else {
                                return fits::image_hdu;
                            }
                        } else if (xtension == "BINTABLE" || xtension == "TABLE") {
                            return fits::table_hdu;
                        } else {
                            vif_check(false, "unknown XTENSION value '", xtension, "'");
                        }
                    }
                }

                return fits::null_hdu;
            }

            void reach_hdu(uint_t hdu) {
                check_is_open_();

                uint_t nhdu = hdu_count();
                if (rights_ == read_only) {
                    vif_check(hdu < nhdu, "requested HDU does not exists in this "
                        "FITS file (", hdu, " vs. ", nhdu, ")");
                } else {
                    if (hdu >= nhdu) {
                        // Create missing HDUs to be able to reach the one requested
                        long naxes = 0;
                        for (uint_t i = nhdu; i <= hdu; ++i) {
                            fits_insert_img(fptr_, impl::fits_impl::traits<float>::image_type,
                                0, &naxes, &status_);
                            fits::vif_check_cfitsio(status_,
                                "could not create new HDUs to reach HDU "+to_string(hdu));
                        }

                        fits_set_hdustruc(fptr_, &status_);
                        fits::vif_check_cfitsio(status_,
                            "could not create new HDUs to reach HDU "+to_string(hdu));
                    }
                }

                fits_movabs_hdu(fptr_, hdu+1, nullptr, &status_);
                fits::vif_check_cfitsio(status_, "could not reach HDU "+to_string(hdu));

                update_internal_state();
            }

            // Return the number of dimensions of a FITS file
            // Note: will return 0 for FITS tables
            uint_t axis_count() const {
                check_is_open_();

                int naxis;
                fits_get_img_dim(fptr_, &naxis, &status_);
                fits::vif_check_cfitsio(status_, "could not get the number of axis in HDU");
                return naxis;
            }

            // Return the dimensions of a FITS image
            // Note: will return an empty array for FITS tables
            vec1u image_dims() const {
                check_is_open_();

                uint_t naxis = axis_count();
                vec1u dims(naxis);
                if (naxis != 0) {
                    std::vector<long> naxes(naxis);
                    fits_get_img_size(fptr_, naxis, naxes.data(), &status_);
                    fits::vif_check_cfitsio(status_, "could not get image size of HDU");
                    for (uint_t i : range(naxis)) {
                        dims.safe[i] = naxes[naxis-1-i];
                    }
                }

                return dims;
            }

        protected :

            const file_type type_ = generic_file;
            const access_right rights_ = read_write;

            std::string filename_;
            fitsfile* fptr_ = nullptr;
            mutable int status_ = 0;
            fits::table_format format_ = fits::table_format::column_oriented;
        };

        class output_file_base : public virtual file_base {
        public :
            output_file_base(file_type type, access_right rights) : file_base(type, rights) {}

            output_file_base(file_type type, const std::string& filename, access_right rights) :
                file_base(type, filename, rights) {}

            output_file_base(output_file_base&& in) = default;
            ~output_file_base() = default;

            void flush() {
                if (!fptr_) return;

                fits_flush_file(fptr_, &status_);
                fits::vif_check_cfitsio(status_, "could not flush data for '"+filename_+"'");
            }

            void flush_buffer() {
                if (!fptr_) return;

                // Will be faster than flush(), but only updates data, not all the keywords
                fits_flush_buffer(fptr_, 0, &status_);
                fits::vif_check_cfitsio(status_, "could not flush data for '"+filename_+"'");
            }

            void write_header(const fits::header& hdr) {
                check_is_open_();

                std::size_t nentry = hdr.size()/80 + 1;
                fits_set_hdrsize(fptr_, nentry, &status_);

                for (uint_t i : range(nentry)) {
                    std::string entry = hdr.substr(i*80, std::min(std::size_t(80), hdr.size() - i*80));
                    if (begins_with(entry, "END ")) continue;

                    // Skip if it is an internal FITS keyword
                    std::size_t eqpos = entry.find_first_of("=");
                    if (eqpos != entry.npos) {
                        std::string nam = trim(entry.substr(0, eqpos));
                        if (nam == "SIMPLE" || nam == "BITPIX" || begins_with(nam, "NAXIS") ||
                            nam == "EXTEND" || nam == "XTENSION" || nam == "EXTNAME" ||
                            nam == "PCOUNT" || nam == "GCOUNT") {
                            continue;
                        }
                    }

                    fits_write_record(fptr_, entry.c_str(), &status_);
                    fits::vif_check_cfitsio(status_, "could not write header");
                }
            }

            template<typename T, typename enable = typename
                std::enable_if<!is_string<meta::decay_t<T>>::value>::type>
            void write_keyword(const std::string& name, const T& value,
                const std::string& com = "") {
                check_is_open_();

                fits_update_key(fptr_, traits<meta::decay_t<T>>::ttype,
                    name.c_str(), const_cast<T*>(&value),
                    const_cast<char*>(com.c_str()), &status_);
                fits::vif_check_cfitsio(status_, "could not write keyword '"+name+"'");
            }

            template<typename T, typename enable = typename
                std::enable_if<!is_string<meta::decay_t<T>>::value>::type>
            void add_keyword(const std::string& name, const T& value,
                const std::string& com = "") {
                check_is_open_();

                fits_write_key(fptr_, traits<meta::decay_t<T>>::ttype,
                    name.c_str(), const_cast<T*>(&value),
                    const_cast<char*>(com.c_str()), &status_);
                fits::vif_check_cfitsio(status_, "could not write keyword '"+name+"'");
            }

            void write_keyword(const std::string& name, const std::string& value,
                const std::string& com = "") {
                check_is_open_();

                if (name.size() > 8) {

                }

                fits_update_key(fptr_, TSTRING, name.c_str(),
                    const_cast<char*>(value.c_str()),
                    const_cast<char*>(com.c_str()), &status_);
                fits::vif_check_cfitsio(status_, "could not write keyword '"+name+"'");
            }

            void add_keyword(const std::string& name, const std::string& value,
                const std::string& com = "") {
                check_is_open_();

                fits_write_key(fptr_, TSTRING, name.c_str(),
                    const_cast<char*>(value.c_str()),
                    const_cast<char*>(com.c_str()), &status_);
                fits::vif_check_cfitsio(status_, "could not write keyword '"+name+"'");
            }

            void remove_keyword(const std::string& name) {
                check_is_open_();

                fits_delete_key(fptr_, name.c_str(), &status_);
                status_ = 0; // don't care if the keyword doesn't exist or other errors
            }

            void remove_hdu() {
                check_is_open_();

                fits_delete_hdu(fptr_, nullptr, &status_);
                fits::vif_check_cfitsio(status_, "could not remove the current HDU");

                update_internal_state();
            }
        };
    }
}

namespace fits {
    // Header keyword manipulation
    template<typename T>
    bool getkey(const fits::header& hdr, const std::string& key, T& v) {
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i : range(nentry)) {
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

    inline bool getkey(const fits::header& hdr, const std::string& key, std::string& v) {
        std::size_t nentry = hdr.size()/80 + 1;
        for (uint_t i : range(nentry)) {
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
        std::string value;
        if (std::is_same<T,double>::value) {
            value = to_string(format::precision(v, 16));
        } else {
            value = to_string(v);
        }
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
        for (uint_t i : range(nentry)) {
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

    struct header_keyword {
        std::string key, value, comment;
        bool novalue = true;
    };

    using parsed_header = vec<1,header_keyword>;

    inline fits::parsed_header parse_header(const fits::header& hdr) {
        fits::parsed_header keys;
        vec1s lines = cut(hdr, 80);
        keys.reserve(lines.size());
        for (uint_t i : range(lines)) {
            if (lines[i].find("HISTORY ") == 0) {
                keys.push_back(header_keyword{});
                auto& k = keys.back();
                k.key = trim(lines[i]);
            } else if (lines[i].find("CONTINUE ") == 0) {
                auto& k = keys.back();
                std::string right = trim(lines[i].substr(10));
                if (!right.empty() && right[0] == '\'') {
                    auto p2 = right.find_first_of('\'', 1);
                    std::string val;
                    if (p2 != std::string::npos) {
                        val = trim(right.substr(1, p2-1));
                        k.comment = trim(right.substr(p2+1));
                    } else {
                        val = trim(right);
                    }


                    if (!k.value.empty() && k.value[0] == '\'') {
                        k.value = k.value.substr(1, k.value.size()-2);
                    }

                    if (k.value.back() == '&') {
                        k.value = k.value.substr(0, k.value.size()-1);
                    }
                    if (val.back() == '&') {
                        val = val.substr(0, val.size()-1);
                    }

                    k.value = '\''+k.value+val+'\'';
                } else {
                    // Don't know what to do, just ignore
                }
            } else {
                keys.push_back(header_keyword{});
                auto& k = keys.back();

                auto p = lines[i].find_first_of("=/");
                if (p == std::string::npos) {
                    k.key = trim(lines[i]);
                } else if (lines[i][p] == '/') {
                    k.key = trim(lines[i].substr(0, p));
                    k.comment = trim(lines[i].substr(p));
                } else if (lines[i][p] == '=') {
                    k.novalue = false;
                    k.key = trim(lines[i].substr(0, p));

                    std::string right = trim(lines[i].substr(p+1));
                    if (!right.empty()) {
                        if (right[0] == '\'') {
                            auto p2 = right.find_first_of('\'', 1);
                            if (p2 != std::string::npos) {
                                ++p2;

                                k.value = trim(right.substr(0, p2));
                                k.comment = trim(right.substr(p2));
                            } else {
                                k.value = trim(right);
                            }
                        } else {
                            uint_t p2 = right.find_first_of('/');
                            if (p2 != std::string::npos) {
                                k.value = trim(right.substr(0, p2));
                                k.comment = trim(right.substr(p2));
                            } else {
                                k.value = trim(right);
                            }
                        }
                    }
                }
            }
        }

        return keys;
    }

    inline fits::header serialize_header(const fits::parsed_header& keys) {
        fits::header hdr;
        hdr.reserve(80*keys.size());
        for (auto& k : keys) {
            std::string entry;
            if (!k.novalue) {
                if (begins_with(k.key, "HIERARCH ")) {
                    entry = k.key+" = "+k.value;
                } else {
                    std::string val;
                    if (k.value[0] == '\'') {
                        val = align_left(k.value, 20, ' ');
                    } else {
                        val = align_right(k.value, 20, ' ');
                    }
                    entry = align_left(k.key, 8, ' ')+"= "+val;
                }
            } else {
                entry = align_left(k.key, 30, ' ');
            }

            if (!k.comment.empty()) {
                entry += " "+k.comment;
            }

            if (entry.size() > 80) {
                entry = entry.substr(0, 80);
            } else if (entry.size() < 80) {
                entry += std::string(80-entry.size(), ' ');
            }

            hdr += entry;
        }

        return hdr;
    }
}
}

#else

// Minimal support
namespace vif {
namespace fits {
    using header = std::string;
}
}

#endif

#endif
