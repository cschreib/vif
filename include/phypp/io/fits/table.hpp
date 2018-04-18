#ifndef PHYPP_IO_FITS_TABLE_HPP
#define PHYPP_IO_FITS_TABLE_HPP

#include "phypp/reflex/reflex_helpers.hpp"
#include "phypp/io/fits/base.hpp"
#include "phypp/math/reduce.hpp"

namespace phypp {
namespace fits {
    // Reading options
    struct table_read_options {
        bool allow_narrow = false;
        bool allow_missing = false;
        bool allow_dim_promote = false;
        uint_t first_row = 0;
        uint_t last_row = npos;

        table_read_options operator | (const table_read_options& o) const {
            table_read_options n = *this;
            if (o.allow_narrow)      n.allow_narrow = true;
            if (o.allow_missing)     n.allow_missing = true;
            if (o.allow_dim_promote) n.allow_dim_promote = true;
            n.first_row = std::min(first_row, o.first_row);
            n.last_row = std::max(last_row, o.last_row);
            return n;
        }
    };

    const table_read_options narrow = []() {
        table_read_options opts;
        opts.allow_narrow = true;
        return opts;
    }();

    const table_read_options missing = []() {
        table_read_options opts;
        opts.allow_missing = true;
        return opts;
    }();

    const table_read_options dim_promote = []() {
        table_read_options opts;
        opts.allow_dim_promote = true;
        return opts;
    }();

    const table_read_options permissive = []() {
        table_read_options opts;
        opts.allow_narrow = true;
        opts.allow_missing = true;
        opts.allow_dim_promote = true;
        return opts;
    }();

    const auto row = [](uint_t r) {
        table_read_options opts;
        opts.first_row = r;
        opts.last_row = r+1;
        return opts;
    };

    const auto rows = [](uint_t r1, uint_t r2) {
        table_read_options opts;
        opts.first_row = r1;
        opts.last_row = r2;
        return opts;
    };
}

namespace impl {
    namespace fits_impl {
        template<typename T>
        struct data_type {
            using type = T;
        };

        template<std::size_t D, typename T>
        struct data_type<vec<D,T>> {
            using type = typename vec<D,T>::rtype;
        };


        template<typename T>
        struct data_dim {
            static constexpr const uint_t value = 1;
        };

        template<std::size_t D, typename T>
        struct data_dim<vec<D,T>> {
            static constexpr const uint_t value = D;
        };

        template<std::size_t D>
        struct data_dim<vec<D,std::string>> {
            static constexpr const uint_t value = D+1;
        };

        // Traits to identify types that can be read from a FITS table
        template<typename T>
        struct is_readable_column_type : meta::is_any_type_of<T, meta::type_list<
            bool, char, int_t, uint_t, float, double, std::string
        >> {};

        template<std::size_t D, typename T>
        struct is_readable_column_type<vec<D,T>> : meta::is_any_type_of<T, meta::type_list<
            bool, char, int_t, uint_t, float, double, std::string
        >> {};

        template<typename T>
        struct is_readable_column_type<reflex::struct_t<T>> : std::true_type {};

        template<std::size_t D, typename T>
        struct is_readable_column_type<vec<D,T*>> : std::false_type {};

        template<typename T>
        struct is_readable_column_type<impl::named_t<T>> : is_readable_column_type<meta::decay_t<T>> {};
    }
}

namespace fits {
    // Data type to describe the content of a column
    struct column_info {
        std::string name;
        enum type_t {
            string, boolean, byte, integer, float_simple, float_double
        } type;
        vec1u dims;
        uint_t length = 1;
    };

    // FITS input table (read only)
    class input_table : public virtual impl::fits_impl::file_base {
    public :
        explicit input_table(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::table_file, filename, impl::fits_impl::read_only) {}
        explicit input_table(const std::string& filename, uint_t hdu) :
            impl::fits_impl::file_base(impl::fits_impl::table_file, filename, impl::fits_impl::read_only) {
            reach_hdu(hdu);
        }

        input_table(input_table&&) noexcept = default;
        input_table(const input_table&) noexcept = delete;
        input_table& operator = (input_table&&) noexcept = delete;
        input_table& operator = (const input_table&&) noexcept = delete;

        vec<1,column_info> read_column_info() const {
            vec<1,column_info> cols;

            fits::header hdr = read_header();

            int ncol;
            fits_get_num_cols(fptr_, &ncol, &status_);
            fits::phypp_check_cfitsio(status_, "could not read number of columns in HDU");

            for (uint_t c : range(ncol)) {
                column_info ci;

                char name[80];
                char type = 0;
                long repeat = 0;
                fits_get_bcolparms(fptr_, c+1, name, nullptr, &type, &repeat,
                    nullptr, nullptr, nullptr, nullptr, &status_);
                fits::phypp_check_cfitsio(status_, "could not parameters of column "+to_string(c+1)+" in HDU");

                const uint_t max_dim = 256;
                long axes[max_dim];
                int naxis = 0;
                fits_read_tdim(fptr_, c+1, max_dim, &naxis, axes, &status_);
                fits::phypp_check_cfitsio(status_, "could not parameters of column '"+std::string(name)+"'");

                if (!getkey(hdr, "TTYPE"+to_string(c+1), ci.name)) {
                    continue;
                }

                ci.name = trim(ci.name);

                switch (type) {
                case 'A' : {
                    ci.type = column_info::string;
                    break;
                }
                case 'B' : {
                    vec<1,char> v(repeat);
                    char def = impl::fits_impl::traits<char>::def();
                    int null;
                    fits_read_col(fptr_, impl::fits_impl::traits<char>::ttype, c+1, 1, 1, repeat,
                        &def, v.data.data(), &null, &status_);
                    fits::phypp_check_cfitsio(status_, "could not read column '"+std::string(name)+"'");

                    if (min(v) == 0 && max(v) <= 1) {
                        ci.type = column_info::boolean;
                    } else {
                        ci.type = column_info::byte;
                    }
                    break;
                }
                case 'S' : {
                    ci.type = column_info::byte;
                    break;
                }
                case 'K' :
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
                            ci.dims.push_back(1);
                            break;
                        }

                        continue;
                    }

                    ci.dims.push_back(axes[i]);
                }

                ci.dims = reverse(ci.dims);

                cols.push_back(ci);
            }

            return cols;
        }

        class read_sentry {
            friend input_table;

            mutable bool checked = false;
            bool good = true;
            int status = 0;
            std::string errmsg;

            read_sentry(const read_sentry&) = delete;
            read_sentry& operator=(const read_sentry&) = delete;

            read_sentry() = default;

            explicit read_sentry(const input_table* parent, std::string err) {
                status = parent->cfitsio_status();
                good = false;

                if (status != 0) {
                    char txt[FLEN_STATUS];
                    fits_get_errstatus(status, txt);
                    errmsg = "error: cfitsio: "+std::string(txt)+"\nerror: "+err;
                } else {
                    errmsg = std::move(err);
                }

                errmsg += "\nnote: reading '"+parent->filename()+"'";
            }

        public :
            ~read_sentry() {
                if (!checked) {
                    phypp_check(good, errmsg);
                }
            }

            read_sentry(read_sentry&& s) noexcept :
                checked(s.checked), good(s.good), errmsg(std::move(s.errmsg)) {
                s.checked = true;
            }

            read_sentry& operator=(read_sentry&& s) noexcept {
                checked = s.checked;
                good = s.good;
                errmsg = std::move(s.errmsg);
                s.checked = true;
                return *this;
            }

            bool is_good() const {
                checked = true;
                return good;
            }

            const std::string& reason() const {
                checked = true;
                return errmsg;
            }

            explicit operator bool() const {
                return is_good();
            }
        };

        static constexpr uint_t max_column_dims = 256;

    private :

        template<std::size_t Dim, typename Type,
            typename enable = typename std::enable_if<!std::is_same<Type,std::string>::value>::type>
        void read_column_impl_(const table_read_options& opts, vec<Dim,Type>& v,
            const std::string& cname, int cid, long naxis, const std::array<long,max_column_dims>& naxes,
            long repeat, long nrow, bool colfits) const {

            if (v.empty()) return;

            long firstrow = 1, firstelem = 1;
            if (colfits) {
                firstelem = opts.first_row*(repeat/naxes[naxis-1]) + 1;
            } else {
                firstrow = opts.first_row + 1;
            }

            long nelem = v.size();
            Type def = impl::fits_impl::traits<Type>::def();
            int null;
            fits_read_col(
                fptr_, impl::fits_impl::traits<Type>::ttype, cid, firstrow, firstelem, nelem, &def,
                v.data.data(), &null, &status_
            );
            fits::phypp_check_cfitsio(status_, "could not read column '"+cname+"'");
        }

        template<typename Type>
        void read_column_impl_(const table_read_options&, Type& v, const std::string& cname,  int cid,
            long naxis, const std::array<long,max_column_dims>& naxes, long, long, bool) const {

            Type def = impl::fits_impl::traits<Type>::def();
            int null;
            fits_read_col(
                fptr_, impl::fits_impl::traits<Type>::ttype, cid, 1, 1, 1, &def,
                reinterpret_cast<typename impl::fits_impl::traits<Type>::dtype*>(&v), &null, &status_
            );
            fits::phypp_check_cfitsio(status_, "could not read column '"+cname+"'");
        }

        template<std::size_t Dim>
        void read_column_impl_(const table_read_options& opts, vec<Dim,std::string>& v,
            const std::string& cname, int cid, long naxis, const std::array<long,max_column_dims>& naxes,
            long repeat, long nrow, bool colfits) const {

            if (v.empty()) return;

            // NB: cfitsio doesn't seem to like reading empty strings
            if (naxes[0] == 0) {
                return;
            }

            long firstrow = 1, firstelem = 1;
            if (colfits) {
                firstelem = opts.first_row*(repeat/naxes[naxis-1]/naxes[0]) + 1;
            } else {
                firstrow = opts.first_row + 1;
            }

            char** buffer = new char*[v.size()];
            for (uint_t i : range(v)) {
                buffer[i] = new char[naxes[0]+1];
            }

            long nelem = v.size();

            char def = '\0';
            int null;
            fits_read_col(
                fptr_, impl::fits_impl::traits<std::string>::ttype, cid, firstrow, firstelem, nelem, &def,
                buffer, &null, &status_
            );
            fits::phypp_check_cfitsio(status_, "could not read column '"+cname+"'");

            for (uint_t i : range(v)) {
                v.safe[i] = trim(std::string(buffer[i]));
                delete[] buffer[i];
            }

            delete[] buffer;
        }

        void read_column_impl_(const table_read_options&, std::string& v, const std::string& cname,
            int cid,long naxis, const std::array<long,max_column_dims>& naxes, long repeat, long nrow,
            bool colfits) const {

            // NB: cfitsio doesn't seem to like reading empty strings
            if (repeat == 0) {
                v.clear();
                return;
            }

            char** buffer = new char*[1];
            buffer[0] = new char[repeat+1];

            char def = '\0';
            int null;

            // Hack to make it work...
            fptr_->Fptr->tableptr[cid-1].twidth = repeat;

            fits_read_col(
                fptr_, TSTRING, cid, 1, 1, 1, &def,
                buffer, &null, &status_
            );
            fits::phypp_check_cfitsio(status_, "could not read column '"+cname+"'");

            v = trim(std::string(buffer[0]));
            delete[] buffer[0];
            delete[] buffer;
        }

        template<std::size_t Dim, typename Type,
            typename enable = typename std::enable_if<!std::is_same<Type,std::string>::value>::type>
        void read_column_resize_(vec<Dim,Type>& v, long naxis,
            const std::array<long,max_column_dims>& naxes) const {

            for (uint_t i : range(Dim)) {
                if (i < (Dim-naxis)) {
                    v.dims[i] = 1;
                } else {
                    v.dims[i] = naxes[Dim-1-i];
                }
            }

            v.resize();
        }

        template<typename Type>
        void read_column_resize_(Type& v, long naxis,
            const std::array<long,max_column_dims>& naxes) const {}

        template<std::size_t Dim>
        void read_column_resize_(vec<Dim,std::string>& v, long naxis,
            const std::array<long,max_column_dims>& naxes) const {

            for (uint_t i : range(Dim)) {
                if (i < (Dim+1-naxis)) {
                    v.dims[i] = 1;
                } else {
                    v.dims[i] = naxes[Dim-i];
                }
            }

            v.resize();
        }

        void read_column_resize_(std::string& v, long naxis,
            const std::array<long,max_column_dims>& naxes) const {}

        template<typename T>
        bool read_column_check_type_(const table_read_options& opts, int type) const {
            if (opts.allow_narrow) {
                return impl::fits_impl::traits<T>::is_convertible_narrow(type);
            } else {
                return impl::fits_impl::traits<T>::is_convertible(type);
            }
        }

        template<std::size_t Dim, typename T, typename enable =
            typename std::enable_if<!std::is_same<T,std::string>::value>::type>
        bool read_column_check_dim_(const table_read_options& opts, vec<Dim,T>&, uint_t vdim,
            long naxis, long repeat) const {

            if (opts.allow_dim_promote) {
                return uint_t(naxis) <= vdim;
            } else {
                return uint_t(naxis) == vdim;
            }
        }

        template<typename T>
        bool read_column_check_dim_(const table_read_options& opts, T&, uint_t vdim,
            long naxis, long repeat) const {

            // Early exit to refuse loading 1D vectors into scalars if the number of
            // element is greater than one
            if (naxis == 1 && repeat > 1) return false;

            if (opts.allow_dim_promote) {
                return uint_t(naxis) <= vdim;
            } else {
                return uint_t(naxis) == vdim;
            }
        }

        template<std::size_t Dim>
        bool read_column_check_dim_(const table_read_options& opts, vec<Dim,std::string>&,
            uint_t vdim, long naxis, long repeat) const {

            // Early exit to allow loading 1D byte into vec1s
            if (Dim == 1 && naxis == 1) return true;

            if (opts.allow_dim_promote) {
                return uint_t(naxis) <= vdim;
            } else {
                return uint_t(naxis) == vdim;
            }
        }

        bool read_column_check_dim_(const table_read_options& opts, std::string&, uint_t vdim,
            long naxis, long repeat) const {

            if (opts.allow_dim_promote) {
                return uint_t(naxis) <= vdim;
            } else {
                return uint_t(naxis) == vdim;
            }
        }

        bool read_column_check_rows_(table_read_options& opts,
            std::array<long,max_column_dims>& axes, long naxis) const {
            if (opts.last_row == npos) {
                opts.last_row = axes[naxis-1];
            }

            if (opts.first_row == npos) {
                opts.first_row = 0;
            }

            if (opts.last_row < opts.first_row) {
                std::swap(opts.last_row, opts.first_row);
            }

            if (opts.last_row <= uint_t(axes[naxis-1])) {
                axes[naxis-1] = opts.last_row - opts.first_row;
                return true;
            } else {
                return false;
            }
        }

        struct do_read_struct_ {
            const input_table* tbl;
            const table_read_options& opts;
            std::string base;

            template<typename P>
            void operator () (reflex::member_t& m, P&& v) {
                tbl->read_column(opts, base+to_upper(m.name), std::forward<P>(v));
            }
        };

        template<typename T>
        read_sentry read_column_(table_read_options opts,
            const std::string& tcolname, T& value, std::false_type) const {

            static_assert(impl::fits_impl::is_readable_column_type<typename std::decay<T>::type>::value,
                "this value cannot be read from a FITS file");

            // Check if column exists
            std::string colname = to_upper(tcolname);
            int cid;
            fits_get_colnum(fptr_, CASEINSEN, const_cast<char*>(colname.c_str()), &cid, &status_);
            if (status_ != 0) {
                status_ = 0;
                fits_clear_errmsg();
                if (opts.allow_missing) {
                    return read_sentry{};
                } else {
                    return read_sentry{this, "cannot find collumn '"+colname+"'"};
                }
            }

            // Collect data on the output type
            using vtype = typename impl::fits_impl::data_type<T>::type;
            const uint_t vdim = impl::fits_impl::data_dim<T>::value;

            // Check if type match
            int type;
            long repeat, width;
            fits_get_coltype(fptr_, cid, &type, &repeat, &width, &status_);
            fits::phypp_check_cfitsio(status_, "could not read type of column '"+tcolname+"'");
            if (!read_column_check_type_<vtype>(opts, type)) {
                return read_sentry{this, "wrong type for column '"+colname+"' "
                    "(expected "+pretty_type_t(vtype)+", got "+impl::fits_impl::type_to_string_(type)+")"};
            }

            // Check if dimensions match
            std::array<long,max_column_dims> axes;
            int naxis = 0;
            fits_read_tdim(fptr_, cid, max_column_dims, &naxis, axes.data(), &status_);
            fits::phypp_check_cfitsio(status_, "could not read dimensions of column '"+tcolname+"'");

            // Support ASCII tables with string columns
            std::string extension;
            if (read_keyword("XTENSION", extension) && extension == "TABLE") {
                std::string tform;
                if (read_keyword("TFORM"+to_string(cid), tform) && tform[0] == 'A') {
                    from_string(erase_begin(tform, "A"), axes[0]);
                }
            }

            // Support loading row-oriented FITS tables (and ASCII)
            uint_t nrow = 1;
            bool colfits = true;
            if (read_keyword("NAXIS2", nrow) && nrow > 1) {
                colfits = false;
                if (naxis == 1 && axes[naxis-1] == 1) {
                    axes[naxis-1] = nrow;
                } else {
                    axes[naxis] = nrow;
                    ++naxis;
                }
            }

            if (!read_column_check_dim_(opts, value, vdim, naxis, repeat)) {
                return read_sentry{this, "wrong dimension for column '"+colname+"' "
                    "(expected "+to_string(vdim)+", got "+to_string(naxis)+")"};
            }

            auto raxes = axes;
            if (!read_column_check_rows_(opts, raxes, naxis)) {
                std::string saxes = "{";
                for (uint_t i : range(naxis)) {
                    saxes += (i != 0 ? ", " : "") + to_string(axes[i]);
                }
                saxes += "}";

                return read_sentry{this, "wrong dimension for column '"+colname+"' "
                    "(asking for "+(opts.first_row == opts.last_row ?
                        "row "+to_string(opts.first_row) :
                        "rows "+
                        (opts.first_row == npos ? "0" : to_string(opts.first_row))
                        +" to "+
                        (opts.last_row == npos ? to_string(axes[naxis-1]) : to_string(opts.last_row))
                    )+", in dimensions "+saxes+")"};
            }

            // Resize vector
            read_column_resize_(value, naxis, raxes);

            // Read
            if (nrow != 0) {
                read_column_impl_(opts, value, tcolname, cid, naxis, axes, repeat, nrow, colfits);
            }

            return read_sentry{};
        }

        template<typename T>
        read_sentry read_column_(const table_read_options& opts,
            const std::string& colname, reflex::struct_t<T> value, std::true_type) const {

            #ifdef NO_REFLECTION
            static_assert(!std::is_same<T,T>::value,
                "this function requires reflection capabilities (NO_REFLECTION=0)");
            #endif

            do_read_struct_ run{this, opts, to_upper(colname)+"."};
            reflex::foreach_member(value, run);

            return read_sentry{};
        }

        template<typename T>
        read_sentry read_column_(const table_read_options& opts,
            const std::string& colname, T& value, std::true_type) const {

            #ifdef NO_REFLECTION
            static_assert(!std::is_same<T,T>::value,
                "this function requires reflection capabilities (NO_REFLECTION=0)");
            #endif

            do_read_struct_ run{this, opts, to_upper(colname)+"."};
            reflex::foreach_member(reflex::wrap(value), run);

            return read_sentry{};
        }

    public :

        template<typename T>
        read_sentry read_column(const table_read_options& opts,
            const std::string& tcolname, T&& value) const {
            return read_column_(opts, tcolname, std::forward<T>(value), reflex::enabled<meta::decay_t<T>>{});
        }

        template<typename T>
        read_sentry read_column(const std::string& tcolname, T&& value) const {
            return read_column(table_read_options{}, tcolname, std::forward<T>(value));
        }

    private :

        void read_columns_impl_(const table_read_options&) const {
            // Nothing more to do
        }

        template<typename T, typename ... Args>
        void read_columns_impl_(const table_read_options& opts,
            const std::string& tcolname, T& value, Args&& ... args) const {

            read_column(opts, tcolname, value);
            read_columns_impl_(opts, std::forward<Args>(args)...);
        }

    public :

        template<typename ... Args, typename enable = typename std::enable_if<
            !std::is_same<typename std::decay<meta::first_type<meta::type_list<Args...>>>::type, impl::ascii_impl::macroed_t>::value &&
            (sizeof...(Args) > 1)>::type>
        void read_columns(const table_read_options& opts, Args&& ... args) const {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;
            using filtered_first  = meta::filter_type_list<meta::bool_list<true,false>, arg_list>;
            using filtered_second = meta::filter_type_list<meta::bool_list<false,true>, arg_list>;

            static_assert(
                meta::are_all_true<meta::binary_first_apply_type_to_bool_list<
                filtered_first, std::is_convertible, std::string>>::value &&
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                filtered_second, impl::fits_impl::is_readable_column_type>>::value,
                "arguments must be a sequence of 'column name', 'readable value'");

            // Read
            read_columns_impl_(opts, std::forward<Args>(args)...);
        }

        template<typename ... Args, typename enable = typename std::enable_if<
            !std::is_same<typename std::decay<meta::first_type<meta::type_list<Args...>>>::type, impl::ascii_impl::macroed_t>::value &&
            !std::is_same<typename std::decay<meta::first_type<meta::type_list<Args...>>>::type, table_read_options>::value &&
            (sizeof...(Args) > 1)>::type>
        void read_columns(Args&& ... args) const {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;
            using filtered_first  = meta::filter_type_list<meta::bool_list<true,false>, arg_list>;
            using filtered_second = meta::filter_type_list<meta::bool_list<false,true>, arg_list>;

            static_assert(
                meta::are_all_true<meta::binary_first_apply_type_to_bool_list<
                filtered_first, std::is_convertible, std::string>>::value &&
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                filtered_second, impl::fits_impl::is_readable_column_type>>::value,
                "arguments must be a sequence of 'column name', 'readable value'");

            // Read
            read_columns_impl_(table_read_options{}, std::forward<Args>(args)...);
        }

    private :

        void read_columns_impl_(const table_read_options&, impl::ascii_impl::macroed_t, const std::string&) const {
            // Nothing more to do
        }

        template<typename T, typename ... Args>
        void read_columns_impl_(const table_read_options& opts, impl::ascii_impl::macroed_t,
            std::string names, T& value, Args&& ... args) const {

            std::string tcolname = impl::ascii_impl::pop_macroed_name(names);
            read_column(opts, impl::ascii_impl::bake_macroed_name(tcolname), value);
            read_columns_impl_(opts, impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }
        template<typename T, typename ... Args>
        void read_columns_impl_(const table_read_options& opts, impl::ascii_impl::macroed_t,
            std::string names, const impl::named_t<T>& value, Args&& ... args) const {

            impl::ascii_impl::pop_macroed_name(names);
            read_column(opts, value.name, value.obj);
            read_columns_impl_(opts, impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

        template<typename ... Args>
        void read_columns_impl_(impl::ascii_impl::macroed_t, const std::string& names, Args&& ... args) const {
            read_columns_impl_(table_read_options{}, impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    public :

        template<typename ... Args>
        void read_columns(const table_read_options& opts, impl::ascii_impl::macroed_t, const std::string& names,
            Args&& ... args) const {

            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;

            static_assert(
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                arg_list, impl::fits_impl::is_readable_column_type>>::value,
                "arguments must be a sequence of readable values");

            // Read
            read_columns_impl_(opts, impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

        template<typename ... Args>
        void read_columns(impl::ascii_impl::macroed_t, const std::string& names, Args&& ... args) const {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;

            static_assert(
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                arg_list, impl::fits_impl::is_readable_column_type>>::value,
                "arguments must be a sequence of readable values");

            // Read
            read_columns_impl_(table_read_options{}, impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    public :

        template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
        void read_columns(const table_read_options& opts, T& t) {
            reflex::foreach_member(reflex::wrap(t), do_read_struct_{this, opts, ""});
        }

        template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
        void read_columns(T& t) {
            read_columns(table_read_options{}, t);
        }
    };
}

    // Traits to identify types that can be written into a FITS table
namespace impl {
    namespace fits_impl {
        template<typename T>
        struct is_writable_column_type : meta::is_any_type_of<T, meta::type_list<
            bool, char, int_t, uint_t, float, double, std::string
        >> {};

        template<std::size_t D, typename T>
        struct is_writable_column_type<vec<D,T>> : meta::is_any_type_of<typename vec<D,T>::rtype, meta::type_list<
            bool, char, int_t, uint_t, float, double, std::string
        >> {};

        template<typename T>
        struct is_writable_column_type<reflex::struct_t<T>> : std::true_type {};

        template<typename T>
        struct is_writable_column_type<impl::named_t<T>> : is_writable_column_type<meta::decay_t<T>> {};
    }
}

namespace fits {
    // Output FITS table (write only, overwrites existing files)
    class output_table : public impl::fits_impl::output_file_base {
    public :

        explicit output_table(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::table_file, filename, impl::fits_impl::write_only),
            impl::fits_impl::output_file_base(impl::fits_impl::table_file, filename, impl::fits_impl::write_only) {}

        output_table(output_table&&) noexcept = default;
        output_table(const output_table&) = delete;
        output_table& operator = (output_table&&) noexcept = delete;
        output_table& operator = (const output_table&&) = delete;

    protected :

        void create_table_() {
            uint_t nhdu = hdu_count();
            if (nhdu > 0) {
                // There are HDUs in this file already
                auto type = hdu_type();
                if (type != fits::table_hdu) {
                    // The current HDU is not a table
                    phypp_check(type == fits::null_hdu || type == fits::empty_hdu,
                        "cannot create table, there is already data in this HDU");

                    // It is empty, so delete it and create a new table extension
                    uint_t hdu = current_hdu();
                    if (hdu == nhdu-1) {
                        // We are at the last HDU, so just delete and create
                        fits_delete_hdu(fptr_, nullptr, &status_);
                        fits::phypp_check_cfitsio(status_, "could not create table HDU");
                        fits_insert_btbl(fptr_, 1, 0, 0, 0, 0, nullptr, 0, &status_);
                        fits::phypp_check_cfitsio(status_, "could not create table HDU");
                    } else {
                        phypp_check(hdu != 0, "cannot write tables in the primary HDU");

                        // Delete current
                        fits_delete_hdu(fptr_, nullptr, &status_);
                        fits::phypp_check_cfitsio(status_, "could not create table HDU");
                        // Move back because CFITSIO will insert *after* current HDU
                        fits_movrel_hdu(fptr_, -1, nullptr, &status_);
                        fits::phypp_check_cfitsio(status_, "could not create table HDU");
                        // Add the new extension
                        fits_insert_btbl(fptr_, 1, 0, 0, 0, 0, nullptr, 0, &status_);
                        fits::phypp_check_cfitsio(status_, "could not create table HDU");

                    }
                }
            } else {
                // No HDU yet, just create the primary array and a table extension
                fits_insert_btbl(fptr_, 1, 0, 0, 0, 0, nullptr, 0, &status_);
                fits::phypp_check_cfitsio(status_, "could not create table HDU");
            }
        }

    private :

        template<std::size_t Dim, typename Type,
            typename enable = typename std::enable_if<
                !std::is_same<Type,std::string>::value && !std::is_pointer<Type>::value>::type>
        void write_column_impl_(const vec<Dim,Type>& value, const std::array<long,Dim>&,
            const std::string& tcolname, int cid) {
            // cfitsio doesn't like writing empty columns
            if (value.empty()) return;

            fits_write_col(fptr_, impl::fits_impl::traits<Type>::ttype, cid, 1, 1, value.size(),
                const_cast<typename vec<Dim,Type>::dtype*>(value.data.data()), &status_);
            fits::phypp_check_cfitsio(status_, "could not write column '"+tcolname+"'");
        }

        template<typename Type>
        void write_column_impl_(const Type& value, const std::array<long,0>&,
            const std::string& tcolname, int cid) {
            using dtype = typename impl::fits_impl::traits<Type>::dtype;
            fits_write_col(fptr_, impl::fits_impl::traits<Type>::ttype, cid, 1, 1, 1,
                const_cast<dtype*>(reinterpret_cast<const dtype*>(&value)), &status_);
            fits::phypp_check_cfitsio(status_, "could not write column '"+tcolname+"'");
        }

        template<std::size_t Dim>
        void write_column_impl_(const vec<Dim,std::string>& value,
            const std::array<long,Dim+1>& dims, const std::string& tcolname, int cid) {

            const uint_t nmax = dims[0];

            // cfitsio doesn't like writing empty string or columns
            if (nmax == 0 || value.empty()) return;

            char** buffer = new char*[value.size()];
            for (uint_t i : range(value)) {
                buffer[i] = new char[nmax+1];
                strcpy(buffer[i], value[i].c_str());
            }

            fits_write_col(
                fptr_, impl::fits_impl::traits<std::string>::ttype, cid, 1, 1,
                value.size(), buffer, &status_
            );
            fits::phypp_check_cfitsio(status_, "could not write column '"+tcolname+"'");

            for (uint_t i : range(value)) {
                delete[] buffer[i];
            }

            delete[] buffer;
        }

        void write_column_impl_(const std::string& value, const std::array<long,1>&,
            const std::string& tcolname, int cid) {

            // NB: for some reason cfitsio crashes on empty string
            if (value.empty()) return;

            char** buffer = new char*[1];
            buffer[0] = const_cast<char*>(value.c_str());
            fits_write_col(fptr_, impl::fits_impl::traits<std::string>::ttype, cid, 1, 1, 1,
                buffer, &status_);
            fits::phypp_check_cfitsio(status_, "could not write column '"+tcolname+"'");

            delete[] buffer;
        }

        template<std::size_t Dim, typename Type, std::size_t DD>
        void write_column_impl_(const vec<Dim,Type*>& value, const std::array<long,DD>& dims,
            const std::string& tcolname, int cid) {
            write_column_impl_(value.concretise(), dims, tcolname, cid);
        }

        template<std::size_t Dim, typename Type>
        std::array<long,Dim> write_column_get_dims_(const vec<Dim,Type>& value) {
            std::array<long,Dim> dims;
            for (uint_t i : range(Dim)) {
                dims[i] = value.dims[Dim-1-i];
            }

            return dims;
        }

        template<typename T>
        std::array<long,0> write_column_get_dims_(const T&) {
            return {};
        }

        template<std::size_t Dim, typename T>
        std::array<long,Dim+1> write_column_get_dims_stringv_(const vec<Dim,T>& value) {
            std::array<long,Dim+1> dims;

            std::size_t nmax = 0;
            for (auto& s : value) {
                if (s.size() > nmax) {
                    nmax = s.size();
                }
            }

            dims[0] = nmax;
            for (uint_t i : range(Dim)) {
                dims[i+1] = value.dims[Dim-1-i];
            }

            return dims;
        }

        template<std::size_t Dim>
        std::array<long,Dim+1> write_column_get_dims_(const vec<Dim,std::string>& value) {
            return write_column_get_dims_stringv_(value);
        }

        template<std::size_t Dim>
        std::array<long,Dim+1> write_column_get_dims_(const vec<Dim,std::string*>& value) {
            return write_column_get_dims_stringv_(value);
        }

        std::array<long,1> write_column_get_dims_(const std::string& value) {
            return {{long(value.size())}};
        }

        template<typename Type, std::size_t Dim>
        std::string write_column_make_tform_(meta::type_list<Type>, const std::array<long,Dim>& dims) {
            uint_t size = 1;
            for (uint_t d : dims) {
                size *= d;
            }

            return to_string(size)+impl::fits_impl::traits<Type>::tform;
        }

        template<std::size_t Dim>
        std::string write_column_make_tform_(meta::type_list<std::string>, const std::array<long,Dim>& dims) {
            uint_t size = 1;
            for (uint_t d : dims) {
                size *= d;
            }

            return to_string(size)+impl::fits_impl::traits<std::string>::tform+to_string(dims[0]);
        }


        template<std::size_t Dim>
        void write_column_write_tdim_(const std::array<long,Dim>& dims, const std::string& tcolname,
            int cid) {
            fits_write_tdim(fptr_, cid, Dim, const_cast<long*>(dims.data()), &status_);
            fits::phypp_check_cfitsio(status_, "could not write TDIM for column '"+tcolname+
                "' (dims "+to_string(dims)+")");
        }

        void write_column_write_tdim_(const std::array<long,0>&, const std::string&, int) {
            // Nothing to do
        }

        struct do_write_struct_ {
            output_table* tbl;
            std::string base;

            template<typename P>
            void operator () (const reflex::member_t& m, P&& v) {
                tbl->write_column(base+to_upper(m.name), std::forward<P>(v));
            }
        };

        template<typename T>
        void write_column_(const std::string& tcolname, const T& value, std::false_type) {
            static_assert(impl::fits_impl::is_writable_column_type<typename std::decay<T>::type>::value,
                "this variable cannot be written into a FITS file");

            int cid;
            fits_get_num_cols(fptr_, &cid, &status_);
            fits::phypp_check_cfitsio(status_, "could not get number of columns in HDU");
            ++cid;

            // Collect data on the input type
            using vtype = typename impl::fits_impl::data_type<T>::type;
            const auto dims = write_column_get_dims_(value);

            // Create empty column
            std::string colname = to_upper(tcolname);
            std::string tform = write_column_make_tform_(meta::type_list<vtype>{}, dims);
            fits_insert_col(
                fptr_, cid, const_cast<char*>(colname.c_str()),
                const_cast<char*>(tform.c_str()), &status_
            );
            fits::phypp_check_cfitsio(status_, "could not create column '"+tcolname+"'");

            // Set TDIM keyword (if needed)
            write_column_write_tdim_(dims, tcolname, cid);

            // Write
            write_column_impl_(value, dims, tcolname, cid);
        }

        template<typename T>
        void write_column_(const std::string& colname, reflex::struct_t<T> value, std::true_type) {
            #ifdef NO_REFLECTION
            static_assert(!std::is_same<T,T>::value,
                "this function requires reflection capabilities (NO_REFLECTION=0)");
            #endif

            do_write_struct_ run{this, to_upper(colname)+"."};
            reflex::foreach_member(value, run);
        }

        template<typename T>
        void write_column_(const std::string& colname, const T& value, std::true_type) {
            #ifdef NO_REFLECTION
            static_assert(!std::is_same<T,T>::value,
                "this function requires reflection capabilities (NO_REFLECTION=0)");
            #endif

            do_write_struct_ run{this, to_upper(colname)+"."};
            reflex::foreach_member(reflex::wrap(value), run);
        }

    public :

        template<typename T>
        void write_column(const std::string& colname, T&& value) {
            create_table_();
            write_column_(colname, std::forward<T>(value), reflex::enabled<meta::decay_t<T>>{});
        }

    private :

        void write_columns_impl_() {
            // Nothing to do
        }

        template<typename T, typename ... Args>
        void write_columns_impl_(const std::string& tcolname, const T& value, Args&& ... args) {
            write_column(tcolname, value);
            write_columns_impl_(std::forward<Args>(args)...);
        }

        void write_columns_impl_(impl::ascii_impl::macroed_t, std::string) {
            // Nothing to do
        }

        template<typename T, typename ... Args>
        void write_columns_impl_(impl::ascii_impl::macroed_t, std::string names, T& value, Args&& ... args) {
            std::string tcolname = impl::ascii_impl::pop_macroed_name(names);
            write_column(impl::ascii_impl::bake_macroed_name(tcolname), value);
            write_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

        template<typename T, typename ... Args>
        void write_columns_impl_(impl::ascii_impl::macroed_t, std::string names, const phypp::impl::named_t<T>& value,
            Args&& ... args) {

            impl::ascii_impl::pop_macroed_name(names);
            write_column(value.name, value.obj);
            write_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    public :

        template<typename ... Args, typename enable = typename std::enable_if<
            !std::is_same<typename std::decay<meta::first_type<meta::type_list<Args...>>>::type, impl::ascii_impl::macroed_t>::value &&
            (sizeof...(Args) > 1)>::type>
        void write_columns(Args&& ... args) {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;
            using filtered_first  = meta::filter_type_list<meta::bool_list<true,false>, arg_list>;
            using filtered_second = meta::filter_type_list<meta::bool_list<false,true>, arg_list>;

            static_assert(
                meta::are_all_true<meta::binary_first_apply_type_to_bool_list<
                filtered_first, std::is_convertible, std::string>>::value &&
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                filtered_second, impl::fits_impl::is_writable_column_type>>::value,
                "arguments must be a sequence of 'column name', 'writable value'");

            // Write
            create_table_();
            write_columns_impl_(std::forward<Args>(args)...);
        }

        template<typename ... Args>
        void write_columns(impl::ascii_impl::macroed_t, const std::string& names, Args&& ... args) {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;

            static_assert(
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                arg_list, impl::fits_impl::is_writable_column_type>>::value,
                "arguments must be a sequence of writable values");

            // Write
            create_table_();
            write_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    public :

        template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
        void write_columns(T& t) {
            create_table_();
            reflex::foreach_member(reflex::wrap(t), do_write_struct_{this, ""});
        }
    };

    // Input/output FITS table (read & write, modifies existing files)
    class table : public output_table, public input_table {
    public :
        explicit table(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::table_file, filename, impl::fits_impl::read_write),
            output_table(filename), input_table(filename) {}

        explicit table(const std::string& filename, uint_t hdu) :
            impl::fits_impl::file_base(impl::fits_impl::table_file, filename, impl::fits_impl::read_write),
            output_table(filename), input_table(filename) {
            reach_hdu(hdu);
        }

        table(table&&) noexcept = default;
        table(const table&) = delete;
        table& operator = (table&&) noexcept = delete;
        table& operator = (const table&&) = delete;

        void remove_column(const std::string& tcolname) {
            int cid;
            std::string colname = to_upper(tcolname);
            fits_get_colnum(fptr_, CASEINSEN, const_cast<char*>(colname.c_str()), &cid, &status_);
            if (status_ == 0) {
                fits_delete_col(fptr_, cid, &status_);
            }

            status_ = 0; // ignore errors
        }

        template<typename ... Args>
        void remove_columns(Args&& ... args) {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;

            static_assert(
                meta::are_all_true<meta::binary_first_apply_type_to_bool_list<
                arg_list, std::is_convertible, std::string>>::value,
                "arguments must be a sequence of column names");

            (void)std::initializer_list<int>{(remove_column(args), 0)...};
        }

        template<typename T>
        void update_column(const std::string& tcolname, const T& value) {
            create_table_();
            remove_column(tcolname);
            write_column(tcolname, value);
        }

    private :

        void update_columns_impl_() {}

        template<typename T, typename ... Args>
        void update_columns_impl_(const std::string& tcolname, const T& value, Args&& ... args) {
            update_column(tcolname, value);
            update_columns_impl_(std::forward<Args>(args)...);
        }

        void update_columns_impl_(impl::ascii_impl::macroed_t, std::string) {
            // Nothing to do
        }

        template<typename T, typename ... Args>
        void update_columns_impl_(impl::ascii_impl::macroed_t, std::string names, T& value, Args&& ... args) {
            std::string tcolname = impl::ascii_impl::pop_macroed_name(names);
            update_column(impl::ascii_impl::bake_macroed_name(tcolname), value);
            update_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

        template<typename T, typename ... Args>
        void update_columns_impl_(impl::ascii_impl::macroed_t, std::string names, const phypp::impl::named_t<T>& value,
            Args&& ... args) {

            impl::ascii_impl::pop_macroed_name(names);
            update_column(value.name, value.obj);
            update_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    public :

        template<typename ... Args, typename enable = typename std::enable_if<
            !std::is_same<typename std::decay<meta::first_type<meta::type_list<Args...>>>::type, impl::ascii_impl::macroed_t>::value &&
            (sizeof...(Args) > 1)>::type>
        void update_columns(Args&& ... args) {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;
            using filtered_first  = meta::filter_type_list<meta::bool_list<true,false>, arg_list>;
            using filtered_second = meta::filter_type_list<meta::bool_list<false,true>, arg_list>;

            static_assert(
                meta::are_all_true<meta::binary_first_apply_type_to_bool_list<
                filtered_first, std::is_convertible, std::string>>::value &&
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                filtered_second, impl::fits_impl::is_writable_column_type>>::value,
                "arguments must be a sequence of 'column name', 'readable value'");

            create_table_();
            update_columns_impl_(std::forward<Args>(args)...);
        }

        template<typename ... Args>
        void update_columns(impl::ascii_impl::macroed_t, std::string names, Args&& ... args) {
            // Check types of arguments
            using arg_list = meta::type_list<typename std::decay<Args>::type...>;

            static_assert(
                meta::are_all_true<meta::unary_apply_type_to_bool_list<
                arg_list, impl::fits_impl::is_writable_column_type>>::value,
                "arguments must be a sequence of 'column name', 'readable value'");

            create_table_();
            update_columns_impl_(impl::ascii_impl::macroed_t{}, names, std::forward<Args>(args)...);
        }

    private :

        struct do_update_struct_ {
            table* tbl;

            template<typename P>
            void operator () (const reflex::member_t& m, P&& v) {
                tbl->update_column(to_upper(m.name), std::forward<P>(v));
            }
        };

    public :

        template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
        void update_columns(T& t) {
            create_table_();
            reflex::foreach_member(reflex::wrap(t), do_update_struct_{this});
        }
    };
}
}

#endif
