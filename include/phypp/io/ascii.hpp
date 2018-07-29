#ifndef PHYPP_IO_ASCII_HPP
#define PHYPP_IO_ASCII_HPP

#include <fstream>
#include <tuple>
#include <stdexcept>
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/range.hpp"
#include "phypp/math/base.hpp"

namespace phypp {
namespace ascii {
    struct exception : std::runtime_error {
        exception(const std::string& w) : std::runtime_error(w) {}
    };

    namespace input_format {
        struct options {
            bool auto_skip    = true;
            char skip_pattern = '#';
            uint_t skip_first = 0;
            std::string delim = " \t";
            bool delim_single = false;

            options() = default;
            options(bool sk, char sp, uint_t sf, const std::string& d, bool ds) :
                auto_skip(sk), skip_pattern(sp), skip_first(sf), delim(d), delim_single(ds) {}
        };
    }

    template<typename T, typename ... Args>
    auto columns(std::size_t n, T& t, Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }

    template<typename T, typename ... Args>
    auto columns(std::size_t n, const T& t, const Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }
}

namespace impl {
    namespace ascii_impl {
        using placeholder_t = phypp::impl::placeholder_t;

        struct line_splitter_t {
            std::string line;
            uint_t pos = 0;

            std::string delim = " \t";
            bool delim_single = false;

            bool skip_word() {
                if (delim_single) {
                    if (pos == line.npos) return false;

                    if (pos != 0) {
                        pos += delim.size();
                    }

                    uint_t p1 = line.find(delim, pos);
                    if (p1 == line.npos) {
                        pos = p1;
                    } else {
                        pos = p1 + delim.size();
                    }

                    pos = p1;
                } else {
                    uint_t p0 = line.find_first_not_of(delim, pos);
                    if (p0 == line.npos) return false;

                    pos = line.find_first_of(delim, p0);
                }

                return true;
            }

            bool next_word(std::string& sub) {
                if (delim_single) {
                    if (pos == line.npos) return false;

                    if (pos != 0) {
                        pos += delim.size();
                    }

                    uint_t p1 = line.find(delim, pos);
                    if (p1 == line.npos) {
                        sub = line.substr(pos);
                        pos = p1;
                    } else {
                        sub = line.substr(pos, p1-pos);
                        pos = p1 + delim.size();
                    }

                    pos = p1;
                } else {
                    uint_t p0 = line.find_first_not_of(delim, pos);
                    if (p0 == line.npos) return false;

                    pos = line.find_first_of(delim, p0);
                    if (pos == line.npos) {
                        sub = line.substr(p0);
                    } else {
                        sub = line.substr(p0, pos-p0);
                    }
                }

                return true;
            }

            void reset() {
                pos = 0;
            }
        };

        template<typename T>
        struct is_tuple : std::false_type {};
        template<typename ... Args>
        struct is_tuple<std::tuple<Args...>> : std::true_type {};

        template<typename T>
        struct is_other : std::integral_constant<bool,
            !meta::is_vec<T>::value && !std::is_same<T,impl::placeholder_t>::value &&
            !is_tuple<T>::value> {};

        inline void read_table_resize_(std::size_t n) {}

        template<typename T, typename ... Args>
        void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args);
        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args);
        template<typename ... Args>
        void read_table_resize_(std::size_t n, impl::placeholder_t, Args& ... args);
        template<typename T, typename ... Args, typename enable = typename std::enable_if<is_other<T>::value>::type>
        void read_table_resize_(std::size_t n, T&, Args& ... args);

        template<typename T, typename ... Args>
        void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args) {
            v.resize(n);
            read_table_resize_(n, args...);
        }

        inline void read_table_resize_cols_(std::size_t n, std::size_t m) {}

        template<typename Type, typename ... VArgs>
        void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args);
        template<typename ... VArgs>
        void read_table_resize_cols_(std::size_t n, std::size_t m, impl::placeholder_t, VArgs&... args);

        template<typename Type, typename ... VArgs>
        void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args) {
            v.resize(n, m);
            read_table_resize_cols_(n, m, args...);
        }

        template<typename ... VArgs>
        void read_table_resize_cols_(std::size_t n, std::size_t m, impl::placeholder_t, VArgs&... args) {
            read_table_resize_cols_(n, m, args...);
        }

        template<typename U, typename ... VArgs, std::size_t ... S>
        void read_table_resize_cols_i_(std::size_t n, std::size_t m, std::tuple<U,VArgs&...>& v, meta::seq_t<S...>) {
            read_table_resize_cols_(n, m, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args) {
            read_table_resize_cols_i_(n, std::get<0>(v), v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            read_table_resize_(n, args...);
        }

        template<typename ... Args>
        void read_table_resize_(std::size_t n, impl::placeholder_t, Args& ... args) {
            read_table_resize_(n, args...);
        }

        template<typename T, typename ... Args, typename enable>
        void read_table_resize_(std::size_t n, T&, Args& ... args) {
            phypp_check(n <= 1, "cannot read multiple values into a scalar variable");
            read_table_resize_(n, args...);
        }

        template<typename T>
        void read_value_(line_splitter_t& spl, std::size_t i, std::size_t j, T& v) {
            std::string fallback;
            if (!spl.next_word(fallback)) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }

            if (!from_string(fallback, v)) {
                throw ascii::exception("cannot extract value '"+fallback+"' from file, wrong type for l."+
                    to_string(i+1)+":"+to_string(j+1)+" (expected '"+pretty_type(T())+"'):\n"+spl.line);
            }
        }

        void read_value_(line_splitter_t& spl, std::size_t i, std::size_t j, std::string& v) {
            if (!spl.next_word(v)) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }
        }

        inline void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j) {}

        template<typename T, typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, vec<1,T>& v, Args& ... args);
        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::tuple<U,VArgs&...> v, Args& ... args);
        template<typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, impl::placeholder_t, Args& ... args);
        template<typename T, typename ... Args, typename enable = typename std::enable_if<is_other<T>::value>::type>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, T& v, Args& ... args);

        template<typename T, typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, vec<1,T>& v, Args& ... args) {
            read_value_(spl, i, j, v.safe[i]);
            read_table_(spl, i, ++j, args...);
        }

        inline void read_table_cols_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k) {}

        template<typename T, typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args);
        template<typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k, impl::placeholder_t, VArgs&... args);

        template<typename T, typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args) {
            read_value_(spl, i, j, v.safe(i,k));
            read_table_cols_(spl, i, ++j, k, args...);
        }

        template<typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k, impl::placeholder_t, VArgs&... args) {
            if (!spl.skip_word()) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }

            read_table_cols_(spl, i, ++j, k, args...);
        }

        template<typename U, typename ... VArgs, std::size_t ... S>
        void read_table_cols_i_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::size_t k, std::tuple<U,VArgs&...>& v, meta::seq_t<S...>) {
            read_table_cols_(spl, i, j, k, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, std::tuple<U,VArgs&...> v, Args& ... args) {
            std::size_t n = std::get<0>(v);
            for (std::size_t k : range(n)) {
                read_table_cols_i_(spl, i, j, k, v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            }

            read_table_(spl, i, j, args...);
        }

        template<typename ... Args>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, impl::placeholder_t, Args& ... args) {
            if (!spl.skip_word()) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }

            read_table_(spl, i, ++j, args...);
        }

        template<typename T, typename ... Args, typename enable>
        void read_table_(line_splitter_t& spl, std::size_t i, std::size_t& j, T& v, Args& ... args) {
            read_value_(spl, i, j, v);
            read_table_(spl, i, ++j, args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void read_table(const std::string& name, const input_format::options& opts, Args&& ... args) {
        phypp_check(file::exists(name), "cannot open file '"+name+"'");

        try {
            impl::ascii_impl::line_splitter_t spl;
            spl.delim = opts.delim;
            spl.delim_single = opts.delim_single;

            std::ifstream file(name);

            // Count number of lines
            std::size_t n = 0; {
                while (!file.eof()) {
                    std::getline(file, spl.line);

                    auto p = spl.line.find_first_not_of(" \t");
                    if (p == spl.line.npos || (opts.auto_skip && spl.line[p] == opts.skip_pattern)) {
                        continue;
                    }

                    ++n;
                }

                n -= opts.skip_first;
            }

            // Resize all vectors
            impl::ascii_impl::read_table_resize_(n, args...);

            // Read data
            file.clear();
            file.seekg(0);

            uint_t i = 0;
            uint_t to_skip = opts.skip_first;
            while (!file.eof()) {
                spl.reset();
                std::getline(file, spl.line);

                auto p = spl.line.find_first_not_of(" \t");
                if (p == spl.line.npos || (opts.auto_skip && spl.line[p] == opts.skip_pattern)) {
                    continue;
                }

                if (to_skip > 0) {
                    --to_skip;
                    continue;
                }

                std::size_t j = 0;
                impl::ascii_impl::read_table_(spl, i, j, args...);
                ++i;
            }
        } catch (ascii::exception& e) {
            phypp_check(false, std::string(e.what())+" (reading "+name+")");
        }
    }

    namespace input_format {
        static const options standard = options{};
        static const options csv      = options{true, '#', 0, ",", true};
    }

    template<typename T, typename ... Args, typename enable = typename std::enable_if<
        !std::is_same<typename std::decay<T>::type, input_format::options>::value
    >::type>
    void read_table(const std::string& name, T&& t, Args&& ... args) {
        read_table(name, input_format::standard, std::forward<T>(t), std::forward<Args>(args)...);
    }
}

namespace impl {
    namespace ascii_impl {
        inline void write_table_check_size_(std::size_t n, std::size_t i) {}

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_check_size_(std::size_t& n, std::size_t i, const std::tuple<U,VArgs...>& v,
            const Args& ... args);

        template<typename F, typename ... Args,
            typename enable = typename std::enable_if<meta::is_format_tag<F>::value>::type>
        void write_table_check_size_(std::size_t& n, std::size_t i, const F& f,
            const Args& ... args);

        template<std::size_t Dim, typename Type, typename ... Args>
        void write_table_check_size_(std::size_t& n, std::size_t i, const vec<Dim,Type>& v,
            const Args& ... args) {

            if (n == 0) {
                n = v.dims[0];
            }

            if (v.dims[0] != n) {
                throw ascii::exception("incorrect dimension for column "+to_string(i)+" ("+
                    to_string(v.dims[0])+" vs "+to_string(n)+")");
            }

            write_table_check_size_(n, i+1, args...);
        }

        template<typename U, typename ... VArgs, std::size_t ... S>
        void write_table_check_size_tuple_(std::size_t& n, std::size_t i, const std::tuple<U,VArgs...>& v,
            meta::seq_t<S...>) {
            write_table_check_size_(n, i, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_check_size_(std::size_t& n, std::size_t i, const std::tuple<U,VArgs...>& v,
            const Args& ... args) {

            write_table_check_size_tuple_(n, i, v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            write_table_check_size_(n, i+sizeof...(VArgs), args...);
        }

        template<typename F, typename ... Args, typename enable>
        void write_table_check_size_(std::size_t& n, std::size_t i, const F& f,
            const Args& ... args) {

            write_table_check_size_(n, i, f.obj, args...);
        }

        inline void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j) {
            file << '\n';
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const std::tuple<U,VArgs...>& v, const Args& ... args);

        template<typename F, typename ... Args,
            typename enable = typename std::enable_if<meta::is_format_tag<F>::value>::type>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const F& v, const Args& ... args);

        template<std::size_t D, typename Type, typename ... Args>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const vec<D,Type>& v, const Args& ... args);

        template<typename Type, typename F, typename ... Args>
        void write_table_do_impl_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const vec<1,Type>& v, F&& fmt, const Args& ... args) {
                if (j == 0) {
                    file << std::string(sep.size(), ' ');
                } else {
                    file << sep;
                }

                std::string s = fmt(v.safe[i]);
                if (s.size() < cwidth) {
                    file << std::string(cwidth - s.size(), ' ');
                }

                file << s;

                write_table_do_(file, cwidth, sep, i, j+1, args...);
        }

        template<typename Type, typename F, typename ... Args>
        void write_table_do_impl_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const vec<2,Type>& v, F&& fmt, const Args& ... args) {

            for (uint_t k : range(v.dims[1])) {
                if (j == 0) {
                    file << std::string(sep.size(), ' ');
                } else {
                    file << sep;
                }

                std::string s = fmt(v(i,k));
                if (s.size() < cwidth) {
                    file << std::string(cwidth - s.size(), ' ');
                }

                file << s;
                ++j;
            }

            write_table_do_(file, cwidth, sep, i, j, args...);
        }

        template<std::size_t D, typename Type, typename ... Args>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const vec<D,Type>& v, const Args& ... args) {

            using DType = typename std::decay<decltype(v[0])>::type;
            write_table_do_impl_(file, cwidth, sep, i, j, v, [](const DType& t) {
                return to_string(t);
            }, args...);
        }

        template<typename F, typename ... Args, typename enable>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const F& v, const Args& ... args) {

            using DType = typename std::decay<decltype(v.obj[0])>::type;
            write_table_do_impl_(file, cwidth, sep, i, j, v.obj, [&](const DType& t) {
                return to_string(v.forward(t));
            }, args...);
        }

        inline void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j) {}

        template<typename Type, typename ... Args>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const vec<2,Type>& v, const Args& ... args);
        template<typename F, typename ... Args>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const F& v, const Args& ... args);
        template<typename U, typename ... VArgs, std::size_t ... S>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const std::tuple<U,VArgs...>& v, meta::seq_t<S...>);

        template<typename Type, typename F, typename ... Args>
        void write_table_do_tuple_impl_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const vec<2,Type>& v, F&& fmt,
            const Args& ... args) {

            if (j == 0) {
                file << std::string(sep.size(), ' ');
            } else {
                file << sep;
            }

            std::string s = fmt(v.safe(i,k));
            if (s.size() < cwidth) {
                file << std::string(cwidth - s.size(), ' ');
            }

            file << s;

            write_table_do_tuple_(file, cwidth, sep, i, k, j+1, args...);
        }

        template<typename Type, typename ... Args>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const vec<2,Type>& v, const Args& ... args) {

            using DType = typename std::decay<decltype(v[0])>::type;
            write_table_do_tuple_impl_(file, cwidth, sep, i, k, j, v, [](const DType& t) {
                return to_string(t);
            }, args...);
        }

        template<typename F, typename ... Args>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const F& v, const Args& ... args) {

            using DType = typename std::decay<decltype(v.obj[0])>::type;
            write_table_do_tuple_impl_(file, cwidth, sep, i, k, j, v.obj, [&](const DType& t) {
                return to_string(v.forward(t));
            }, args...);
        }

        template<typename U, typename ... VArgs, std::size_t ... S>
        void write_table_do_tuple_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t k, std::size_t j, const std::tuple<U,VArgs...>& v, meta::seq_t<S...>) {

            write_table_do_tuple_(file, cwidth, sep, i, k, j, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
            std::size_t i, std::size_t j, const std::tuple<U,VArgs...>& v, const Args& ... args) {

            uint_t m = std::get<0>(v); // number of column groups
            uint_t o = sizeof...(VArgs); // number of columns in a group
            for (uint_t k : range(m)) {
                write_table_do_tuple_(file, cwidth, sep, i, k, j+k*o, v,
                    typename meta::gen_seq<1, sizeof...(VArgs)>::type()
                );
            }

            write_table_do_(file, cwidth, sep, i, j+m*o, args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void write_table(const std::string& filename, std::size_t cwidth, const Args& ... args) {
        std::size_t n = 0, t = 0;
        impl::ascii_impl::write_table_check_size_(n, t, args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        try {
            for (std::size_t i : range(n)) {
                impl::ascii_impl::write_table_do_(file, cwidth, "", i, 0, args...);
            }
        } catch (ascii::exception& e) {
            phypp_check(false, std::string(e.what())+" (writing "+filename+")");
        }
    }

    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, const vec1s& nhdr,
        const Args& ... args) {

        std::size_t n = 0, t = 0;
        impl::ascii_impl::write_table_check_size_(n, t, args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        try {
            // Write header
            std::string hdr = collapse(align_right(nhdr, cwidth));
            if (hdr.size() > 1 && hdr[0] == ' ') hdr = hdr.substr(1);
            file << "#" << hdr << "\n#\n";

            for (std::size_t i : range(n)) {
                impl::ascii_impl::write_table_do_(file, cwidth, "", i, 0, args...);
            }
        } catch (ascii::exception& e) {
            phypp_check(false, std::string(e.what())+" (writing "+filename+")");
        }
    }
}

namespace impl {
    namespace ascii_impl {
        struct macroed_t {};

        inline std::string bake_macroed_name(std::string t) {
            // Just remove the leading object before the '.' operator
            auto p = t.find_first_of('.');
            if (p != t.npos) {
                t = t.substr(p+1);
            }
            return t;
        }

        inline std::string pop_macroed_name(std::string& s) {
            static const std::array<char,3> opening = {{'(', '[', '{'}};
            static const std::array<char,3> closing = {{')', ']', '}'}};

            std::string name;

            // Locate the next comma (',') that is not part of some nested function call
            // or constructor. Need a small parser...
            std::vector<uint_t> stack;
            uint_t p0 = s.find_first_not_of(" \t");
            uint_t pos = p0;
            for (; pos < s.size(); ++pos) {
                if (stack.empty() && s[pos] == ',') break;
                for (uint_t k : range(opening)) {
                    if (s[pos] == opening[k]) {
                        stack.push_back(k);
                        break;
                    } else if (s[pos] == closing[k]) {
                        phypp_check(stack.back() == k,
                            "library bug: error parsing macroed name list: '", s, "'");
                        stack.pop_back();
                        break;
                    }
                }
            }

            uint_t p1 = pos;
            if (pos != s.size()) {
                ++p1;
            }

            --pos;

            // Extract the trimmed name
            pos = s.find_last_not_of(" \t", pos);
            name = s.substr(p0, pos+1-p0);

            // Remove that element from the name list
            s = s.substr(p1);

            return name;
        }

        inline vec1s split_macroed_names(std::string s) {
            vec1s names;

            do {
                names.push_back(pop_macroed_name(s));
            } while (!s.empty());

            return names;
        }

        inline void write_table_hdr_fix_2d_(vec1s& vnames, uint_t i) {}

        template<typename T, typename ... Args>
        void write_table_hdr_fix_2d_(vec1s& vnames, uint_t i, const T&, const Args& ... args);
        template<typename T, typename ... Args>
        void write_table_hdr_fix_2d_(vec1s& vnames, uint_t i, const vec<2,T>& v, const Args& ... args);

        template<typename T, typename ... Args>
        void write_table_hdr_fix_2d_(vec1s& vnames, uint_t i, const T&, const Args& ... args) {
            write_table_hdr_fix_2d_(vnames, i+1, args...);
        }

        template<typename T, typename ... Args>
        void write_table_hdr_fix_2d_(vec1s& vnames, uint_t i, const vec<2,T>& v, const Args& ... args) {
            std::string base = vnames[i];

            vec1s ncols = base+"_"+to_stringa(uindgen(v.dims[1])+1);
            vnames.data.insert(vnames.data.begin()+i+1, ncols.begin()+1, ncols.end());
            vnames[i] = ncols[0];
            vnames.dims[0] += v.dims[1]-1;

            write_table_hdr_fix_2d_(vnames, i+v.dims[1], args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, impl::ascii_impl::macroed_t,
        const std::string& names, const Args& ... args) {
        vec1s vnames = impl::ascii_impl::split_macroed_names(names);
        for (auto& s : vnames) {
            s = impl::ascii_impl::bake_macroed_name(s);
        }

        impl::ascii_impl::write_table_hdr_fix_2d_(vnames, 0, args...);

        write_table_hdr(filename, cwidth, vnames, args...);
    }

    #define ftable(...) phypp::impl::ascii_impl::macroed_t(), #__VA_ARGS__, __VA_ARGS__

    template<typename ... Args>
    void write_table_csv(const std::string& filename, const Args& ... args) {
        std::size_t n = 0, t = 0;
        impl::ascii_impl::write_table_check_size_(n, t,args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        for (std::size_t i : range(n)) {
            impl::ascii_impl::write_table_do_(file, 0, ",", i, 0, args...);
        }
    }
}
}

#endif

