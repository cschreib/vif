#ifndef VIF_IO_ASCII_HPP
#define VIF_IO_ASCII_HPP

#include <fstream>
#include <tuple>
#include <stdexcept>
#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/math/base.hpp"

namespace vif {
namespace ascii {
    // ASCII exception type
    struct exception : std::runtime_error {
        exception(const std::string& w) : std::runtime_error(w) {}
    };

    // Control table layout for read_table()
    struct input_format {
        bool auto_skip    = true;
        std::string skip_pattern = "#";
        uint_t skip_first = 0;
        std::string delim = " \t";
        bool delim_single = false;

        input_format() = default;
        input_format(bool sk, const std::string sp, uint_t sf, const std::string& d, bool ds) :
            auto_skip(sk), skip_pattern(sp), skip_first(sf), delim(d), delim_single(ds) {}

        static input_format standard() {
            return input_format{};
        }

        static input_format csv() {
            return input_format{true, "#", 0, ",", true};
        }
    };

    // Control table layout for write_table()
    struct output_format {
        bool auto_width   = true;
        uint_t min_width  = 0;
        std::string delim = " ";
        std::string header_chars = "# ";
        vec1s header;

        output_format() = default;
        output_format(bool aw, uint_t mw, const std::string& d) :
            auto_width(aw), min_width(mw), delim(d) {}

        static output_format standard() {
            return output_format{};
        }
        static output_format csv() {
            return output_format{false, 0, ","};
        }
    };

    // Helper function to merge multiple columns into one vector
    template<typename T, typename ... Args>
    auto columns(uint_t n, T& t, Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }

    template<typename T, typename ... Args>
    auto columns(uint_t n, const T& t, const Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }
}

namespace impl {
    namespace ascii_impl {
        using placeholder_t = vif::impl::placeholder_t;

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

        inline void read_table_resize_(uint_t n) {}

        template<typename T, typename ... Args>
        void read_table_resize_(uint_t n, vec<1,T>& v, Args& ... args);
        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_resize_(uint_t n, std::tuple<U,VArgs&...> v, Args& ... args);
        template<typename ... Args>
        void read_table_resize_(uint_t n, impl::placeholder_t, Args& ... args);
        template<typename T, typename ... Args, typename enable = typename std::enable_if<is_other<T>::value>::type>
        void read_table_resize_(uint_t n, T&, Args& ... args);

        template<typename T, typename ... Args>
        void read_table_resize_(uint_t n, vec<1,T>& v, Args& ... args) {
            v.resize(n);
            read_table_resize_(n, args...);
        }

        inline void read_table_resize_cols_(uint_t n, uint_t m) {}

        template<typename Type, typename ... VArgs>
        void read_table_resize_cols_(uint_t n, uint_t m, vec<2,Type>& v, VArgs&... args);
        template<typename ... VArgs>
        void read_table_resize_cols_(uint_t n, uint_t m, impl::placeholder_t, VArgs&... args);

        template<typename Type, typename ... VArgs>
        void read_table_resize_cols_(uint_t n, uint_t m, vec<2,Type>& v, VArgs&... args) {
            v.resize(n, m);
            read_table_resize_cols_(n, m, args...);
        }

        template<typename ... VArgs>
        void read_table_resize_cols_(uint_t n, uint_t m, impl::placeholder_t, VArgs&... args) {
            read_table_resize_cols_(n, m, args...);
        }

        template<typename U, typename ... VArgs, uint_t ... S>
        void read_table_resize_cols_i_(uint_t n, uint_t m, std::tuple<U,VArgs&...>& v, meta::seq_t<S...>) {
            read_table_resize_cols_(n, m, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_resize_(uint_t n, std::tuple<U,VArgs&...> v, Args& ... args) {
            read_table_resize_cols_i_(n, std::get<0>(v), v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            read_table_resize_(n, args...);
        }

        template<typename ... Args>
        void read_table_resize_(uint_t n, impl::placeholder_t, Args& ... args) {
            read_table_resize_(n, args...);
        }

        template<typename T, typename ... Args, typename enable>
        void read_table_resize_(uint_t n, T&, Args& ... args) {
            vif_check(n <= 1, "cannot read multiple values into a scalar variable");
            read_table_resize_(n, args...);
        }

        template<typename T>
        void read_value_(line_splitter_t& spl, uint_t i, uint_t j, T& v) {
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

        inline void read_value_(line_splitter_t& spl, uint_t i, uint_t j, std::string& v) {
            if (!spl.next_word(v)) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }
        }

        inline void read_table_(line_splitter_t& spl, uint_t i, uint_t& j) {}

        template<typename T, typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, vec<1,T>& v, Args& ... args);
        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, std::tuple<U,VArgs&...> v, Args& ... args);
        template<typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, impl::placeholder_t, Args& ... args);
        template<typename T, typename ... Args, typename enable = typename std::enable_if<is_other<T>::value>::type>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, T& v, Args& ... args);

        template<typename T, typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, vec<1,T>& v, Args& ... args) {
            read_value_(spl, i, j, v.safe[i]);
            read_table_(spl, i, ++j, args...);
        }

        inline void read_table_cols_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k) {}

        template<typename T, typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k, vec<2,T>& v, VArgs&... args);
        template<typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k, impl::placeholder_t, VArgs&... args);

        template<typename T, typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k, vec<2,T>& v, VArgs&... args) {
            read_value_(spl, i, j, v.safe(i,k));
            read_table_cols_(spl, i, ++j, k, args...);
        }

        template<typename ... VArgs>
        void read_table_cols_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k, impl::placeholder_t, VArgs&... args) {
            if (!spl.skip_word()) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }

            read_table_cols_(spl, i, ++j, k, args...);
        }

        template<typename U, typename ... VArgs, uint_t ... S>
        void read_table_cols_i_(line_splitter_t& spl, uint_t i, uint_t& j, uint_t k, std::tuple<U,VArgs&...>& v, meta::seq_t<S...>) {
            read_table_cols_(spl, i, j, k, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, std::tuple<U,VArgs&...> v, Args& ... args) {
            uint_t n = std::get<0>(v);
            for (uint_t k : range(n)) {
                read_table_cols_i_(spl, i, j, k, v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            }

            read_table_(spl, i, j, args...);
        }

        template<typename ... Args>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, impl::placeholder_t, Args& ... args) {
            if (!spl.skip_word()) {
                throw ascii::exception("cannot extract value from file, too few columns on line l."+
                    to_string(i+1));
            }

            read_table_(spl, i, ++j, args...);
        }

        template<typename T, typename ... Args, typename enable>
        void read_table_(line_splitter_t& spl, uint_t i, uint_t& j, T& v, Args& ... args) {
            read_value_(spl, i, j, v);
            read_table_(spl, i, ++j, args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void read_table(const std::string& name, const input_format& opts, Args&& ... args) {
        vif_check(file::exists(name), "cannot open file '"+name+"'");

        try {
            impl::ascii_impl::line_splitter_t spl;
            spl.delim = opts.delim;
            spl.delim_single = opts.delim_single;

            std::ifstream file(name);

            // Count number of lines
            uint_t n = 0; {
                while (!file.eof()) {
                    std::getline(file, spl.line);

                    auto p = spl.line.find_first_not_of(" \t");
                    if (p == spl.line.npos || (opts.auto_skip && spl.line.find(opts.skip_pattern) == p)) {
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
                if (p == spl.line.npos || (opts.auto_skip && spl.line.find(opts.skip_pattern) == p)) {
                    continue;
                }

                if (to_skip > 0) {
                    --to_skip;
                    continue;
                }

                uint_t j = 0;
                impl::ascii_impl::read_table_(spl, i, j, args...);
                ++i;
            }
        } catch (ascii::exception& e) {
            vif_check(false, std::string(e.what())+" (reading "+name+")");
        }
    }

    template<typename T, typename ... Args, typename enable = typename std::enable_if<
        !std::is_same<typename std::decay<T>::type, input_format>::value
    >::type>
    void read_table(const std::string& name, T&& t, Args&& ... args) {
        read_table(name, input_format::standard(), std::forward<T>(t), std::forward<Args>(args)...);
    }
}

namespace impl {
    namespace ascii_impl {
        struct file_writer {
            std::ofstream out;
            std::string delim;

            uint_t j = 0;
            vec1u cwidth;

            vec1s header;
            std::string header_chars;

            void end_line() {
                out << '\n';
                j = 0;
            }

            void write(std::string s) {
                if (j != 0) {
                    out << delim;
                }

                out << align_right(std::move(s), cwidth[j]);

                ++j;
            }

            void write_header() {
                if (header.empty()) return;

                out << header_chars;
                for (uint_t i : range(header)) {
                    uint_t hw = cwidth[i];
                    if (i == 0 && hw >= header_chars.size()) {
                        hw -= header_chars.size();
                    }

                    if (i != 0) {
                        out << delim;
                    }

                    out << align_right(header[i], hw);
                }

                out << "\n";
            }
        };

        struct cache_writer {
            vec2s data;
            uint_t k = 0;

            uint_t j = 0;
            vec1u cwidth;

            void end_line() {
                j = 0;
            }

            void write(std::string s) {
                data.safe[k] = align_right(std::move(s), cwidth[j]);
                ++k;
                ++j;
            }
        };

        inline void write_table_check_size_(uint_t& r, uint_t& c) {}

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_check_size_(uint_t& r, uint_t& c, const std::tuple<U,VArgs...>& v,
            const Args& ... args);

        template<typename F, typename ... Args,
            typename enable = typename std::enable_if<meta::is_format_tag<F>::value>::type>
        void write_table_check_size_(uint_t& r, uint_t& c, const F& f,
            const Args& ... args);

        template<uint_t Dim, typename Type, typename ... Args>
        void write_table_check_size_(uint_t& r, uint_t& c, const vec<Dim,Type>& v,
            const Args& ... args) {

            static_assert(Dim <= 2, "cannot write vectors of more than two dimensions");

            if (r == 0) {
                r = v.dims[0];
            }

            if (v.dims[0] != r) {
                throw ascii::exception("incorrect dimension for column "+to_string(c)+" ("+
                    to_string(v.dims[0])+" vs "+to_string(r)+")");
            }

            if (Dim == 2) {
                c += v.dims[1];
            } else {
                ++c;
            }

            write_table_check_size_(r, c, args...);
        }

        template<typename U, typename ... VArgs, uint_t ... S>
        void write_table_check_size_tuple_(uint_t& r, uint_t& c, const std::tuple<U,VArgs...>& v,
            meta::seq_t<S...>) {

            write_table_check_size_(r, c, std::get<S>(v)...);
        }

        template<typename U, typename ... VArgs, typename ... Args>
        void write_table_check_size_(uint_t& r, uint_t& c, const std::tuple<U,VArgs...>& v,
            const Args& ... args) {

            write_table_check_size_tuple_(r, c, v, typename meta::gen_seq<1, sizeof...(VArgs)>::type());
            write_table_check_size_(r, c, args...);
        }

        template<typename F, typename ... Args, typename enable>
        void write_table_check_size_(uint_t& r, uint_t& c, const F& f, const Args& ... args) {
            write_table_check_size_(r, c, f.obj, args...);
        }

        template<typename O>
        void write_table_do_(O&& out, uint_t) {
            out.end_line();
        }

        template<typename O, typename U, typename ... VArgs, typename ... Args>
        void write_table_do_(O&& out, uint_t i, const std::tuple<U,VArgs...>& v,
            const Args& ... args);

        template<typename O, typename F, typename ... Args,
            typename enable = typename std::enable_if<meta::is_format_tag<F>::value>::type>
        void write_table_do_(O&& out, uint_t i, const F& v, const Args& ... args);

        template<typename O, uint_t D, typename Type, typename ... Args>
        void write_table_do_(O&& out, uint_t i, const vec<D,Type>& v, const Args& ... args);

        template<typename O, typename Type, typename F, typename ... Args>
        void write_table_do_impl_(O&& out, uint_t i, const vec<1,Type>& v, F&& fmt,
            const Args& ... args) {

            out.write(fmt(v.safe[i]));

            write_table_do_(out, i, args...);
        }

        template<typename O, typename Type, typename F, typename ... Args>
        void write_table_do_impl_(O&& out, uint_t i, const vec<2,Type>& v, F&& fmt,
            const Args& ... args) {

            for (uint_t k : range(v.dims[1])) {
                out.write(fmt(v(i,k)));
            }

            write_table_do_(out, i, args...);
        }

        template<typename O, uint_t D, typename Type, typename ... Args>
        void write_table_do_(O&& out, uint_t i, const vec<D,Type>& v, const Args& ... args) {
            using DType = typename std::decay<decltype(v[0])>::type;
            write_table_do_impl_(out, i, v, [](const DType& t) {
                return to_string(t);
            }, args...);
        }

        template<typename O, typename F, typename ... Args, typename enable>
        void write_table_do_(O&& out, uint_t i, const F& v, const Args& ... args) {
            using DType = typename std::decay<decltype(v.obj[0])>::type;
            write_table_do_impl_(out, i, v.obj, [&](const DType& t) {
                return to_string(v.forward(t));
            }, args...);
        }

        template<typename O>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k) {}

        template<typename O, typename Type, typename ... Args>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const vec<2,Type>& v,
            const Args& ... args);
        template<typename O, typename F, typename ... Args>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const F& v, const Args& ... args);
        template<typename O, typename U, typename ... VArgs, uint_t ... S>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const std::tuple<U,VArgs...>& v,
            meta::seq_t<S...>);

        template<typename O, typename Type, typename F, typename ... Args>
        void write_table_do_tuple_impl_(O&& out, uint_t i, uint_t k, const vec<2,Type>& v, F&& fmt,
            const Args& ... args) {

            out.write(fmt(v.safe(i,k)));
            write_table_do_tuple_(out, i, k, args...);
        }

        template<typename O, typename Type, typename ... Args>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const vec<2,Type>& v,
            const Args& ... args) {

            using DType = typename std::decay<decltype(v[0])>::type;
            write_table_do_tuple_impl_(out, i, k, v, [](const DType& t) {
                return to_string(t);
            }, args...);
        }

        template<typename O, typename F, typename ... Args>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const F& v, const Args& ... args) {
            using DType = typename std::decay<decltype(v.obj[0])>::type;
            write_table_do_tuple_impl_(out, i, k, v.obj, [&](const DType& t) {
                return to_string(v.forward(t));
            }, args...);
        }

        template<typename O, typename U, typename ... VArgs, uint_t ... S>
        void write_table_do_tuple_(O&& out, uint_t i, uint_t k, const std::tuple<U,VArgs...>& v,
            meta::seq_t<S...>) {

            write_table_do_tuple_(out, i, k, std::get<S>(v)...);
        }

        template<typename O, typename U, typename ... VArgs, typename ... Args>
        void write_table_do_(O&& out, uint_t i, const std::tuple<U,VArgs...>& v,
            const Args& ... args) {

            uint_t m = std::get<0>(v); // number of column groups
            for (uint_t k : range(m)) {
                write_table_do_tuple_(out, i, k, v,
                    typename meta::gen_seq<1, sizeof...(VArgs)>::type()
                );
            }

            write_table_do_(out, i, args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void write_table(const std::string& filename, const output_format& opts,
        const Args& ... args) {

        impl::ascii_impl::file_writer file;
        file.out.open(filename);
        file.delim = opts.delim;
        file.header = opts.header;
        file.header_chars = opts.header_chars;

        // Check we can write to the file
        vif_check(file.out.is_open(), "could not open file "+filename+" to write data");

        // Compute number of rows and check that all vectors have the same size
        uint_t r = 0, c = 0;
        impl::ascii_impl::write_table_check_size_(r, c, args...);

        if (!opts.header.empty()) {
            vif_check(opts.header.size() == c, "mismatch between dimensions of header and number "
                "of columns to write (", opts.header.size(), " vs. ", c, ")");
        }

        file.cwidth = replicate(opts.min_width, c);

        try {
            if (opts.auto_width) {
                // Pre-serialize data to compute column width
                impl::ascii_impl::cache_writer cache;
                cache.data.resize(r, c);
                cache.cwidth = file.cwidth;
                for (uint_t i : range(r)) {
                    impl::ascii_impl::write_table_do_(cache, i, args...);
                }

                // Compute width
                for (uint_t j : range(c)) {
                    for (auto& s : cache.data.stride(_,j)) {
                        cache.cwidth.safe[j] = std::max(cache.cwidth.safe[j], s.size());
                    }
                }

                // Increase width if header is larger
                if (!opts.header.empty()) {
                    for (uint_t j : range(c)) {
                        uint_t hs = opts.header[j].size();
                        if (j == 0) {
                            hs += opts.header_chars.size();
                        }
                        cache.cwidth.safe[j] = std::max(cache.cwidth.safe[j], hs);
                    }
                }

                // Go back to writing to file
                file.cwidth = cache.cwidth;

                // Write header
                file.write_header();

                // Write pre-serialized data
                for (uint_t i : range(r)) {
                    for (uint_t j : range(c)) {
                        file.write(cache.data.safe(i,j));
                    }

                    file.end_line();
                }
            } else {
                // Write directly to file

                // Write header
                file.write_header();

                // Write data
                for (uint_t i : range(r)) {
                    impl::ascii_impl::write_table_do_(file, i, args...);
                }
            }
        } catch (ascii::exception& e) {
            vif_check(false, std::string(e.what())+" (writing "+filename+")");
        }
    }

    template<typename T, typename ... Args, typename enable = typename std::enable_if<
        !std::is_same<typename std::decay<T>::type, output_format>::value
    >::type>
    void write_table(const std::string& filename, const T& t,
        const Args& ... args) {
        write_table(filename, output_format::standard(), t, args...);
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
                        vif_check(stack.back() == k,
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

            vec1s ncols = base+"_"+to_string_vector(indgen<uint_t>(v.dims[1])+1);
            vnames.data.insert(vnames.data.begin()+i+1, ncols.begin()+1, ncols.end());
            vnames[i] = ncols[0];
            vnames.dims[0] += v.dims[1]-1;

            write_table_hdr_fix_2d_(vnames, i+v.dims[1], args...);
        }
    }
}

namespace ascii {
    template<typename ... Args>
    void write_table(const std::string& filename, output_format opts, impl::ascii_impl::macroed_t,
        const char* names, const Args& ... args) {

        opts.header = impl::ascii_impl::split_macroed_names(names);
        for (auto& s : opts.header) {
            s = impl::ascii_impl::bake_macroed_name(s);
        }

        impl::ascii_impl::write_table_hdr_fix_2d_(opts.header, 0, args...);

        write_table(filename, opts, args...);
    }

    template<typename ... Args>
    void write_table(const std::string& filename, impl::ascii_impl::macroed_t,
        const char* names, const Args& ... args) {
        write_table(filename, output_format::standard(), impl::ascii_impl::macroed_t{}, names, args...);
    }

    #define ftable(...) vif::impl::ascii_impl::macroed_t(), #__VA_ARGS__, __VA_ARGS__
}
}

#endif
