#ifndef PHYPP_FILE_HPP
#define PHYPP_FILE_HPP

// Emulate directory and file listing windows functions
// Note : Taken directly from Ogre3D
// http://www.ogre3d.org/

#include <dirent.h>
#include <unistd.h>
#include <fnmatch.h>
#include <cstring>
#include <cstdio>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include <tuple>
#include "phypp/vec.hpp"
#include "phypp/math.hpp"
#include "phypp/error.hpp"

namespace file {
    namespace impl {
        struct _finddata_t {
            char *name;
            int attrib;
            unsigned long size;
        };

        struct _find_search_t {
            char *pattern;
            char *curfn;
            char *directory;
            int dirlen;
            DIR *dirfd;
        };

        #define _A_NORMAL 0x00  /* Normalfile-Noread/writerestrictions */
        #define _A_RDONLY 0x01  /* Read only file */
        #define _A_HIDDEN 0x02  /* Hidden file */
        #define _A_SYSTEM 0x04  /* System file */
        #define _A_SUBDIR 0x10  /* Subdirectory */
        #define _A_ARCH   0x20  /* Archive file */

        int _findclose(long id);
        int _findnext(long id, _finddata_t *data);

        inline long _findfirst(const char *pattern, _finddata_t *data) {
            _find_search_t *fs = new _find_search_t;
            fs->curfn = NULL;
            fs->pattern = NULL;

            const char *mask = strrchr(pattern, '/');
            if (mask) {
                fs->dirlen = mask - pattern;
                mask++;
                fs->directory = static_cast<char*>(malloc(fs->dirlen + 1));
                memcpy(fs->directory, pattern, fs->dirlen);
                fs->directory[fs->dirlen] = 0;
            } else {
                mask = pattern;
                fs->directory = strdup(".");
                fs->dirlen = 1;
            }

            fs->dirfd = opendir(fs->directory);
            if (!fs->dirfd) {
                _findclose(reinterpret_cast<long>(fs));
                return -1;
            }

            if (strcmp(mask, "*.*") == 0) {
                mask += 2;
            }

            fs->pattern = strdup(mask);

            if (_findnext(reinterpret_cast<long>(fs), data) < 0) {
                _findclose(reinterpret_cast<long>(fs));
                return -1;
            }

            return reinterpret_cast<long>(fs);
        }

        inline int _findnext(long id, _finddata_t *data) {
            _find_search_t *fs = reinterpret_cast<_find_search_t*>(id);

            dirent *entry;
            for (;;) {
                if (!(entry = readdir(fs->dirfd))) {
                    return -1;
                }

                if (fnmatch(fs->pattern, entry->d_name, 0) == 0) {
                    break;
                }
            }

            if (fs->curfn) {
                free(fs->curfn);
            }

            data->name = fs->curfn = strdup(entry->d_name);

            size_t namelen = strlen(entry->d_name);
            char *xfn = new char[fs->dirlen + 1 + namelen + 1];
            sprintf(xfn, "%s/%s", fs->directory, entry->d_name);

            struct stat stat_buf;
            if (stat(xfn, &stat_buf)) {
                data->attrib = _A_NORMAL;
                data->size = 0;
            } else {
                if (S_ISDIR(stat_buf.st_mode)) {
                    data->attrib = _A_SUBDIR;
                } else {
                    data->attrib = _A_NORMAL;
                }

                data->size = stat_buf.st_size;
            }

            delete[] xfn;

            if (data->name [0] == '.') {
                data->attrib |= _A_HIDDEN;
            }

            return 0;
        }

        inline int _findclose(long id) {
            int ret;
            _find_search_t *fs = reinterpret_cast<_find_search_t*>(id);

            ret = fs->dirfd ? closedir(fs->dirfd) : 0;
            free(fs->pattern);
            free(fs->directory);
            if (fs->curfn)
                free(fs->curfn);
            delete fs;

            return ret;
        }
    }

    inline bool exists(const std::string& file) {
        if (file.empty()) {
            return false;
        }

        std::ifstream f(file.c_str());
        return f.is_open();
    }

    inline bool is_older(const std::string& file1, const std::string& file2) {
        struct stat st1, st2;
        if (::stat(file1.c_str(), &st1) != 0) return false;
        if (::stat(file2.c_str(), &st2) != 0) return false;
        return std::difftime(st1.st_ctime, st2.st_ctime) < 0.0;
    }

    inline bool is_absolute_path(const std::string& file) {
        auto pos = file.find_first_not_of(" \t");
        return pos != file.npos && file[pos] == '/';
    }

    inline bool copy(const std::string& file_from, const std::string& file_to) {
        std::ifstream src(file_from, std::ios::binary);
        if (!src.is_open()) return false;
        std::ofstream dst(file_to,   std::ios::binary);
        if (!dst.is_open()) return false;
        dst << src.rdbuf();
        return true;
    }

    inline bool remove(const std::string& file) {
        return !file::exists(file) || ::remove(file.c_str()) == 0;
    }

    inline std::string to_string(const std::string& file_name) {
        std::string   dst;
        std::ifstream src(file_name, std::ios::binary);
        if (src) {
            src.seekg(0, std::ios::end);
            dst.resize(src.tellg());
            src.seekg(0, std::ios::beg);
            src.read(&dst[0], dst.size());
            src.close();
        }

        return dst;
    }

    inline vec1s list_directories(const std::string& pattern = "*") {
        vec1s dlist;

        long handle, res;
        struct impl::_finddata_t tagData;

        handle = impl::_findfirst(pattern.c_str(), &tagData);
        res = 0;
        while (handle != -1 && res != -1) {
            if ((tagData.attrib & _A_HIDDEN) != 0) {
                res = impl::_findnext(handle, &tagData);
                continue;
            }

            if ((tagData.attrib & _A_SUBDIR) != 0) {
                std::string s = tagData.name;
                if (s != "." && s != "..") {
                    dlist.data.push_back(s);
                }
            }

            res = impl::_findnext(handle, &tagData);
        }

        if (handle != -1) {
            impl::_findclose(handle);
        }

        dlist.dims[0] = dlist.size();
        return dlist;
    }

    inline vec1s list_files(const std::string& pattern = "*") {
        vec1s flist;

        long handle, res;
        struct impl::_finddata_t tagData;

        handle = impl::_findfirst(pattern.c_str(), &tagData);
        res = 0;
        while (handle != -1 && res != -1) {
            if ((tagData.attrib & _A_HIDDEN) != 0) {
                res = impl::_findnext(handle, &tagData);
                continue;
            }

            if ((tagData.attrib & _A_SUBDIR) == 0) {
                flist.data.push_back(tagData.name);
            }

            res = impl::_findnext(handle, &tagData);
        }

        if (handle != -1) {
            impl::_findclose(handle);
        }

        flist.dims[0] = flist.size();
        return flist;
    }

    inline std::string directorize(const std::string& path) {
        std::string dir = trim(path);
        if (!dir.empty() && dir.back() != '/') {
            dir.push_back('/');
        }

        return dir;
    }

    // Same behavior as 'basename'
    inline std::string get_basename(std::string path) {
        auto pos = path.find_last_of('/');
        if (pos == path.npos) {
            return path;
        } else {
            auto lpos = path.find_last_not_of(' ');
            if (pos == lpos) {
                if (pos == 0) {
                    return "/";
                } else {
                    pos = path.find_last_of('/', pos-1);
                    if (pos == path.npos) {
                        return path.substr(lpos);
                    } else {
                        return path.substr(pos+1, lpos-pos-1);
                    }
                }
            } else {
                return path.substr(pos+1);
            }
        }
    }

    inline std::string remove_extension(std::string s) {
        auto p = s.find_last_of('.');
        if (p == s.npos) return s;
        return s.substr(0u, p);
    }

    inline std::string get_extension(std::string s) {
        auto p = s.find_last_of('.');
        if (p == s.npos) return "";
        return s.substr(p);
    }

    inline std::pair<std::string,std::string> split_extension(std::string s) {
        auto p = s.find_last_of('.');
        if (p == s.npos) return std::make_pair(s, std::string{});
        return std::make_pair(s.substr(0u, p), s.substr(p));
    }

    // Same behavior as 'dirname'
    inline std::string get_directory(const std::string& path) {
        auto pos = path.find_last_of('/');
        if (pos == path.npos) {
            return "./";
        } else if (pos == path.find_last_not_of(' ')) {
            if (pos == 0) {
                return "/";
            } else {
                pos = path.find_last_of('/', pos-1);
                if (pos == path.npos) {
                    return "./";
                } else {
                    return path.substr(0, pos+1);
                }
            }
        } else {
            return path.substr(0, pos+1);
        }
    }

    inline bool mkdir(const std::string& path) {
        if (path.empty()) return true;
        vec1s dirs = split(path, "/");
        std::string tmp;
        for (auto& d : dirs) {
            if (d.empty()) continue;
            if (!tmp.empty() || start_with(path, "/")) tmp += "/";
            tmp += d;
            bool res = (::mkdir(tmp.c_str(), 0775) == 0) || (errno == EEXIST);
            if (!res) return false;
        }
        return true;
    }


#define VECTORIZE(name) \
    template<std::size_t Dim, typename Type, typename ... Args, \
        typename enable = typename std::enable_if< \
            std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type> \
    auto name(const vec<Dim,Type>& v, const Args& ... args) -> \
        vec<Dim,decltype(name(v[0], args...))> { \
        using ntype = decltype(name(v[0], args...)); \
        vec<Dim,ntype> r(v.dims); \
        for (uint_t i : range(v)) { \
            r.safe[i] = name(v.safe[i], args...); \
        } \
        return r; \
    }

    VECTORIZE(directorize)
    VECTORIZE(is_absolute_path)
    VECTORIZE(get_basename)
    VECTORIZE(get_directory)
    VECTORIZE(remove_extension)
    VECTORIZE(get_extension)
    VECTORIZE(split_extension)
    VECTORIZE(mkdir)
    VECTORIZE(exists)
    VECTORIZE(is_older)

#undef VECTORIZE

    inline uint_t find_skip(const std::string& name) {
        phypp_check(file::exists(name), "cannot open file '"+name+"'");

        std::string line;

        std::size_t n = 0;
        std::ifstream file(name.c_str());
        while (std::getline(file, line)) {
            auto p = line.find_first_not_of(" \t");
            if (p != line.npos && line[p] != '#') {
                break;
            }

            ++n;
        }

        return n;
    }

    template<typename T, typename ... Args>
    auto columns(std::size_t n, T& t, Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }

    inline void read_table_resize_(std::size_t n) {}

    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args);
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args);
    template<typename ... Args>
    void read_table_resize_(std::size_t n, placeholder_t, Args& ... args);

    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args) {
        v.resize(n);
        read_table_resize_(n, args...);
    }

    inline void read_table_resize_cols_(std::size_t n, std::size_t m) {}

    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, placeholder_t, VArgs&... args);

    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args) {
        v.resize(n, m);
        read_table_resize_cols_(n, m, args...);
    }

    template<typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, placeholder_t, VArgs&... args) {
        read_table_resize_cols_(n, m, args...);
    }

    template<typename U, typename ... VArgs, std::size_t ... S>
    void read_table_resize_cols_i_(std::size_t n, std::size_t m, std::tuple<U,VArgs&...>& v, seq_t<S...>) {
        read_table_resize_cols_(n, m, std::get<S>(v)...);
    }

    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args) {
        read_table_resize_cols_i_(n, std::get<0>(v), v, typename gen_seq<1, sizeof...(VArgs)>::type());
        read_table_resize_(n, args...);
    }

    template<typename ... Args>
    void read_table_resize_(std::size_t n, placeholder_t, Args& ... args) {
        read_table_resize_(n, args...);
    }

    template<typename I, typename T>
    bool read_value_(I& in, T& v, std::string& fallback) {
        auto pos = in.tellg();
        in >> fallback;

        std::istringstream ss(fallback);
        ss >> v;
        if (ss.fail() || !ss.eof()) {
            in.clear();
            in.seekg(pos);
            return false;
        }

        return true;
    }

    template<typename I, typename T>
    bool read_value_float_(I& in, T& v, std::string& fallback) {
        auto pos = in.tellg();
        in >> fallback;

        std::istringstream ss(fallback);
        ss >> v;
        if (ss.fail()) {
            std::string s = trim(toupper(fallback));
            if (s == "NAN" || s == "+NAN" || s == "-NAN") {
                v = dnan;
                return true;
            } else if (s == "+INF" || s == "INF+" || s == "INF") {
                v = dinf;
                return true;
            } else if (s == "-INF" || s == "INF-") {
                v = -dinf;
                return true;
            } else if (s == "NULL") {
                v = dnan;
                return true;
            }

            in.clear();
            in.seekg(pos);
            return false;
        } else if (!ss.eof()) {
            in.clear();
            in.seekg(pos);
            return false;
        }

        return true;
    }

    template<typename I>
    bool read_value_(I& in, float& v, std::string& fallback) {
        return read_value_float_(in, v, fallback);
    }

    template<typename I>
    bool read_value_(I& in, double& v, std::string& fallback) {
        return read_value_float_(in, v, fallback);
    }

    inline void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j) {}

    template<typename T, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, vec<1,T>& v, Args& ... args);
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, std::tuple<U,VArgs&...> v, Args& ... args);
    template<typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, placeholder_t, Args& ... args);

    template<typename T, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, vec<1,T>& v, Args& ... args) {
        std::string fb;
        if (!read_value_(fs, v[i], fb)) {
            if (fb.empty()) {
                phypp_check(false, "cannot extract value from file, too few columns on line l."+strn(i+1));
            } else {
                phypp_check(false, "cannot extract value '", fb, "' from file, wrong type for l."+
                    strn(i+1)+":"+strn(j+1)+" (expected '"+pretty_type(T())+"'):\n"+fs.str());
            }
        }
        read_table_(fs, i, ++j, args...);
    }

    inline void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k) {}

    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, placeholder_t, VArgs&... args);

    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args) {
        std::string fb;
        if (!read_value_(fs, v(i,k), fb)) {
            if (fb.empty()) {
                phypp_check(false, "cannot extract value from file, too few columns on line l."+strn(i+1));
            } else {
                phypp_check(false, "cannot extract value '", fb, "' from file, wrong type for l."+
                    strn(i+1)+":"+strn(j+1)+" (expected '"+pretty_type(T())+"'):\n"+fs.str());
            }
        }
        read_table_cols_(fs, i, ++j, k, args...);
    }

    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, placeholder_t, VArgs&... args) {
        std::string s;
        if (!(fs >> s)) {
            phypp_check(!fs.eof(), "cannot extract value from file, "
                "too few columns on line l."+strn(i+1));
        }
        read_table_cols_(fs, i, ++j, k, args...);
    }

    template<typename U, typename ... VArgs, std::size_t ... S>
    void read_table_cols_i_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, std::tuple<U,VArgs&...>& v, seq_t<S...>) {
        read_table_cols_(fs, i, j, k, std::get<S>(v)...);
    }

    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, std::tuple<U,VArgs&...> v, Args& ... args) {
        std::size_t n = std::get<0>(v);
        for (std::size_t k = 0; k < n; ++k) {
            read_table_cols_i_(fs, i, j, k, v, typename gen_seq<1, sizeof...(VArgs)>::type());
        }

        read_table_(fs, i, j, args...);
    }

    template<typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j, placeholder_t, Args& ... args) {
        phypp_check(!fs.eof(), "cannot extract value at l."+strn(i+1)+":"+strn(j+1)+" from file, "
            "too few columns on line l."+strn(i+1));
        std::string s;
        fs >> s;
        read_table_(fs, i, ++j, args...);
    }

    template<typename ... Args>
    void read_table(const std::string& name, std::size_t skip, Args&& ... args) {
        phypp_check(file::exists(name), "cannot open file '"+name+"'");

        std::string line;

        std::size_t n = 0;
        std::size_t rpos = 0; {
            std::ifstream file(name.c_str());
            while (!file.eof()) {
                std::getline(file, line);

                if (n < skip) {
                    ++n;
                    if (n == skip) {
                        rpos = file.tellg();
                    }
                } else {
                    if (line.find_first_not_of(" \t") != line.npos) {
                        ++n;
                    }
                }
            }

            if (skip > n) return;

            n -= skip;
        }

        std::ifstream file(name.c_str());
        file.seekg(rpos);

        read_table_resize_(n, args...);

        for (std::size_t i = 0; i < n; ++i) {
            do {
                std::getline(file, line);
            } while (line.find_first_not_of(" \t") == line.npos);

            std::istringstream fs(line);
            std::size_t j = 0;
            read_table_(fs, i, j, args...);
        }
    }

    inline void write_table_check_size_(std::size_t n, std::size_t i) {}

    template<std::size_t Dim, typename Type, typename ... Args>
    void write_table_check_size_(std::size_t& n, std::size_t i, const vec<Dim,Type>& v,
        const Args& ... args) {

        if (n == 0) {
            n = v.dims[0];
        }

        phypp_check(v.dims[0] == n, "incorrect dimension for column "+strn(i)+" ("+
            strn(v.dims[0])+" vs "+strn(n)+")");

        write_table_check_size_(n, i+1, args...);
    }

    inline void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
        std::size_t i, std::size_t j) {
        file << '\n';
    }

    template<typename Type, typename ... Args>
    void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
        std::size_t i, std::size_t j, const vec<2,Type>& v, const Args& ... args);

    template<typename Type, typename ... Args>
    void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
        std::size_t i, std::size_t j, const vec<1,Type>& v, const Args& ... args) {

        if (j == 0) {
            file << std::string(sep.size(), ' ');
        } else {
            file << sep;
        }

        std::string s = strn(v[i]);
        if (s.size() < cwidth) {
            file << std::string(cwidth - s.size(), ' ');
        }

        file << s;

        write_table_do_(file, cwidth, sep, i, j+1, args...);
    }

    template<typename Type, typename ... Args>
    void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
        std::size_t i, std::size_t j, const vec<2,Type>& v, const Args& ... args) {

        for (uint_t k = 0; k < v.dims[1]; ++k) {
            if (j == 0) {
                file << std::string(sep.size(), ' ');
            } else {
                file << sep;
            }

            std::string s = strn(v(i,k));
            if (s.size() < cwidth) {
                file << std::string(cwidth - s.size(), ' ');
            }

            file << s;
            ++j;
        }

        write_table_do_(file, cwidth, sep, i, j, args...);
    }

    template<typename ... Args>
    void write_table(const std::string& filename, std::size_t cwidth, const Args& ... args) {
        std::size_t n = 0, t = 0;
        write_table_check_size_(n, t, args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, "", i, 0, args...);
        }
    }

    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, const vec1s& nhdr,
        const Args& ... args) {

        std::size_t n = 0, t = 0;
        write_table_check_size_(n, t, args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        // Write header
        std::string hdr = collapse(align_right(nhdr, cwidth));
        if (hdr.size() > 1 && hdr[0] == ' ') hdr = hdr.substr(1);
        file << "#" << hdr << "\n#\n";

        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, "", i, 0, args...);
        }
    }

    struct macroed_t {};
    #define ftable(...) file::macroed_t(), #__VA_ARGS__, __VA_ARGS__

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

        vec1s ncols = base+"_"+strna(uindgen(v.dims[1])+1);
        vnames.data.insert(vnames.data.begin()+i+1, ncols.begin()+1, ncols.end());
        vnames[i] = ncols[0];
        vnames.dims[0] += v.dims[1]-1;

        write_table_hdr_fix_2d_(vnames, i+v.dims[1], args...);
    }

    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, macroed_t,
        const std::string& names, const Args& ... args) {
        vec1s vnames = file::split_macroed_names(names);
        for (auto& s : vnames) {
            s = file::bake_macroed_name(s);
        }

        write_table_hdr_fix_2d_(vnames, 0, args...);

        write_table_hdr(filename, cwidth, vnames, args...);
    }

    template<typename ... Args>
    void write_table_csv(const std::string& filename, std::size_t cwidth, const Args& ... args) {
        std::size_t n = 0, t = 0;
        write_table_check_size_(n, t,args...);

        std::ofstream file(filename);
        phypp_check(file.is_open(), "could not open file "+filename+" to write data");

        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, ",", i, 0, args...);
        }
    }
}

inline bool fork(const std::string& cmd) {
    return system((cmd+" &").c_str()) == 0;
}

inline bool spawn(const std::string& cmd) {
    return system(cmd.c_str()) == 0;
}

#endif

