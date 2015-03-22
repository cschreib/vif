#ifndef FILE_HPP
#define FILE_HPP

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
        int _findnext(long id, struct _finddata_t *data);

        long _findfirst(const char *pattern, struct _finddata_t *data) {
            _find_search_t *fs = new _find_search_t;
            fs->curfn = NULL;
            fs->pattern = NULL;

            const char *mask = strrchr(pattern, '/');
            if (mask) {
                fs->dirlen = mask - pattern;
                mask++;
                fs->directory = (char *)malloc(fs->dirlen + 1);
                memcpy(fs->directory, pattern, fs->dirlen);
                fs->directory[fs->dirlen] = 0;
            } else {
                mask = pattern;
                fs->directory = strdup(".");
                fs->dirlen = 1;
            }

            fs->dirfd = opendir(fs->directory);
            if (!fs->dirfd) {
                _findclose((long)fs);
                return -1;
            }

            if (strcmp(mask, "*.*") == 0) {
                mask += 2;
            }

            fs->pattern = strdup(mask);

            if (_findnext((long)fs, data) < 0) {
                _findclose((long)fs);
                return -1;
            }

            return (long)fs;
        }

        int _findnext(long id, struct _finddata_t *data) {
            _find_search_t *fs = (_find_search_t*)id;

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

        int _findclose(long id) {
            int ret;
            _find_search_t *fs = (_find_search_t *)id;

            ret = fs->dirfd ? closedir(fs->dirfd) : 0;
            free(fs->pattern);
            free(fs->directory);
            if (fs->curfn)
                free(fs->curfn);
            delete fs;

            return ret;
        }
    }

    bool exists(const std::string& file) {
        if (file.empty()) {
            return false;
        }

        std::ifstream f(file.c_str());
        return f.is_open();
    }

    bool is_older(const std::string& file1, const std::string& file2) {
        struct stat st1, st2;
        if (::stat(file1.c_str(), &st1) != 0) return false;
        if (::stat(file2.c_str(), &st2) != 0) return false;
        return std::difftime(st1.st_ctime, st2.st_ctime) < 0.0;
    }

    void copy(const std::string& file_from, const std::string& file_to) {
        std::ifstream src(file_from, std::ios::binary);
        std::ofstream dst(file_to,   std::ios::binary);
        dst << src.rdbuf();
    }

    vec1s list_directories(const std::string& path = "") {
        vec1s dlist;

        long handle, res;
        struct impl::_finddata_t tagData;

        std::string pattern;
        if (path.empty()) {
            pattern = "*";
        } else {
            pattern = path + "/*";
        }

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

    vec1s list_files(const std::string& pattern = "*") {
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

    std::string directorize(const std::string& path) {
        std::string dir = trim(path);
        if (!dir.empty() && dir.back() != '/') {
            dir.push_back('/');
        }

        return dir;
    }

    // Same behavior as 'basename'
    std::string get_basename(std::string path) {
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

    // Same behavior as 'dirname'
    std::string get_directory(const std::string& path) {
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

    bool mkdir(const std::string& path) {
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
        vec<Dim,ntype> r = arr<ntype>(v.dims); \
        for (uint_t i = 0; i < v.size(); ++i) { \
            r.safe[i] = name(v.safe[i], args...); \
        } \
        return r; \
    }

    VECTORIZE(directorize)
    VECTORIZE(get_basename)
    VECTORIZE(get_directory)
    VECTORIZE(mkdir)
    VECTORIZE(exists)
    VECTORIZE(is_older)

#undef VECTORIZE

    template<typename T, typename ... Args>
    auto columns(std::size_t n, T& t, Args& ... args) ->
        decltype(std::tuple_cat(std::make_tuple(n), std::tie(t, args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(t, args...));
    }

    void read_table_resize_(std::size_t n) {}

    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args);
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args);
    template<typename ... Args>
    void read_table_resize_(std::size_t n, placeholder_t, Args& ... args);

    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec<1,T>& v, Args& ... args) {
        v = arr<T>(n);
        read_table_resize_(n, args...);
    }

    void read_table_resize_cols_(std::size_t n, std::size_t m) {}

    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, placeholder_t, VArgs&... args);

    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec<2,Type>& v, VArgs&... args) {
        v = arr<Type>(n, m);
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
        if (!(in >> v)) {
            in.clear();
            in.seekg(pos);
            in >> fallback;
            return false;
        }

        return true;
    }

    template<typename I>
    bool read_value_(I& in, double& v, std::string& fallback) {
        // std::istream will fail to extract a float/double value where the string is
        // 'INF', 'NAN', 'NULL', or any other text
        auto pos = in.tellg();
        if (!(in >> v)) {
            in.clear();
            in.seekg(pos);
            std::string s;
            in >> s;
            s = trim(toupper(s));
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

            fallback = s;
            return false;
        }

        return true;
    }

    template<typename I>
    bool read_value_(I& in, float& v, std::string& fallback) {
        // std::istream will fail to excract a float value if the exponent is too large
        // To fix this, we first try with float, and on failure, try again with double.
        // If exctraction is successful, then we just assign the double to the float, making it
        // infinite, else there is another error.
        auto pos = in.tellg();
        if (!(in >> v)) {
            in.clear();
            in.seekg(pos);
            double d;
            if (!read_value_(in, d, fallback)) {
                return false;
            } else {
                v = d;
                return true;
            }
        }

        return true;
    }

    void read_table_(std::istringstream& fs, std::size_t i, std::size_t& j) {}

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
            phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
            phypp_check(false, "cannot extract value '", fb, "' from file, wrong type for l."+
                strn(i)+":"+strn(j)+" (expected '"+std::string(typeid(T).name())+"'):\n"+fs.str());
        }
        read_table_(fs, i, ++j, args...);
    }

    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k) {}

    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, placeholder_t, VArgs&... args);

    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, vec<2,T>& v, VArgs&... args) {
        std::string fb;
        if (!read_value_(fs, v(i,k), fb)) {
            phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
            phypp_check(false, "cannot extract value '", fb, "' from file, wrong type for l."+
                strn(i)+":"+strn(j)+" (expected '"+std::string(typeid(T).name())+"'):\n"+fs.str());
        }
        read_table_cols_(fs, i, ++j, k, args...);
    }

    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t& j, std::size_t k, placeholder_t, VArgs&... args) {
        std::string s;
        if (!(fs >> s)) {
            phypp_check(!fs.eof(), "cannot extract value at l."+strn(i)+":"+strn(j)+" from file, "
                "too few columns");
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
        phypp_check(!fs.eof(), "cannot extract value at l."+strn(i)+":"+strn(j)+" from file, "
            "too few columns");
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

    void write_table_phypp_check_size_(std::size_t n, std::size_t i) {}

    template<std::size_t Dim, typename Type, typename ... Args>
    void write_table_phypp_check_size_(std::size_t& n, std::size_t i, const vec<Dim,Type>& v,
        const Args& ... args) {

        if (n == 0) {
            n = v.dims[0];
        }

        phypp_check(v.dims[0] == n, "incorrect dimension for column "+strn(i)+" ("+
            strn(v.dims[0])+" vs "+strn(n)+")");

        write_table_phypp_check_size_(n, i+1, args...);
    }

    void write_table_do_(std::ofstream& file, std::size_t cwidth, const std::string& sep,
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
        write_table_phypp_check_size_(n, t, args...);

        std::ofstream file(filename);
        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, "", i, 0, args...);
        }
    }

    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, const vec1s& nhdr,
        const Args& ... args) {

        std::size_t n = 0, t = 0;
        write_table_phypp_check_size_(n, t, args...);

        std::ofstream file(filename);

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

    std::string bake_macroed_name(const std::string& s) {
        std::string t = trim(s);
        auto p = t.find_first_of('.');
        if (p != t.npos) {
            t = t.substr(p+1);
        }
        return t;
    }

    template<typename ... Args>
    void write_table_hdr(const std::string& filename, std::size_t cwidth, macroed_t,
        const std::string& names, const Args& ... args) {

        vec1s vnames = split(names.substr(names.find_first_of(')')+1), ",");
        for (auto& s : vnames) {
            s = file::bake_macroed_name(s);
        }

        write_table_hdr(filename, cwidth, vnames, args...);
    }

    template<typename ... Args>
    void write_table_csv(const std::string& filename, std::size_t cwidth, const Args& ... args) {
        std::size_t n = 0, t = 0;
        write_table_phypp_check_size_(n, t,args...);

        std::ofstream file(filename);
        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, ",", i, 0, args...);
        }
    }
}

bool fork(const std::string& cmd) {
    return system((cmd+" &").c_str()) == 0;
}

bool spawn(const std::string& cmd) {
    return system(cmd.c_str()) == 0;
}

#endif

