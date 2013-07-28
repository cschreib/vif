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
#include "vec.hpp"
#include "math.hpp"

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

    template<std::size_t Dim>
    vec_t<Dim,bool> exists(const vec_t<Dim,std::string>& file) {
        vec_t<Dim,bool> r = boolarr(file.dims);
        for (std::size_t i = 0; i < file.size(); ++i) {
            r.data[i] = file_exists(file.data[i]);
        }
        
        return r;
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

    bool mkdir(const std::string& path) {
        return (::mkdir(path.c_str(), 0775) == 0) || (errno == EEXIST);
    }
    
    template<typename U, typename ... Args>
    auto columns(U n, Args& ... args) -> decltype(std::tuple_cat(std::make_tuple(n), std::tie(args...))) {
        return std::tuple_cat(std::make_tuple(n), std::tie(args...));
    }
    
    void read_table_resize_(std::size_t n) {}
    
    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec_t<1,T>& v, Args& ... args);
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_resize_(std::size_t n, std::tuple<U,VArgs&...> v, Args& ... args);
    template<typename ... Args>
    void read_table_resize_(std::size_t n, placeholder_t, Args& ... args);
    
    template<typename T, typename ... Args>
    void read_table_resize_(std::size_t n, vec_t<1,T>& v, Args& ... args) {
        v = arr<T>(n);
        read_table_resize_(n, args...);
    }
    
    void read_table_resize_cols_(std::size_t n, std::size_t m) {}
    
    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec_t<2,Type>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, placeholder_t, VArgs&... args);
    
    template<typename Type, typename ... VArgs>
    void read_table_resize_cols_(std::size_t n, std::size_t m, vec_t<2,Type>& v, VArgs&... args) {
        v = arr<Type>(m, n);
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
    bool read_value_(I& in, T& v) {
        return in >> v;
    }

    template<typename I>
    bool read_value_(I& in, double& v) {
        // std::istream will fail to extract a float/double value where the string is 'INF' or 'NAN'
        auto pos = in.tellg();
        if (!(in >> v)) {
            in.clear();
            in.seekg(pos);
            std::string s;
            in >> s;
            s = trim(toupper(s));
            if (s == "NAN") {
                v = dnan;
                return true;
            } else if (s == "+INF" || s == "INF+" || s == "INF") {
                v = dinf;
                return true;
            } else if (s == "-INF" || s == "INF-") {
                v = -dinf;
                return true;
            }
            return false;
        }

        return true;
    }

    template<typename I>
    bool read_value_(I& in, float& v) {
        // std::istream will fail to excract a float value if the exponent is too large
        // To fix this, we first try with float, and on failure, try again with double.
        // If exctraction is successful, then we just assign the double to the float, making it
        // infinite, else there is another error.
        auto pos = in.tellg();
        if (!(in >> v)) {
            in.clear();
            in.seekg(pos);
            double d;
            if (!read_value_(in, d)) {
                return false;
            } else {
                v = d;
                return true;
            }
        }

        return true;
    }

    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j) {}
    
    template<typename T, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, vec_t<1,T>& v, Args& ... args);
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, std::tuple<U,VArgs&...> v, Args& ... args);
    template<typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, placeholder_t, Args& ... args);
    
    template<typename T, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, vec_t<1,T>& v, Args& ... args) {
        phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
        if (!read_value_(fs, v[i])) {
            phypp_check(false, "cannot extract value from file, wrong type for l."+strn(i)+":"+strn(j)+" (expected '"+std::string(typeid(T).name())+"'):\n"+fs.str());
        }
        read_table_(fs, i, j+1, args...);
    }
    
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k) {}
    
    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k, vec_t<2,T>& v, VArgs&... args);
    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k, placeholder_t, VArgs&... args);
    
    template<typename T, typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k, vec_t<2,T>& v, VArgs&... args) {
        phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
        if (!read_value_(fs, v(k,i))) {
            phypp_check(false, "cannot extract value from file, wrong type for l."+strn(i)+":"+strn(j)+" (expected '"+std::string(typeid(T).name())+"'):\n"+fs.str());
        }
        read_table_cols_(fs, i, j+1, k, args...);
    }
    
    template<typename ... VArgs>
    void read_table_cols_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k, placeholder_t, VArgs&... args) {
        phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
        std::string s;
        fs >> s;
        read_table_cols_(fs, i, j+1, k, args...);
    }

    template<typename U, typename ... VArgs, std::size_t ... S>
    void read_table_cols_i_(std::istringstream& fs, std::size_t i, std::size_t j, std::size_t k, std::tuple<U,VArgs&...>& v, seq_t<S...>) {
        read_table_cols_(fs, i, j, k, std::get<S>(v)...);
    }
    
    template<typename U, typename ... VArgs, typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, std::tuple<U,VArgs&...> v, Args& ... args) {
        std::size_t n = std::get<0>(v);
        for (std::size_t k = 0; k < n; ++k) {
            read_table_cols_i_(fs, i, j, k, v, typename gen_seq<1, sizeof...(VArgs)>::type());
        }
        
        read_table_(fs, i, j, args...);
    }
    
    template<typename ... Args>
    void read_table_(std::istringstream& fs, std::size_t i, std::size_t j, placeholder_t, Args& ... args) {
        phypp_check(!fs.eof(), "cannot extract value from file, too few columns");
        std::string s;
        fs >> s;
        read_table_(fs, i, j+1, args...);
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
            read_table_(fs, i, 0, args...);
        }
    }
    
    void write_table_phypp_check_size_(std::size_t n) {}
    
    template<typename Type, typename ... Args>
    void write_table_phypp_check_size_(std::size_t& n, const vec_t<1,Type>& v, const Args& ... args) {
        if (n == 0) {
            n = v.size();
        }
        
        assert(v.size() == n);
        
        write_table_phypp_check_size_(n, args...);
    }
    
    void write_table_do_(std::ofstream& file, std::size_t cwidth, std::size_t i) {
        file << '\n';
    }
    
    template<typename Type, typename ... Args>
    void write_table_do_(std::ofstream& file, std::size_t cwidth, std::size_t i, const vec_t<1,Type>& v, const Args& ... args) {
        std::string s = " "+strn(v[i]);
        if (s.size() < cwidth) {
            s += std::string(cwidth - s.size(), ' ');
        }
        file << s;
        
        write_table_do_(file, cwidth, i, args...);
    }
    
    template<typename ... Args>
    void write_table(const std::string& filename, std::size_t cwidth, const Args& ... args) {
        std::size_t n = 0;
        write_table_phypp_check_size_(n, args...);
    
        std::ofstream file(filename);
        for (std::size_t i = 0; i < n; ++i) {
            write_table_do_(file, cwidth, i, args...);
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

