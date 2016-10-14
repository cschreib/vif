#ifndef PHYPP_IO_FILESYSTEM_HPP
#define PHYPP_IO_FILESYSTEM_HPP

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
#include <stdexcept>
#include "phypp/core/vec.hpp"
#include "phypp/core/error.hpp"
#include "phypp/math/math.hpp"

namespace phypp {
namespace impl {
    namespace file_impl {
        // Emulate directory and file listing windows functions
        // Note : Taken directly from Ogre3D
        // http://www.ogre3d.org/
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
}

namespace file {
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
        struct impl::file_impl::_finddata_t tagData;

        handle = impl::file_impl::_findfirst(pattern.c_str(), &tagData);
        res = 0;
        while (handle != -1 && res != -1) {
            if ((tagData.attrib & _A_HIDDEN) != 0) {
                res = impl::file_impl::_findnext(handle, &tagData);
                continue;
            }

            if ((tagData.attrib & _A_SUBDIR) != 0) {
                std::string s = tagData.name;
                if (s != "." && s != "..") {
                    dlist.data.push_back(s);
                }
            }

            res = impl::file_impl::_findnext(handle, &tagData);
        }

        if (handle != -1) {
            impl::file_impl::_findclose(handle);
        }

        dlist.dims[0] = dlist.size();
        return dlist;
    }

    inline vec1s list_files(const std::string& pattern = "*") {
        vec1s flist;

        long handle, res;
        struct impl::file_impl::_finddata_t tagData;

        handle = impl::file_impl::_findfirst(pattern.c_str(), &tagData);
        res = 0;
        while (handle != -1 && res != -1) {
            if ((tagData.attrib & _A_HIDDEN) != 0) {
                res = impl::file_impl::_findnext(handle, &tagData);
                continue;
            }

            if ((tagData.attrib & _A_SUBDIR) == 0) {
                flist.data.push_back(tagData.name);
            }

            res = impl::file_impl::_findnext(handle, &tagData);
        }

        if (handle != -1) {
            impl::file_impl::_findclose(handle);
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
}
    
    inline bool fork(const std::string& cmd) {
        return system((cmd+" &").c_str()) == 0;
    }

    inline bool spawn(const std::string& cmd) {
        return system(cmd.c_str()) == 0;
    }
}


#endif

