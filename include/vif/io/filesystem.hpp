#ifndef VIF_IO_FILESYSTEM_HPP
#define VIF_IO_FILESYSTEM_HPP

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
#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/utility/string.hpp"

namespace vif {
namespace file {
    inline bool exists(const std::string& path) {
        if (path.empty()) {
            return false;
        }

        std::ifstream f(path.c_str());
        return f.is_open();
    }

    inline bool is_older(const std::string& file1, const std::string& file2) {
        struct stat st1, st2;
        if (::stat(file1.c_str(), &st1) != 0) return false;
        if (::stat(file2.c_str(), &st2) != 0) return false;
        return std::difftime(st1.st_ctime, st2.st_ctime) < 0.0;
    }

    inline bool is_absolute_path(const std::string& path) {
        auto pos = path.find_first_not_of(" \t");
        return pos != path.npos && path[pos] == '/';
    }

    inline bool copy(const std::string& file_from, const std::string& file_to) {
        std::ifstream src(file_from, std::ios::binary);
        if (!src.is_open()) return false;
        std::ofstream dst(file_to,   std::ios::binary);
        if (!dst.is_open()) return false;
        dst << src.rdbuf();
        return true;
    }

    inline bool remove(const std::string& path) {
        return !file::exists(path) || ::remove(path.c_str()) == 0;
    }

    inline bool move(const std::string& from, const std::string& to) {
        return ::rename(from.c_str(), to.c_str()) == 0;
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

    class explorer {
        std::string directory;
        std::string pattern;
        DIR* dirfd = nullptr;

    public :

        struct file_data {
            std::string full_path;
            std::string name;
            uint_t size;
            bool is_hidden = false;
            bool is_dir = false;
        };

        explorer() = default;
        explorer(const explorer&) = delete;
        explorer(explorer&&) = delete;
        explorer& operator=(const explorer&) = delete;
        explorer& operator=(explorer&&) = delete;

        explicit explorer(const std::string& dir, const std::string& p = "") {
            open(dir, p);
        }

        ~explorer() {
            close();
        }

        void close() {
            directory = "";
            pattern = "";

            if (dirfd) {
                closedir(dirfd);
            };
        }

        bool open(const std::string& dir, const std::string& p = "") {
            close();

            dirfd = opendir(dir.c_str());
            if (!dirfd) {
                return false;
            }

            directory = dir;
            pattern = p;
            if (pattern == "*") pattern = "";

            return true;
        }

        bool find_next(file_data& f) {
            if (!dirfd) return false;

            dirent* entry;

            if (pattern.empty()) {
                entry = readdir(dirfd);
            } else {
                while ((entry = readdir(dirfd))) {
                    if (fnmatch(pattern.c_str(), entry->d_name, 0) == 0) {
                        break;
                    }
                }
            }

            if (!entry) return false;

            f.name = entry->d_name;
            f.full_path = directory + "/" + f.name;

            struct stat stat_buf;
            if (stat(f.full_path.c_str(), &stat_buf)) {
                f.is_dir = false;
                f.size = 0;
            } else {
                f.is_dir = S_ISDIR(stat_buf.st_mode);
                f.size = stat_buf.st_size;
            }

            f.is_hidden = (f.name.size() > 1 && f.name[0] == '.');

            return true;
        }
    };

    inline vec1s list_directories(const std::string& dir, const std::string& pattern = "*") {
        vec1s dlist;

        explorer e;
        if (!e.open(dir, pattern)) {
            return dlist;
        }

        explorer::file_data f;
        while (e.find_next(f)) {
            if (!f.is_hidden && f.is_dir && f.name != "." && f.name != "..") {
                dlist.push_back(f.name);
            }
        }

        return dlist;
    }

    inline vec1s list_files(const std::string& dir, const std::string& pattern = "*") {
        vec1s flist;

        explorer e;
        if (!e.open(dir, pattern)) {
            return flist;
        }

        explorer::file_data f;
        while (e.find_next(f)) {
            if (!f.is_hidden && !f.is_dir) {
                flist.push_back(f.name);
            }
        }

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
            if (!tmp.empty() || begins_with(path, "/")) tmp += "/";
            tmp += d;
            bool res = (::mkdir(tmp.c_str(), 0775) == 0) || (errno == EEXIST);
            if (!res) return false;
        }
        return true;
    }

    VIF_VECTORIZE(directorize)
    VIF_VECTORIZE(is_absolute_path)
    VIF_VECTORIZE(get_basename)
    VIF_VECTORIZE(get_directory)
    VIF_VECTORIZE(remove_extension)
    VIF_VECTORIZE(get_extension)
    VIF_VECTORIZE(split_extension)
    VIF_VECTORIZE(mkdir)
    VIF_VECTORIZE(exists)
    VIF_VECTORIZE(is_older)
}
}


#endif

