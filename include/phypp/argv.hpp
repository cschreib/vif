#ifndef ARGV_HPP
#define ARGV_HPP

#include "phypp/vec.hpp"
#include "phypp/string.hpp"
#include <fstream>

#define arg_list(...) #__VA_ARGS__, __VA_ARGS__

template<typename T, typename U,
    typename enable = typename std::enable_if<!std::is_same<T, vec_t<1,std::string>>::value>::type>
bool read_args_n2T_(T& t, const U& u) {
    t = u;
    return true;
}

template<typename T,
    typename enable = typename std::enable_if<!std::is_same<T, vec_t<1,std::string>>::value>::type>
bool read_args_n2T_(T& t, const std::string& s) {
    return from_string(s, t);
}

template<typename U>
bool read_args_n2T_(vec_t<1,std::string>& t, const U& u) {
    t = strn(u);
    return true;
}

template<typename U>
bool read_args_n2T_(vec_t<1,std::string>& t, const std::string& u) {
    t = u;
    return true;
}

template<typename U>
bool read_args_n2T_(std::string& t, const U& u) {
    t = strn(u);
    return true;
}

bool read_args_n2T_(std::string& t, const std::string& s) {
    t = s;
    return true;
}

template<typename T, typename enable = typename std::enable_if<!is_vec<T>::value>::type>
void read_args_impl_(const std::string& arg, bool& read, bool& valid, const std::string& name, T& t) {
    auto p = arg.find_first_of('=');
    if (p == arg.npos) {
        if (arg != name) {
            return;
        }

        read = true;
        valid = read_args_n2T_(t, 1);
    } else {
        std::string aname = trim(arg.substr(0, p));
        if (aname != name) {
            return;
        }

        read = true;
        std::string value = trim(trim(arg.substr(p+1)), "'\"");
        valid = read_args_n2T_(t, value);
    }
}

template<typename T>
void read_args_impl_(const std::string& arg, bool& read, bool& valid, const std::string& name, vec_t<1,T>& t) {
    auto p = arg.find_first_of('=');
    if (p == arg.npos) {
        if (arg != name) {
            return;
        }

        read = true;
        valid = read_args_n2T_(t, 1);
    } else {
        std::string aname = trim(arg.substr(0, p));
        if (aname != name) {
            return;
        }

        read = true;
        std::string value = trim(arg.substr(p+1));
        if (!value.empty() && value.front() == '[' && value.back() == ']') {
            value.erase(0,1); value.pop_back();
            vec1s vals = split(value, ",");
            t.clear();
            t.reserve(n_elements(vals));
            rtype_t<T> tmp;
            for (auto& s : vals) {
                bool v = read_args_n2T_(tmp, trim(trim(s), "'\""));
                if (!v) {
                    t.clear();
                    valid = false;
                    return;
                } else {
                    t.push_back(tmp);
                }
            }

            t.data.shrink_to_fit();
            valid = true;
        } else {
            rtype_t<T> v;
            valid = read_args_n2T_(v, trim(trim(value), "'\""));
            if (valid) {
                t.clear();
                t.push_back(v);
            }
        }
    }
}

void read_args_(const vec1s& argv, vec1b& read, vec1b& valid, const std::string& names) {}

template<typename T, typename ... Args>
void read_args_(const vec1s& argv, vec1b& read, vec1b& valid, const std::string& names, T& t,
    Args&& ... args) {

    std::size_t pos = names.find_first_of(',');
    std::string name = trim(names.substr(0, pos));
    auto p = name.find_last_of('.');
    if (p != name.npos) {
        name = trim(name.substr(p+1));
    }

    vec1u idm = where(find(argv, name) != npos);
    for (auto& i : idm) {
        read_args_impl_(argv[i], read[i], valid[i], name, t);
        if (!valid[i]) {
            warning("could not convert '", name, "' argument to type ", typeid(T).name());
            break;
        }
    }

    if (pos != names.npos) {
        read_args_(argv, read, valid, names.substr(pos+1), std::forward<Args>(args)...);
    } else if (sizeof...(Args) != 0) {
        error("read_args: too few names provided");
        note("please use the arg_list() macro and make sure that all variables are "
            "placed there");
        throw std::logic_error("read_args: too few names provided");
    }
}

template<typename T, typename ... Args>
void read_args_(const vec1s& argv, vec1b& read, vec1b& valid, const std::string& names,
    named_t<T> t, Args&& ... args) {

    std::size_t pos = names.find_first_of(')');
    if (pos != names.npos) ++pos;

    vec1u idm = where(find(argv, t.name) != npos);
    for (auto& i : idm) {
        read_args_impl_(argv[i], read[i], valid[i], t.name, t.obj);
        if (!valid[i]) {
            warning("could not convert '", t.name, "' argument to type ", typeid(T).name());
            break;
        }
    }

    if (pos != names.npos && pos != names.size()) {
        read_args_(argv, read, valid, names.substr(pos+1), std::forward<Args>(args)...);
    }
}

struct program_arguments {
    program_arguments(int argc, char* argv[]) {
        if (argc <= 1) return;

        uint_t narg = argc;
        argv_ = strarr(narg-1);
        for (uint_t i = 1; i < narg; ++i) {
            argv_[i-1] = trim(argv[i]);
        }

        read_.resize(narg-1);
        valid_ = !read_;
    }

    explicit program_arguments(const vec1s& args) {
        argv_ = args;
        read_.resize(args.size());
        valid_ = !read_;
    }

    ~program_arguments() {
        vec1u idm = where(!read_);
        for (auto& i : idm) {
            warning("unrecognized program argument '", argv_[i],"'");
        }
    }

    void merge(const program_arguments& pa) {
        append(argv_, pa.argv_);
        append(read_, pa.read_);
        append(valid_, pa.valid_);
    }

    template<typename ... Args>
    void read(const std::string& names, Args&& ... args) {
        read_args_(argv_, read_, valid_, names, std::forward<Args>(args)...);
    }

private :
    vec1s argv_;
    vec1b read_;
    vec1b valid_;
};

template<typename ... Args>
void read_args(uint_t argc, char* argv[], const std::string& names, Args&& ... args) {
    program_arguments pa(argc, argv);
    pa.read(names, std::forward<Args>(args)...);
}

void save_args(const std::string& file, const std::string& pname, int argc, char* argv[]) {
    std::ofstream cmd(file);
    cmd << pname;
    for (int_t i = 1; i < argc; ++i) {
        cmd << " " << argv[i];
    }
}

std::string make_cmd(int argc, char* argv[]) {
    std::string cmd = argv[0];
    for (int i = 1; i < argc; ++i) {
        std::string tmp = argv[i];
        if (tmp.find_first_of(" \t") != tmp.npos) {
            uint_t p = tmp.find_first_of("=");
            if (p == tmp.npos) {
                tmp = "\"" + tmp + "\"";
            } else {
                tmp = tmp.substr(0, p+1) + "\"" + tmp.substr(p+1) + "\"";
            }
        }
        cmd += " "+tmp;
    }

    return cmd;
}

#endif
