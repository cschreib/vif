#ifndef ARGV_HPP
#define ARGV_HPP

#include "vec.hpp"
#include "string.hpp"

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
            t.data.reserve(n_elements(vals));
            T tmp;
            for (auto& s : vals) {
                bool v = read_args_n2T_(tmp, trim(trim(s), "'\""));
                if (!v) {
                    t.data.clear();
                    valid = false;
                    return;
                } else {
                    t.data.push_back(tmp);
                }
            }

            t.dims[0] = t.data.size();
            t.data.shrink_to_fit();
            valid = true;
        } else {
            valid = read_args_n2T_(t, trim(trim(value), "'\""));
        }
    }
}

void read_args_(const vec1s& argv, vec1b& read, vec1b& valid, const std::string& names) {}

template<typename T, typename ... Args>
void read_args_(const vec1s& argv, vec1b& read, vec1b& valid, const std::string& names, T& t, Args&& ... args) {
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
    }
}

template<typename ... Args>
void read_args(uint_t argc, char* argv[], const std::string& names, Args&& ... args) {
    if (argc <= 1) return;

    vec1s sargv = strarr(argc-1);
    for (uint_t i = 1; i < argc; ++i) {
        sargv[i-1] = trim(argv[i]);
    }

    vec1b read  = boolarr(argc-1);
    vec1b valid = !read;

    read_args_(sargv, read, valid, names, std::forward<Args>(args)...);

    vec1u idm = where(!read);
    for (auto& i : idm) {
        warning("unrecognized program argument '", sargv[i],"'");
    }
}

#endif
