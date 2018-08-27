#ifndef VIF_UTILITY_OS_HPP
#define VIF_UTILITY_OS_HPP

#include <string>
#include <cstdlib>
#include "vif/core/string_conversion.hpp"

namespace vif {
    template <typename T = std::string, typename U>
    T system_var(const std::string& name, U def) {
        char* v = getenv(name.c_str());
        if (!v) return def;

        T ret;
        if (!from_string(v, ret)) return def;
        return ret;
    }

    inline std::string system_var(const std::string& name, std::string def) {
        char* v = getenv(name.c_str());
        if (!v) return def;
        return v;
    }

    inline std::string system_var(const std::string& name, const char* def) {
        return system_var(name, std::string(def));
    }

    inline bool fork(const std::string& cmd) {
        return system((cmd+" &").c_str()) == 0;
    }

    inline bool spawn(const std::string& cmd) {
        return system(cmd.c_str()) == 0;
    }
}

#endif
