#ifndef PRINT_HPP
#define PRINT_HPP

#include <string>
#include <iostream>
#include <array>
#include "phypp/reflex.hpp"

template<typename T>
void print(const T& t) {
    std::cout << reflex::wrap(t) << std::endl;
}

template<typename T, typename ... Args>
void print(const T& t, Args&&... args) {
    std::cout << t;
    print(std::forward<Args>(args)...);
}

template<typename ... Args>
void error(Args&& ... args) {
    print("error: ", std::forward<Args>(args)...);
}

template<typename ... Args>
void warning(Args&& ... args) {
    print("warning: ", std::forward<Args>(args)...);
}

template<typename ... Args>
void note(Args&& ... args) {
    print("note: ", std::forward<Args>(args)...);
}

template<typename T>
bool prompt(const std::string& msg, T& value, const std::string& err_msg = "") {
    bool ok = false;

    do {
        std::cout << msg;
        std::cin >> value;

        if (!std::cin.fail()) {
            if (std::cin.eof()) {
                ok = true;
            } else {
                std::string rem;
                std::cin >> rem;
                ok = rem.find_first_not_of(" \t") == rem.npos;
            }
        } else {
            ok = false;
        }

        if (!ok) {
            std::cin.clear();
            std::cin.ignore();
            if (!err_msg.empty()) {
                error(err_msg);
            }
        }
    } while (!ok);

    return ok;
}

#define phypp_check(value, ...) \
    if (!(value)) { \
        error(__FILE__, ":", __LINE__, ": ", __VA_ARGS__); \
        assert(false); \
    }

#endif

