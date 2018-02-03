#ifndef PHYPP_CORE_STRING_CONVERSION_HPP
#define PHYPP_CORE_STRING_CONVERSION_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <typeinfo>
#include "phypp/core/vec.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/meta.hpp"

namespace std {
    template<typename O, typename T, std::size_t N>
    O& operator << (O& o, const std::array<T,N>& v) {
        o << "{";
        for (phypp::uint_t i : phypp::range(N)) {
            if (i != 0) o << ", ";
            o << v[i];
        }
        o << "}";

        return o;
    }
}


namespace phypp {
    namespace impl {
        template<typename T>
        struct named_t {
            T&          obj;
            std::string name;
        };
    }

    template<typename T>
    impl::named_t<T> name(T& obj, const std::string& str) {
        return {obj, str};
    }

    template<typename T>
    std::string strn(const T& t) {
        std::ostringstream ss;
        ss << t;
        return ss.str();
    }

    inline std::string strn(double t) {
        std::ostringstream ss;
        ss.precision(12);
        ss << t;
        return ss.str();
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,std::string> strna(const vec<Dim,Type>& v) {
        vec<Dim,std::string> s(v.dims);
        for (uint_t i : range(v)) {
            s.safe[i] = strn(v.safe[i]);
        }

        return s;
    }

    template<typename T>
    std::string strn_sci(const T& t) {
        std::ostringstream ss;
        ss << std::scientific << t;
        return ss.str();
    }


    template<std::size_t Dim, typename Type>
    vec<Dim,std::string> strna_sci(const vec<Dim,Type>& v) {
        vec<Dim,std::string> s(v.dims);
        for (uint_t i : range(v)) {
            s.safe[i] = strn_sci(v.safe[i]);
        }

        return s;
    }

    template<typename T>
    bool from_string(const std::string& s, T& t) {
        std::istringstream ss(s);
        ss >> t;

        if (!ss.fail()) {
            if (ss.eof()) return true;
            std::string rem;
            ss >> rem;
            return rem.find_first_not_of(" \t") == rem.npos;
        } else {
            return false;
        }
    }

    template<std::size_t Dim = 1, typename T = std::string, typename O,
    typename enable = typename std::enable_if<
        std::is_same<meta::rtype_t<T>, std::string>::value &&
        meta::is_compatible_output_type<vec<Dim,T>,O>::value
    >::type>
    vec<Dim,bool> from_string(const vec<Dim,T>& s, O&& t) {
        vec<Dim,bool> res(s.dims);
        meta::resize_or_check(t, s.dims);
        for (uint_t i : range(s)) {
            res.safe[i] = from_string(s.safe[i], t.safe[i]);
        }

        return res;
    }

    namespace impl {
        template<typename T>
        std::string pretty_type_(meta::type_list<T>) {
            return typeid(T).name();
        };

        template<typename T>
        std::string pretty_type_(meta::type_list<T&>) {
            return pretty_type_(meta::type_list<T>{})+"&";
        };

        template<typename T>
        std::string pretty_type_(meta::type_list<T*>) {
            return pretty_type_(meta::type_list<T>{})+"*";
        };

        template<typename T>
        std::string pretty_type_(meta::type_list<const T>) {
            return "const "+pretty_type_(meta::type_list<T>{});
        };

        inline std::string pretty_type_(meta::type_list<char>) {
            return "char";
        };

        inline std::string pretty_type_(meta::type_list<short>) {
            return "short";
        };

        inline std::string pretty_type_(meta::type_list<int>) {
            return "int";
        };

        inline std::string pretty_type_(meta::type_list<long>) {
            return "long";
        };

        inline std::string pretty_type_(meta::type_list<long long>) {
            return "llong";
        };

        inline std::string pretty_type_(meta::type_list<unsigned char>) {
            return "uchar";
        };

        inline std::string pretty_type_(meta::type_list<unsigned short>) {
            return "ushort";
        };

        inline std::string pretty_type_(meta::type_list<unsigned int>) {
            return "uint";
        };

        inline std::string pretty_type_(meta::type_list<unsigned long>) {
            return "ulong";
        };

        inline std::string pretty_type_(meta::type_list<unsigned long long>) {
            return "ullong";
        };

        inline std::string pretty_type_(meta::type_list<bool>) {
            return "bool";
        };

        inline std::string pretty_type_(meta::type_list<float>) {
            return "float";
        };

        inline std::string pretty_type_(meta::type_list<double>) {
            return "double";
        };

        inline std::string pretty_type_(meta::type_list<std::string>) {
            return "string";
        };

        template<std::size_t Dim, typename T>
        std::string pretty_type_(meta::type_list<vec<Dim,T>>) {
            return "vec<"+strn(Dim)+","+pretty_type_(meta::type_list<T>{})+">";
        };
    }

    #define pretty_type(x) phypp::impl::pretty_type_(meta::type_list<decltype(x)>{})
    #define pretty_type_t(x) phypp::impl::pretty_type_(meta::type_list<x>{})
}

#endif

