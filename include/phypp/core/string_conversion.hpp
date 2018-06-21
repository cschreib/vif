#ifndef PHYPP_CORE_STRING_CONVERSION_HPP
#define PHYPP_CORE_STRING_CONVERSION_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <array>
#include <typeinfo>
#include <cctype>
#include <limits>
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

    namespace impl {
        struct format_base_t {};

        template<typename F, typename T>
        struct format_t : format_base_t {
            using object_type = typename std::decay<T>::type;

            template<typename U, typename ... Args>
            explicit format_t(U&& o, Args&& ... args) :
                fmt(std::forward<Args>(args)...), obj(std::forward<U>(o)) {}

            F fmt;
            const T& obj;

            template<typename O>
            void apply(O& o) const {
                fmt.apply(o);
            }

            template<typename U>
            format_t<F,U> forward(U&& u) const {
                return format_t<F,U>(std::forward<U>(u), fmt);
            }
        };

        struct format_scientific_t {
            template<typename O>
            void apply(O& o) const {
                o << std::scientific;
            }
        };

        struct format_precision_t {
            uint_t pre = 0;

            explicit format_precision_t(uint_t p) : pre(p) {}

            template<typename O>
            void apply(O& o) const {
                o << std::setprecision(pre);
            }
        };
    }

    namespace meta {
        template<typename T>
        struct is_format_tag : std::is_base_of<impl::format_base_t, typename std::decay<T>::type> {};
    }

    namespace impl {
        template<typename O, typename F,
            typename enable = typename std::enable_if<meta::is_format_tag<F>::value>::type>
        O& operator << (O& o, const F& v) {
            std::ios old(nullptr);
            old.copyfmt(o);

            v.apply(o);
            o << v.obj;

            o.copyfmt(old);

            return o;
        }
    }

    namespace format {
        template<typename T>
        impl::format_t<impl::format_scientific_t,T> scientific(T&& obj) {
            return impl::format_t<impl::format_scientific_t,T>(std::forward<T>(obj));
        }

        template<typename T>
        impl::format_t<impl::format_precision_t,T> precision(T&& obj, uint_t p) {
            return impl::format_t<impl::format_precision_t,T>(std::forward<T>(obj), p);
        }
    }

    template<typename T>
    std::string to_string(const T& t) {
        std::ostringstream ss;
        ss << t;
        return ss.str();
    }

    inline std::string to_string(const std::string& t) {
        return t;
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,std::string> to_string_vector(const vec<Dim,Type>& v) {
        vec<Dim,std::string> s(v.dims);
        for (uint_t i : range(v)) {
            s.safe[i] = to_string(v.safe[i]);
        }

        return s;
    }

    template<typename F, typename enable = typename std::enable_if<
        meta::is_format_tag<F>::value && meta::is_vec<typename F::object_type>::value>::type>
    auto to_string_vector(const F& f) ->
        vec<meta::vec_dim<typename F::object_type>::value,std::string> {
        vec<meta::vec_dim<typename F::object_type>::value,std::string> s(f.obj.dims);
        for (uint_t i : range(f.obj)) {
            s.safe[i] = to_string(f.forward(f.obj.safe[i]));
        }

        return s;
    }

namespace impl {
    template<typename T>
    bool from_string_fallback(const std::string& s, T& t, std::true_type) {
        std::string su = s;
        for (auto& c : su) c = std::toupper(c);

        if (su == "NAN" || su == "+NAN" || su == "-NAN") {
            t = std::numeric_limits<T>::quiet_NaN();
            return true;
        } else if (su == "+INF" || su == "INF+" || su == "INF") {
            t = std::numeric_limits<T>::infinity();
            return true;
        } else if (su == "-INF" || su == "INF-") {
            t = -std::numeric_limits<T>::infinity();
            return true;
        }

        return false;
    }

    template<typename T>
    bool from_string_fallback(const std::string&, T&, std::false_type) {
        return false;
    }
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
            // Try special strings
            return impl::from_string_fallback(s, t, std::is_floating_point<T>{});
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
            return "vec<"+to_string(Dim)+","+pretty_type_(meta::type_list<T>{})+">";
        };
    }

    #define pretty_type(x) phypp::impl::pretty_type_(meta::type_list<decltype(x)>{})
    #define pretty_type_t(x) phypp::impl::pretty_type_(meta::type_list<x>{})
}

#endif

