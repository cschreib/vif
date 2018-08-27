#ifndef VIF_TEST_UNIT_TEST_HPP
#define VIF_TEST_UNIT_TEST_HPP

#include "vif/core/vec.hpp"
#include "vif/core/print.hpp"
#include "vif/core/string_conversion.hpp"
#include "vif/math/math.hpp"

namespace vif {
    uint_t tested = 0u;
    uint_t failed = 0u;

    bool reduce_and(bool b) {
        return b;
    }

    template<std::size_t Dim>
    bool reduce_and(const vec<Dim,bool>& b) {
        bool res = true;
        for (bool v : b) {
            res = res && v;
        }
        return res;
    }

    namespace impl {
        template<typename T1>
        struct is_float_ : std::false_type {
        };
        template<>
        struct is_float_<float> : std::true_type {
        };
        template<>
        struct is_float_<double> : std::true_type {
        };
        template<std::size_t D, typename T>
        struct is_float_<vec<D,T>> : is_float_<typename std::remove_pointer<T>::type> {
        };
    }

    template<typename T1, typename T2>
    struct is_float : std::integral_constant<bool, impl::is_float_<T1>::value || impl::is_float_<T2>::value> {
    };

    namespace impl {
        template<typename T1, typename T2>
        auto is_same_(const T1& v1, const T2& v2, std::false_type) -> decltype(v1 == v2) {
            return v1 == v2;
        }

        template<typename T1, typename T2>
        auto is_same_(const T1& v1, const T2& v2, std::true_type) -> decltype(v1 == v2) {
            return abs(v1 - v2) <= std::numeric_limits<float>::epsilon() ||
                (is_nan(v1) && is_nan(v2)) || is_nan(v1 - v2);
        }
    }

    template<typename T1, typename T2>
    auto is_same(const T1& v1, const T2& v2) -> decltype(v1 == v2) {
        return impl::is_same_(v1, v2, is_float<T1,T2>{});
    }

    static bool check_show_line = false;

    #define check_base_(cond, msg, file, line) do { \
        if (check_show_line) print(std::string(file)+":"+to_string(line)); \
        if (!cond) { \
            print(msg); \
            ++failed; \
        } \
        ++tested; \
    } while (false)

    #define check_base(cond, msg) check_base_(reduce_and(cond), msg, __FILE__, __LINE__)

    #define check(t, s) \
        check_base(is_same(t, s), "  failed: "+std::string(#t)+" = "+to_string(t)+" != "+to_string(s))
}

#endif
