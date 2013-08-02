#ifndef REFLEX_HPP
#define REFLEX_HPP

#include <string>

template<std::size_t Dim, typename Type>
struct vec_t;

using int_t = std::ptrdiff_t;
using uint_t = std::size_t;

struct rgb;

namespace reflex {
    template<typename ... Args>
    struct type_list {
        template<std::size_t N, std::size_t I, typename ... Args2>
        struct get_t;

        template<std::size_t N, std::size_t I, typename T, typename ... Args2>
        struct get_t<N, I, T, Args2...> {
            using type = typename get_t<N, I+1, Args2...>::type;
        };

        template<std::size_t N, typename T, typename ... Args2>
        struct get_t<N, N, T, Args2...> {
            using type = T;
        };

        template<std::size_t N, typename enable =
            typename std::enable_if<N < sizeof...(Args)>::type>
        using get = typename get_t<N, 0, Args...>::type;

        static const std::size_t count = sizeof...(Args);
    };

    template<typename ... Args>
    type_list<Args...> make_types_(Args&& ...);

    template<typename ... Args>
    type_list<typename std::decay<Args>::type...> make_types_decay_(Args&& ...);

    struct data_t;

    struct member_t {
        std::string name;
        data_t* parent = nullptr;
        void* value;

        std::string full_name() const;
    };

    struct data_t {
        data_t() = default;
        data_t(data_t&&) = default;
        data_t(const data_t&) = delete;
        data_t& operator = (data_t&&) = default;
        data_t& operator = (const data_t&) = delete;

        data_t* parent = nullptr;
        std::string name;
        std::vector<member_t> members;

        std::string full_name() const;
    };

    std::string member_t::full_name() const {
        data_t* p = parent;
        std::string str = name;
        while (p) {
            str = p->name + "." + str;
            p = p->parent;
        }

        return str;
    }

    std::string data_t::full_name() const {
        data_t* p = parent;
        std::string str = name;
        while (p) {
            str = p->name + "." + str;
            p = p->parent;
        }

        return str;
    }

    template<typename T>
    struct has_reflex {
        template <typename U> static std::true_type dummy(decltype(((U*)0)->_reflex)*);
        template <typename U> static std::false_type dummy(...);
        static const bool value = decltype(dummy<T>(0))::value;
    };

    template<typename T>
    void do_init_set_name_(member_t& m, T& d, data_t* parent, cte_t<true>) {
        d._reflex.name = m.name;
        d._reflex.parent = parent;
    }

    template<typename T>
    void do_init_set_name_(member_t& m, T& d, data_t* parent, cte_t<false>) {}

    template<std::size_t N, typename ... Args>
    auto& get_value(type_list<Args...> tl, member_t& m) {
        using type = typename type_list<Args...>::template get<N>;
        return *(type*)m.value;
    }

    template<std::size_t N, typename ... Args>
    auto& get_value(type_list<Args...> tl, const member_t& m) {
        using type = typename type_list<Args...>::template get<N>;
        return *(const type*)m.value;
    }

    template<typename ... Args, typename T>
    void do_init_(type_list<Args...> tl, cte_t<0>, T* t, data_t& data) {}

    template<typename ... Args, std::size_t N, typename T>
    void do_init_(type_list<Args...> tl, cte_t<N>, T* t, data_t& data) {
        member_t& m = data.members[N-1];
        m.parent = &t->_reflex;
        auto& v = get_value<N-1>(tl, m);
        do_init_set_name_(m, v, &t->_reflex,
            cte_t<has_reflex<typename std::decay<decltype(v)>::type>::value>()
        );
        do_init_(tl, cte_t<N-1>(), t, data);
    }

    template<typename ... Args, typename T>
    data_t do_init(type_list<Args...> tl, T* t, data_t&& data) {
        data_t d = std::move(data);
        do_init_(tl, cte_t<sizeof...(Args)>(), t, d);
        return d;
    }

    template<typename T>
    struct enabled : std::is_class<T> {};

    template<typename T>
    struct enabled<T*> : std::false_type {};

    template<std::size_t D, typename T>
    struct enabled<vec_t<D,T>> : std::false_type {};

    template<std::size_t D, typename T>
    struct enabled<std::array<T,D>> : std::false_type {};

    template<>
    struct enabled<std::string> : std::false_type {};

    template<>
    struct enabled<rgb> : std::false_type {};

    struct empty_t {
        using _reflex_types = type_list<>;
        data_t _reflex;
    };

    template<typename T, typename U>
    using constify = typename std::conditional<std::is_const<T>::value, const U, U>::type;

    template<typename T>
    struct struct_t {
        constify<T,data_t>& data;
        using member_types = typename T::_reflex_types;
        static const std::size_t member_count = T::_reflex_types::count;
    };

    template<typename T>
    struct is_struct : std::false_type {};

    template<typename T>
    struct is_struct<struct_t<T>> : std::true_type {};

    template<std::size_t N, typename T>
    auto& get_value(struct_t<T>& t) {
        using type = typename struct_t<T>::member_types::template get<N>;
        return *(type*)t.data.members[N].value;
    }

    template<std::size_t N, typename T>
    auto& get_value(const struct_t<T>& t) {
        using type = typename struct_t<T>::member_types::template get<N>;
        return *(const type*)t.data.members[N].value;
    }

    template<bool reflexed>
    struct wrap_t;

    template<>
    struct wrap_t<true> {
        #ifdef REFLECTION_STAGE
        template<typename T>
        static auto wrap(T& t) {
            static empty_t empty;
            return struct_t<empty_t>{empty._reflex};
        }

        #else
        template<typename T>
        static auto wrap(T& t) {
            return struct_t<T>{t._reflex};
        }
        #endif
    };

    template<>
    struct wrap_t<false> {
        template<typename T>
        static T& wrap(T& t) {
            return t;
        }

        template<typename T>
        static T* wrap(T* t) {
            return t;
        }
    };

    template<typename T>
    auto wrap(T&& t) -> decltype(auto) {
        return wrap_t<enabled<typename std::decay<T>::type>::value>::wrap(std::forward<T>(t));
    }
}

#define MAKE_MEMBER(name) reflex::member_t{#name, nullptr, (void*)&name}

#define MEMBERS1(...) using _reflex_types = decltype(reflex::make_types_decay_(__VA_ARGS__))

#define MEMBERS2(name, ...) \
    reflex::data_t _reflex = reflex::do_init(_reflex_types(), this, {nullptr, name, {__VA_ARGS__}})

#define NO_MEMBER(name) \
    using _reflex_types = decltype(reflex::make_types_()); \
    reflex::data_t _reflex = {nullptr, name, {}}

#endif
