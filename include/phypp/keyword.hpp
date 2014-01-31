#ifndef KEYWORD_HPP
#define KEYWORD_HPP

#include "phypp/variadic.hpp"
#include <tuple>

// Check if T = V<X,W> and U = V<Y,W>
template<typename T, typename U>
struct is_same_tp {
    static const bool value = false;
};

template<template<typename,typename> class Tp, typename T, typename U, typename X>
struct is_same_tp<Tp<T,X>, Tp<U,X>> {
    static const bool value = true;
};

// Check if the provided template type T<X> matches one of the other types (regardless of X)
template<typename T, typename U, typename ... Args>
struct is_any_of_tp {
    static const bool value = is_same_tp<T,U>::value || is_any_of_tp<T,Args...>::value;
};

template<typename T, typename U>
struct is_any_of_tp<T,U> {
    static const bool value = is_same_tp<T,U>::value;
};

// Get an element of a tuple by template type (T<X>)
template<typename T, typename Tp, std::size_t N>
auto& tuple_get_tp_(Tp& t, cte_t<true>, cte_t<N>) {
    return std::get<N>(t);
}

template<typename T, typename Tp>
T& tuple_get_tp_(Tp& t, cte_t<false>, cte_t<0>);

template<typename T, typename Tp, std::size_t N>
auto& tuple_get_tp_(Tp& t, cte_t<false>, cte_t<N>) {
    const bool same = is_same_tp<T, typename std::remove_cv<typename std::tuple_element<N-1, Tp>::type>::type>::value;
    static_assert(!(N == 1 && !same), "no such type T in this tuple");
    return tuple_get_tp_<T>(t, cte_t<same>(), cte_t<N-1>());
}

template<typename T, typename ... Args>
auto& tuple_get_tp(std::tuple<Args...>& t) {
    return tuple_get_tp_<T>(t, cte_t<false>(), cte_t<sizeof...(Args)>());
}

template<typename T, typename ... Args>
const auto& tuple_get_tp(const std::tuple<Args...>& t) {
    return tuple_get_tp_<T>(t, cte_t<false>(), cte_t<sizeof...(Args)>());
}

// Merge two tuples together, with matching types only
template<typename ... Args1, typename ... Args2>
void tuple_merge_tp_(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2, cte_t<0>) {
    using type = typename std::tuple_element<0, std::tuple<Args2...>>::type;
    tuple_get_tp<type>(t1) = std::move(std::get<0>(t2));
}

template<typename ... Args1, typename ... Args2, std::size_t N>
void tuple_merge_tp_(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2, cte_t<N>) {
    using type = typename std::tuple_element<N, std::tuple<Args2...>>::type;
    tuple_get_tp<type>(t1) = std::move(std::get<N>(t2));
    tuple_merge_tp_(t1, std::move(t2), cte_t<N-1>());
}

template<typename ... Args1, typename ... Args2>
void tuple_merge_tp(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2) {
    tuple_merge_tp_(t1, std::move(t2), cte_t<sizeof...(Args2)-1>());
}

// Check if a set of types is included into another set of types
template<typename ... Args>
struct check_any_of {
    static constexpr bool check(const std::tuple<>&) {
        return true;
    }

    template<typename T, typename ... Args2>
    static constexpr bool check(const std::tuple<T, Args2...>&) {
        static_assert(is_any_of_tp<T, Args...>::value, "unknown keyword T");
        return check(std::tuple<Args2...>());
    }
};

template<typename T, typename CRTP>
struct keyword_t;

// Class holding all the keywords (wrapper around std::tuple)
template<typename ... Args>
struct keyword_container {
    std::tuple<Args...> vals;

    template<typename ... Keys>
    keyword_container(std::tuple<Keys...>&& kw) {
        check_any_of<Args...>::check(kw);
        tuple_merge_tp(vals, std::move(kw));
    }

    template<typename ... Keys>
    keyword_container(keyword_container<Keys...>&& kw) {
        check_any_of<Args...>::check(kw.vals);
        tuple_merge_tp(vals, std::move(kw.vals));
    }

    template<typename U, typename CRTP2>
    keyword_container<Args..., keyword_t<U,CRTP2>> operator , (keyword_t<U,CRTP2>&& k) {
        return {std::tuple_cat(std::move(vals), std::make_tuple(std::move(k)))};
    }
};

// Helper function to get a keyword from a keyword_container
template<typename T, typename Tp>
auto& keyword_t_(Tp& kw) {
    return tuple_get_tp<T>(kw.vals);
}

// Helper function to create a keyword list
template<typename ... Args>
auto keywords(Args&& ... args) {
    return std::make_tuple(std::forward<Args>(args)...);
}

template<typename T, typename CRTP>
struct keyword_t {
    T value;
    bool provided = false;

    keyword_t() = default;
    keyword_t(keyword_t&&) = default;
    keyword_t& operator= (keyword_t&&) = default;

    template<typename U>
    keyword_t(keyword_t<U,CRTP>&& u) : value(u.value), provided(u.provided) {}

    template<typename U>
    keyword_t(keyword_t<U&,CRTP>&& u) : value(*u.value), provided(u.provided) {}

    T& get() { return value; }
    const T& get() const { return value; }

    template<typename U, typename CRTP2>
    keyword_container<keyword_t, keyword_t<U,CRTP2>> operator , (keyword_t<U,CRTP2>&& k) {
        return {std::make_tuple(std::move(*this), std::move(k))};
    }

    void set(T&& v) {
        provided = true;
        value = std::move(v);
    }
};

template<typename T, typename CRTP>
struct keyword_t<T&, CRTP> {
    T* value;
    bool provided = false;

    keyword_t() = default;
    keyword_t(keyword_t&&) = default;
    keyword_t& operator= (keyword_t&&) = default;

    template<typename U>
    keyword_t(keyword_t<U&,CRTP>&& u) : value(u.value), provided(u.provided) {}

    T& get() { return *value; }
    const T& get() const { return *value; }

    template<typename U, typename CRTP2>
    keyword_container<keyword_t, keyword_t<U,CRTP2>> operator , (keyword_t<U,CRTP2>&& k) {
        return {std::make_tuple(std::move(*this), std::move(k))};
    }

    void set(T& v) {
        provided = true;
        value = &v;
    }
};

#define DECORATE_KW_TYPE(name) name##_t

// Macro to define a new keyword type
#define define_keyword(name) \
    struct DECORATE_KW_TYPE(name); \
    template<typename T> \
    auto name(T&& v) { \
        keyword_t<T, DECORATE_KW_TYPE(name)> r; r.set(std::forward<T>(v)); \
        return r; \
    } \
    auto name(const char* v) { \
        keyword_t<std::string, DECORATE_KW_TYPE(name)> r; r.set(v); \
        return r; \
    } \
    template<std::size_t N> \
    auto name(const char v[N]) { \
        keyword_t<std::string, DECORATE_KW_TYPE(name)> r; r.set(v); \
        return r; \
    } \
    template<typename T> \
    auto name() { \
        return keyword_t<T, DECORATE_KW_TYPE(name)>(); \
    }

// Helper type to construct a keyword_container from a std::tuple
template<typename T>
struct make_kw_container;

template<typename ... Args>
struct make_kw_container<std::tuple<Args...>> {
    using type = keyword_container<Args...>;
};

// Macro to declare the keywords of a function
#define declare_keywords(...) \
    typename make_kw_container<decltype(keywords(__VA_ARGS__))>::type kw__ = keywords(__VA_ARGS__)

// Macro to get a keyword's value inside the function
#define get_keyword(name) \
    keyword_t_<keyword_t<bool, DECORATE_KW_TYPE(name)>>(kw__).get()

// Macro to check if a keyword has been defined
#define keyword_set(name) \
    keyword_t_<keyword_t<bool, DECORATE_KW_TYPE(name)>>(kw__).provided

// Defined keyword list
// Note: since it is not possible to define twice the same keyword, it is better to put them all in
// one place, such as this one.
define_keyword(_btit);
define_keyword(_ltit);
define_keyword(_nth);
define_keyword(_renorm);
define_keyword(_rtit);
define_keyword(_self);
define_keyword(_save);
define_keyword(_thread);
define_keyword(_ttit);
define_keyword(_verbose);
define_keyword(_xtit);
define_keyword(_ytit);

#endif
