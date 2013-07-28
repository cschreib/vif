#ifndef VARIADIC_HPP
#define VARIADIC_HPP

#include <type_traits>
#include <array>
#include <utility>
#include <tuple>
#include <cassert>

// Placeholder variable & type.
static struct placeholder_t {} _;

// Class holding an integer value.
template<std::size_t I>
struct cte_t {};

// Return a particular element of a variadic list
template<std::size_t N, typename T, typename ... Args>
T get_(cte_t<N>, cte_t<N>, T&& t, Args&& ... args) {
    return t;
}

template<std::size_t I, std::size_t N, typename T, typename ... Args>
auto get_(cte_t<I>, cte_t<N>, T&& t, Args&& ... args) -> decltype(get_(cte_t<I+1>(), cte_t<N>(), std::forward<Args>(args)...)) {
    return get_(cte_t<I+1>(), cte_t<N>(), std::forward<Args>(args)...);
}

template<std::size_t N, typename ... Args>
auto get(Args&& ... args) -> decltype(get_(cte_t<1>(), cte_t<N>(), std::forward<Args>(args)...)) {
    return get_(cte_t<0>(), cte_t<N>(), std::forward<Args>(args)...);
}

// Shortcut for std::decay
template<typename T>
using decay_t = typename std::decay<T>::type;

// Class holding a sequence of integer values.
template<std::size_t ...>
struct seq_t {};

// Generate a sequence of integer values from N to D (inclusive).
template<std::size_t N, std::size_t D, std::size_t ... S>
struct gen_seq : public gen_seq<N+1, D, S..., N> {};

template<std::size_t D, std::size_t ... S>
struct gen_seq<D, D, S...> {
  using type = seq_t<S..., D>;
};

// Generate a sequence of integer from 0 to D (exclusive).
template<std::size_t D>
using gen_seq_t = typename gen_seq<0, D-1>::type;

// Check if all the types in a variadic list are of the same given type.
template<typename U, typename T, typename ... Args>
struct are_same {
    static const bool value = std::is_same<U,T>::value ? are_same<U,Args...>::value : false;
};

template<typename U, typename T>
struct are_same<U, T> : public std::is_same<U,T> {};

// Count the number of types in a variadic list that match a given type.
template<typename U, typename T, typename ... Args>
struct count_same {
    static const std::size_t value = std::is_same<U,T>::value + count_same<U,Args...>::value;
};

template<typename U, typename T>
struct count_same<U, T> {
    static const std::size_t value = std::is_same<U,T>::value;
};

// Multiply together an arbitrary list of values.
template<typename U>
U multiply_() {
    return 1;
}

template<typename U, typename T, typename ... Args>
U multiply_(const T& t, Args&& ... args) {
    return t*multiply_<U>(std::forward<Args>(args)...);
}

template<typename T, typename ... Args>
T multiply(const T& t, Args&& ... args) {
    return multiply_<T>(t, std::forward<Args>(args)...);
}

// Check if a type is an std::array
template<typename T>
struct is_array : std::false_type {};
template<typename T, std::size_t N>
struct is_array<std::array<T,N>> : std::true_type {};

// Assign an arbitrary list of values to an array.
template<template<typename> class V, typename T, std::size_t I, typename U, typename ... Args>
void set_array_(V<T>& v, U&& t, Args&& ... args) {
    v[I] = std::move(t);
    set_array_<V,T,I+1>(v, std::forward<Args>(args)...);
}

template<template<typename> class V, typename T, std::size_t I, typename U>
void set_array_(V<T>& v, U&& t) {
    v[I] = std::move(t);
}

template<template<typename> class V, typename T, typename ... Args>
void set_array(V<T>& v, Args&& ... args) {
    assert(v.size() >= sizeof...(Args));
    set_array_<V,T,0>(v, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, std::size_t I, typename U>
void set_array_(std::array<T,N>& v, U&& t) {
    v[I] = static_cast<T>(std::move(t));
}

template<std::size_t N, typename T, std::size_t I, typename U, typename ... Args>
void set_array_(std::array<T,N>& v, U&& t, Args&& ... args) {
    v[I] = static_cast<T>(std::move(t));
    set_array_<N,T,I+1>(v, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, typename ... Args>
void set_array(std::array<T,N>& v, Args&& ... args) {
    static_assert(N >= sizeof...(Args), "too many elements for this array");
    set_array_<N,T,0>(v, std::forward<Args>(args)...);
}

template<typename T>
struct return_type;

template<typename R, typename ... Args>
struct return_type<R (*)(Args...)> {
    using type = R;
};

// Check if the provided type T matches any type in the list
template<typename T, typename U, typename ... Args>
struct is_any_of {
    static const bool value = std::is_same<T,U>::value || is_any_of<T,Args...>::value;
};

template<typename T, typename U>
struct is_any_of<T,U> {
    static const bool value = std::is_same<T,U>::value;
};

// Get a tuple element by type (only available in C++14)
template<typename T, typename ... Args, std::size_t N>
auto& tuple_get_(std::tuple<Args...>& t, cte_t<true>, cte_t<N>) {
    return std::get<N>(t);
}

template<typename T, typename ... Args>
auto& tuple_get_(std::tuple<Args...>& t, cte_t<false>, cte_t<0>);

template<typename T, typename ... Args, std::size_t N>
auto& tuple_get_(std::tuple<Args...>& t, cte_t<false>, cte_t<N>) {
    const bool same = std::is_same<T, typename std::tuple_element<N-1, std::tuple<Args...>>::type>::value;
    static_assert(!(N == 1 && !same), "no such type T in this tuple");
    return tuple_get_<T>(t, cte_t<same>(), cte_t<N-1>());
}

template<typename T, typename ... Args>
auto& tuple_get(std::tuple<Args...>& t) {
    return tuple_get_<T>(t, cte_t<false>(), cte_t<sizeof...(Args)>());
}

// Merge two tuples by matching their types
template<typename ... Args1, typename ... Args2>
void tuple_merge_(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2, cte_t<0>) {
    using type = typename std::tuple_element<0, std::tuple<Args2...>>::type;
    tuple_get<type>(t1) = std::move(std::get<0>(t2));
}

template<typename ... Args1, typename ... Args2, std::size_t N>
void tuple_merge_(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2, cte_t<N>) {
    using type = typename std::tuple_element<N, std::tuple<Args2...>>::type;
    tuple_get<type>(t1) = std::move(std::get<N>(t2));
    tuple_merge_(t1, std::move(t2), cte_t<N-1>());
}

template<typename ... Args1, typename ... Args2>
void tuple_merge(std::tuple<Args1...>& t1, std::tuple<Args2...>&& t2) {
    tuple_merge_(t1, std::move(t2), cte_t<sizeof...(Args2)-1>());
}

// Iterate over the content of a tuple
template<typename T, typename F, std::size_t ... S>
void tuple_for__(T& t, F&& f, seq_t<S...>) {
    auto l = { (f(std::get<S>(t)), 0)... };
}

template<typename T, typename F>
void tuple_for_(T& t, F&& tf, cte_t<0>) {}

template<typename T, typename F, std::size_t N>
void tuple_for_(T& t, F&& tf, cte_t<N>) {
    F& f = tf;
    tuple_for__(t, f, gen_seq_t<N>());
}

template<typename T, typename F>
void tuple_for(T& t, F&& tf) {
    F& f = tf;
    tuple_for_(t, f, cte_t<std::tuple_size<T>::value>());
}

#endif

