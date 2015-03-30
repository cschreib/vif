#ifndef VARIADIC_HPP
#define VARIADIC_HPP

#include <type_traits>
#include <array>
#include <utility>
#include <cassert>

// Class holding an integer value.
template<std::size_t I>
struct cte_t {};

// Return a particular element of a variadic list
template<std::size_t N, typename T, typename ... Args>
auto get_(cte_t<N>, cte_t<N>, T&& t, Args&& ... args) -> decltype(std::forward<T>(t)) {
    return std::forward<T>(t);
}

template<std::size_t I, std::size_t N, typename T, typename ... Args>
auto get_(cte_t<I>, cte_t<N>, T&& t, Args&& ... args) ->
    decltype(get_(cte_t<I+1>(), cte_t<N>(), std::forward<Args>(args)...)) {
    return get_(cte_t<I+1>(), cte_t<N>(), std::forward<Args>(args)...);
}

template<std::size_t N, typename ... Args>
auto get(Args&& ... args) -> decltype(get_(cte_t<0>(), cte_t<N>(), std::forward<Args>(args)...)) {
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

// Simple type list
template<typename ... Args>
struct type_list {};

template<typename T, typename ... Args>
type_list<Args...> pop_front(type_list<T,Args...>) { return type_list<Args...>{}; }

template<typename T, typename ... Args>
type_list<Args...> pop_back(type_list<Args...,T>) { return type_list<Args...>{}; }

// Generate a type indicating that its value is to be repeated N times
template<std::size_t N, typename T>
struct repeated_value {
    T value;
};

template<std::size_t N, typename T>
repeated_value<N,T> repeat(T t) {
    return repeated_value<N,T>{t};
}

// Unwrap repeated_value arguments into individual arguments
template<typename ... IArgs>
struct unfolder;

template<typename F, typename ... Args>
struct deduce_return_type {
    using type = decltype(std::declval<typename std::decay<F>::type>()(std::declval<Args&>()...));
};

template<>
struct unfolder<> {
    template<typename F, typename ... OArgs>
    using out_type = typename deduce_return_type<F,OArgs...>::type;

    template<typename ... OArgs>
    struct out {
        template<typename F>
        static out_type<F,OArgs...> unfold(F&& func, OArgs& ... oargs) {
            return func(oargs...);
        }
    };
};

template<typename T, typename ... IArgs>
struct unfolder<T, IArgs...> {
    template<typename F, typename ... OArgs>
    using out_type = typename unfolder<IArgs...>::template out_type<F, OArgs..., T>;

    template<typename ... OArgs>
    struct out {
        template<typename F>
        static out_type<F,OArgs...> unfold(F&& func, T& t, IArgs& ... iargs, OArgs& ... oargs) {
            using actor = typename unfolder<IArgs...>::template out<OArgs..., T>;
            return actor::unfold(std::forward<F>(func), iargs..., oargs..., t);
        }
    };
};

template<std::size_t N, typename T, typename ... IArgs>
struct unfolder<repeated_value<N,T>, IArgs...> {
    template<typename F, std::size_t I, typename ... OArgs>
    struct out_type_impl;

    template<typename F, typename ... OArgs>
    struct out_type_impl<F, 0, OArgs...> {
        using type = typename unfolder<IArgs...>::template out_type<F, OArgs...>;
    };

    template<typename F, std::size_t I, typename ... OArgs>
    struct out_type_impl {
        using type = typename out_type_impl<F, I-1, OArgs..., T>::type;
    };

    template<typename F, typename ... OArgs>
    using out_type = typename out_type_impl<F, N, OArgs...>::type;

    template<typename ... OArgs>
    struct out {
        template<typename F, typename U>
        static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, cte_t<0>, IArgs& ... iargs,
            OArgs& ... oargs) {

            using actor = typename unfolder<IArgs...>::template out<OArgs...>;
            return actor::unfold(std::forward<F>(func), iargs..., oargs...);
        }

        template<typename F, std::size_t I, typename U>
        static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, cte_t<I>, IArgs& ... iargs,
            OArgs& ... oargs) {

            using actor = typename unfolder::template out<OArgs..., T>;
            return actor::unfold_repeat(std::forward<F>(func), t, cte_t<I-1>(),
                iargs..., oargs..., t.value);
        }

        template<typename F, typename U>
        static out_type<F,OArgs...> unfold(F&& func, U& t, IArgs& ... iargs, OArgs& ... oargs) {
            return unfold_repeat(std::forward<F>(func), t, cte_t<N>(), iargs..., oargs...);
        }
    };
};

template<std::size_t N, typename T, typename ... IArgs>
struct unfolder<const repeated_value<N,T>, IArgs...> {
    template<typename F, std::size_t I, typename ... OArgs>
    struct out_type_impl;

    template<typename F, typename ... OArgs>
    struct out_type_impl<F, 0, OArgs...> {
        using type = typename unfolder<IArgs...>::template out_type<F, OArgs...>;
    };

    template<typename F, std::size_t I, typename ... OArgs>
    struct out_type_impl {
        using type = typename out_type_impl<F, I-1, OArgs..., const T>::type;
    };

    template<typename F, typename ... OArgs>
    using out_type = typename out_type_impl<F, N, OArgs...>::type;

    template<typename ... OArgs>
    struct out {
        template<typename F, typename U>
        static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, cte_t<0>, IArgs& ... iargs,
            OArgs& ... oargs) {

            using actor = typename unfolder<IArgs...>::template out<OArgs...>;
            return actor::unfold(std::forward<F>(func), iargs..., oargs...);
        }

        template<typename F, std::size_t I, typename U>
        static out_type<F,OArgs...> unfold_repeat(F&& func, U& t, cte_t<I>, IArgs& ... iargs,
            OArgs& ... oargs) {

            using actor = typename unfolder::template out<OArgs..., const T>;
            return actor::unfold_repeat(std::forward<F>(func), t, cte_t<I-1>(),
                iargs..., oargs..., t.value);
        }

        template<typename F, typename U>
        static out_type<F,OArgs...> unfold(F&& func, U& t, IArgs& ... iargs, OArgs& ... oargs) {
            return unfold_repeat(std::forward<F>(func), t, cte_t<N>(), iargs..., oargs...);
        }
    };
};

template<typename F, typename ... Args>
typename unfolder<Args...>::template out_type<F> unfold(F&& func, Args& ... args) {
    using actor = typename unfolder<Args...>::template out<>;
    return actor::unfold(std::forward<F>(func), args...);
}

// Check if all the types in a variadic list are of the same given type.
template<typename U, typename T, typename ... Args>
struct are_same :
    std::integral_constant<bool, std::is_same<U,T>::value ? are_same<U,Args...>::value : false> {};

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

// Get the size of an std::array
template<typename T>
struct array_size;

template<std::size_t N, typename T>
struct array_size<std::array<T,N>> {
    static const std::size_t size = N;
};

// Assign an arbitrary list of values to an array.
template<typename T>
struct add_dim : std::integral_constant<std::size_t, 1> {};
template<std::size_t N, typename T>
struct add_dim<std::array<T,N>> : std::integral_constant<std::size_t, N> {};

template<typename T, typename ... Args>
struct dim_total : std::integral_constant<std::size_t,
    add_dim<typename std::decay<T>::type>::value + dim_total<Args...>::value> {};

template<typename T>
struct dim_total<T> : std::integral_constant<std::size_t,
    add_dim<typename std::decay<T>::type>::value> {};

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

template<std::size_t N, typename T, std::size_t I>
void set_array_(std::array<T,N>& v, cte_t<I>) {}

template<std::size_t N, typename T, std::size_t I, typename U, typename ... Args>
void set_array_(std::array<T,N>& v, cte_t<I>, const U& t, Args&& ... args) {
    v[I] = static_cast<T>(t);
    set_array_(v, cte_t<I+1>{}, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, std::size_t I, typename U, std::size_t M, typename ... Args>
void set_array_(std::array<T,N>& v, cte_t<I>, const std::array<U,M>& t, Args&& ... args) {
    for (std::size_t i = 0; i < M; ++i) {
        v[I+i] = static_cast<T>(t[i]);
    }

    set_array_(v, cte_t<I+M>{}, std::forward<Args>(args)...);
}

template<std::size_t N, typename T, typename ... Args>
void set_array(std::array<T,N>& v, Args&& ... args) {
    static_assert(N == dim_total<Args...>::value, "wrong number of elements for this array");
    set_array_(v, cte_t<0>{}, std::forward<Args>(args)...);
}

template<typename T>
struct return_type {
    using type = typename return_type<decltype(&T::operator())>::type;
};

template<typename R, typename ... Args>
struct return_type<R (*)(Args...)> {
    using type = R;
};

template<typename R, typename T, typename ... Args>
struct return_type<R (T::*)(Args...)> {
    using type = R;
};

template<typename R, typename T, typename ... Args>
struct return_type<R (T::*)(Args...) const> {
    using type = R;
};

// Check if the provided type T matches any type in the list
template<typename T, typename U, typename ... Args>
struct is_any_type_of :
    std::integral_constant<bool, std::is_same<T,U>::value || is_any_type_of<T,Args...>::value> {};

template<typename T, typename U>
struct is_any_type_of<T,U> : std::is_same<T,U> {};

// Get the Nth type from a variadic list
template<std::size_t N, std::size_t I, typename ... Args>
struct nth_type_impl;

template<std::size_t N, std::size_t I>
struct nth_type_impl<N,I> {
    static_assert(true, "not enough type in this type list");
};

template<std::size_t I, typename T, typename ... Args>
struct nth_type_impl<I,I,T,Args...> {
    using type = T;
};

template<std::size_t N, std::size_t I, typename T, typename ... Args>
struct nth_type_impl<N,I,T,Args...> : nth_type_impl<N-1,I, Args...> {};

template<std::size_t N, typename ... Args>
using nth_type = typename nth_type_impl<sizeof...(Args), N+1, Args...>::type;

template<typename T, typename C>
using constify = typename std::conditional<std::is_const<C>::value, const T, T>::type;

#endif

