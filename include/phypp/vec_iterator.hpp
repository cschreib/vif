#include "phypp/iterator.hpp"

// Iterator for bool vectors (necessary to avoid std::vector<bool>)
template<typename I>
struct bool_iterator_policy {
    static auto get_obj(I& i) -> decltype(dref<bool>(*i)) {
        return dref<bool>(*i);
    }
    static auto get_ptr(I& i) -> decltype(ptr<bool>(*i)) {
        return ptr<bool>(*i);
    }
};

template<typename T, typename C>
using bool_iterator_base = iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using const_bool_iterator_base = const_iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using reverse_bool_iterator_base = reverse_iterator_base<T,C,bool_iterator_policy<T>>;
template<typename T, typename C>
using const_reverse_bool_iterator_base = const_reverse_iterator_base<T,C,bool_iterator_policy<T>>;

namespace std {
    template<typename T, typename C>
    struct iterator_traits<bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = const bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<reverse_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_reverse_bool_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = const bool;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };
}

// Helper to define the vector iterator types.
// Necessary because of tricks to avoid std::vector<bool>.
template<typename T>
struct vec_iterator_type {
    using vtype = typename T::vtype;
    using iterator = typename vtype::iterator;
    using const_iterator = typename vtype::const_iterator;
    using reverse_iterator = typename vtype::reverse_iterator;
    using const_reverse_iterator = typename vtype::const_reverse_iterator;
};

template<std::size_t Dim>
struct vec_iterator_type<vec<Dim,bool>> {
    using vtype = typename vec<Dim,bool>::vtype;
    using iterator = bool_iterator_base<typename vtype::iterator,vec<Dim,bool>>;
    using const_iterator = const_bool_iterator_base<typename vtype::const_iterator,vec<Dim,bool>>;
    using reverse_iterator = reverse_bool_iterator_base<
        typename vtype::reverse_iterator,vec<Dim,bool>>;
    using const_reverse_iterator = const_reverse_bool_iterator_base<
        typename vtype::const_reverse_iterator,vec<Dim,bool>>;
};

// Iterators for vector view
template<typename I>
struct ptr_iterator_policy {
    static auto get_obj(I& i) -> decltype(**i) {
        return **i;
    }
    static auto get_ptr(I& i) -> decltype(*i) {
        return *i;
    }
};


template<typename T, typename C>
using ptr_iterator_base = iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using const_ptr_iterator_base = const_iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using reverse_ptr_iterator_base = reverse_iterator_base<T,C,ptr_iterator_policy<T>>;
template<typename T, typename C>
using const_reverse_ptr_iterator_base = const_reverse_iterator_base<T,C,ptr_iterator_policy<T>>;

namespace std {
    template<typename T, typename C>
    struct iterator_traits<ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<reverse_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C>
    struct iterator_traits<const_reverse_ptr_iterator_base<T,C>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::remove_pointer<typename iterator_traits<T>::value_type>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };
}
