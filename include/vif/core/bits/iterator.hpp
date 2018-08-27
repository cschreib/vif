#ifndef VIF_INCLUDING_CORE_VEC_BITS
#error this file is not meant to be included separately, include "vif/core/vec.hpp" instead
#endif

namespace vif {
namespace impl {
    // Iterator for normal vectors.
    // Note: The default implementation uses native iterators, so these functions are
    // actually not used in this case! However, these functions can be (and are) used
    // to create your own iterator type, if you do not want to repeat too much code.
    // See below for a few examples.
    struct native_iterator_policy {
        template<typename I>
        void initialize(I& i) const {}
        template<typename I>
        auto get_obj(I& i) const -> decltype(*i) {
            return *i;
        }
        template<typename I>
        auto get_ptr(I& i) const -> decltype(&(*i)) {
            return &(*i);
        }
        template<typename I>
        void increment(I& i) const {
            ++i;
        }
        template<typename I>
        void decrement(I& i) const {
            --i;
        }
        template<typename I, typename U>
        void advance(I& i, U n) const {
            i += n;
        }
        template<typename I>
        auto difference(const I& i, const I& j) const -> decltype(i-j) {
            return i-j;
        }
    };

    // Iterator for bool vectors (necessary to avoid std::vector<bool>)
    struct bool_iterator_policy : native_iterator_policy {
        template<typename I>
        auto get_obj(I& i) const -> decltype(impl::dref<bool>(*i)) {
            return impl::dref<bool>(*i);
        }
        template<typename I>
        auto get_ptr(I& i) const -> decltype(impl::ptr<bool>(*i)) {
            return impl::ptr<bool>(*i);
        }
    };

    // Iterators for vector view
    struct ptr_iterator_policy : native_iterator_policy {
        template<typename I>
        auto get_obj(I& i) const -> decltype(**i) {
            return **i;
        }
        template<typename I>
        auto get_ptr(I& i) const -> decltype(*i) {
            return *i;
        }
    };

    // Default policy for each vector type
    template<typename T>
    struct default_iterator_policy_t {
        using policy = native_iterator_policy;
    };

    template<std::size_t Dim>
    struct default_iterator_policy_t<vec<Dim,bool>> {
        using policy = bool_iterator_policy;
    };

    template<std::size_t Dim, typename T>
    struct default_iterator_policy_t<vec<Dim,T*>> {
        using policy = ptr_iterator_policy;
    };

    template<typename T>
    using default_iterator_policy = typename default_iterator_policy_t<T>::policy;

    // Iterators for strided traversal to iterate on one dimension
    template<typename T>
    struct strided_iterator_policy : default_iterator_policy<T> {
        int_t stride = 1;
        int_t offset = 0;

        strided_iterator_policy() = default;
        strided_iterator_policy(int_t s, int_t o) : stride(s), offset(o) {}

        template<typename I>
        void initialize(I& i) const {
            i += offset;
        }
        template<typename I, typename U>
        void advance(I& i, U n) const {
            i += n*stride;
        }
        template<typename I>
        void increment(I& i) const {
            advance(i, 1);
        }
        template<typename I>
        void decrement(I& i) const {
            advance(i, -1);
        }
        template<typename I>
        auto difference(const I& i, const I& j) const -> decltype(i-j) {
            return (i - j)/stride;
        }
    };

    // Iterators for strided traversal to iterate on one dimension (raw data)
    template<typename T>
    struct raw_strided_iterator_policy : native_iterator_policy {
        int_t stride = 1;
        int_t offset = 0;

        raw_strided_iterator_policy() = default;
        raw_strided_iterator_policy(int_t s, int_t o) : stride(s), offset(o) {}

        template<typename I>
        void initialize(I& i) const {
            i += offset;
        }
        template<typename I, typename U>
        void advance(I& i, U n) const {
            i += n*stride;
        }
        template<typename I>
        void increment(I& i) const {
            advance(i, 1);
        }
        template<typename I>
        void decrement(I& i) const {
            advance(i, -1);
        }
        template<typename I>
        auto difference(const I& i, const I& j) const -> decltype(i-j) {
            return (i - j)/stride;
        }
    };

    // Helper to define the vector iterator types.
    template<typename T, typename P>
    struct vec_iterator_type {
        using vtype = typename T::vtype;
        using iterator = iterator_adaptor<typename vtype::iterator,T,P>;
        using const_iterator = const_iterator_adaptor<typename vtype::const_iterator,T,P>;
        using reverse_iterator = reverse_iterator_adaptor<
            typename vtype::reverse_iterator,T,P>;
        using const_reverse_iterator = const_reverse_iterator_adaptor<
            typename vtype::const_reverse_iterator,T,P>;
    };

    // Default implementation uses native iterators for maximal speed.
    template<typename T>
    struct vec_iterator_type<T,native_iterator_policy> {
        using vtype = typename T::vtype;
        using iterator = typename vtype::iterator;
        using const_iterator = typename vtype::const_iterator;
        using reverse_iterator = typename vtype::reverse_iterator;
        using const_reverse_iterator = typename vtype::const_reverse_iterator;
    };
}
}

namespace std {
    template<typename T, typename C, typename P>
    struct iterator_traits<vif::impl::iterator_adaptor<T,C,P>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::decay<decltype(std::declval<P>().get_obj(std::declval<T&>()))>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C, typename P>
    struct iterator_traits<vif::impl::const_iterator_adaptor<T,C,P>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::decay<decltype(std::declval<P>().get_obj(std::declval<T&>()))>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C, typename P>
    struct iterator_traits<vif::impl::reverse_iterator_adaptor<T,C,P>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::decay<decltype(std::declval<P>().get_obj(std::declval<T&>()))>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };

    template<typename T, typename C, typename P>
    struct iterator_traits<vif::impl::const_reverse_iterator_adaptor<T,C,P>> {
        using difference_type = typename iterator_traits<T>::difference_type;
        using value_type = typename std::decay<decltype(std::declval<P>().get_obj(std::declval<T&>()))>::type;
        using pointer = value_type*;
        using reference = value_type&;
        using iterator_category = typename iterator_traits<T>::iterator_category;
    };
}
