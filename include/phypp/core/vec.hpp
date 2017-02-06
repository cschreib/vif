#ifndef PHYPP_CORE_VEC_HPP
#define PHYPP_CORE_VEC_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <initializer_list>
#include "phypp/core/typedefs.hpp"
#include "phypp/core/range.hpp"
#include "phypp/core/meta.hpp"
#include "phypp/core/error.hpp"
#include "phypp/core/iterator_base.hpp"

namespace phypp {
    namespace impl {
        // Tag type to mark initialization of a reference vector.
        static struct vec_ref_tag_t {}    vec_ref_tag;
        // Tag type to permit copy initialization without data copy.
        static struct vec_nocopy_tag_t {} vec_nocopy_tag;
    }
}

// Helper code is located in separate headers for clarity
#define PHYPP_INCLUDING_CORE_VEC_BITS
#include "phypp/core/bits/helpers.hpp"
#include "phypp/core/bits/iterator.hpp"
#include "phypp/core/bits/access.hpp"
#include "phypp/core/bits/initializer_list.hpp"
#undef PHYPP_INCLUDING_CORE_VEC_BITS

namespace phypp {
    ////////////////////////////////////////////
    //            Generic vector              //
    ////////////////////////////////////////////

    template<std::size_t Dim, typename Type>
    struct vec;

    // Disable dimension zero
    template<typename Type>
    struct vec<0,Type> {};

    // The vector
    template<std::size_t Dim, typename Type>
    struct vec {
        static_assert(!std::is_pointer<Type>::value, "library bug: pointer specialization failed");

        using effective_type = vec;
        using rtype = meta::rtype_t<Type>;
        using dtype = meta::dtype_t<Type>;
        using drtype = meta::dtype_t<Type>;
        using vtype = std::vector<dtype>;
        using dim_type = std::array<std::size_t, Dim>;
        struct comparator {
            constexpr bool operator() (const dtype& t1, const dtype& t2) const {
                return t1 < t2;
            }
            template<typename U>
            constexpr bool operator() (const dtype& t1, const U& t2) const {
                return t1 < t2;
            }
            template<typename U>
            constexpr bool operator() (const U& t1, const dtype& t2) const {
                return t1 < t2;
            }
        };

        vtype    data;
        dim_type dims = {{0}};

        // Default constructor
        vec() : safe(*this) {}

        // Copy constructor
        vec(const vec& v) : data(v.data), dims(v.dims), safe(*this) {}
        // Copy constructor (only dims)
        vec(impl::vec_nocopy_tag_t, const vec& v) : dims(v.dims), safe(*this) {}

        // Move constructor
        vec(vec&& v) : data(std::move(v.data)), dims(v.dims), safe(*this) {
            for (uint_t i : range(Dim)) {
                v.dims[i] = 0;
            }
        }

        // Dimension constructor
        template<typename ... Args, typename enable =
            typename std::enable_if<meta::is_dim_list<Args...>::value>::type>
        explicit vec(Args&& ... d) : safe(*this) {
            static_assert(meta::dim_total<Args...>::value == Dim, "dimension list does not match "
                "the dimensions of this vector");

            impl::set_array(dims, std::forward<Args>(d)...);
            resize();
        }

        // Initializer list constructor
        vec(meta::nested_initializer_list<Dim,meta::dtype_t<Type>> il) : safe(*this) {
            impl::vec_ilist::helper<Dim, Type>::fill(*this, il);
        }

        // Implicit conversion
        template<typename T, typename std::enable_if<!meta::vec_only_explicit_convertible<meta::rtype_t<T>,Type>::value, bool>::type = false>
        vec(const vec<Dim,T>& v) : dims(v.dims), safe(*this) {
            static_assert(meta::vec_implicit_convertible<meta::rtype_t<T>,Type>::value,
                "could not construct vector from non-convertible type");
            data.resize(v.data.size());
            for (uint_t i : range(v)) {
                data[i] = v.safe[i];
            }
        }

        // Explicit conversion
        template<typename T, typename std::enable_if<meta::vec_only_explicit_convertible<meta::rtype_t<T>,Type>::value, bool>::type = false>
        explicit vec(const vec<Dim,T>& v) : dims(v.dims), safe(*this) {
            static_assert(meta::vec_explicit_convertible<meta::rtype_t<T>,Type>::value,
                "could not construct vector from non-convertible type");
            data.resize(v.data.size());
            for (uint_t i : range(v)) {
                data[i] = v.safe[i];
            }
        }

        vec& operator = (meta::nested_initializer_list<Dim,meta::dtype_t<Type>> il) {
            impl::vec_ilist::helper<Dim, Type>::fill(*this, il);
            return *this;
        }

        vec& operator = (const vec& v) {
            data = v.data;
            dims = v.dims;
            return *this;
        }

        vec& operator = (vec&& v) {
            data = std::move(v.data);
            dims = std::move(v.dims);
            return *this;
        }

        template<typename T>
        vec& operator = (const vec<Dim,T>& v) {
            static_assert(meta::vec_implicit_convertible<T,Type>::value, "could not assign vectors of "
                "non-implicitly-convertible types");

            if (v.view_same(*this)) {
                std::vector<dtype> t(v.data.size());
                for (uint_t i : range(v)) {
                    t[i] = v.safe[i];
                }
                dims = v.dims;
                data.resize(v.data.size());
                for (uint_t i : range(v)) {
                    data[i] = t[i];
                }
            } else {
                dims = v.dims;
                data.resize(v.data.size());
                for (uint_t i : range(v)) {
                    data[i] = v.safe[i];
                }
            }

            return *this;
        }

        template<std::size_t D, typename T>
        bool view_same(const vec<D,T>&) const {
            return false;
        }

        template<std::size_t D>
        bool view_same(const vec<D,Type*>& v) const {
            return v.view_same(*this);
        }

        bool empty() const {
            return data.empty();
        }

        uint_t size() const {
            return data.size();
        }

        void resize() {
            uint_t size = 1;
            for (uint_t i : range(Dim)) {
                size *= dims[i];
            }

            data.resize(size);
        }

        template<typename ... Args>
        void resize(Args&& ... d) {
            impl::set_array(dims, d...);
            resize();
        }

        void clear() {
            data.clear();
            for (uint_t i : range(Dim)) {
                dims[i] = 0;
            }
        }

        const Type& back() const {
            static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
            return reinterpret_cast<const Type&>(data.back());
        }

        Type& back() {
            static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
            return const_cast<Type&>(const_cast<const vec&>(*this).back());
        }

        const Type& front() const {
            static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
            return reinterpret_cast<const Type&>(data.front());
        }

        Type& front() {
            static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
            return const_cast<Type&>(const_cast<const vec&>(*this).front());
        }

        void push_back(const Type& t) {
            static_assert(Dim == 1, "cannot call push_back(Type) on multidimensional vectors");
            data.push_back(t);
            ++dims[0];
        }

        void push_back(Type&& t) {
            static_assert(Dim == 1, "cannot call push_back(Type) on multidimensional vectors");
            data.push_back(std::move(t));
            ++dims[0];
        }

        template<typename T = Type, typename enable =
            typename std::enable_if<std::is_convertible<T,Type>::value>::type>
        void push_back(const vec<Dim-1,T>& t) {
            static_assert(Dim > 1, "cannot call push_back(vec<D-1>) on monodimensional vectors");
            if (empty()) {
                dims[0] = 1;
                for (uint_t i = 1; i < Dim; ++i) {
                    dims[i] = t.dims[i-1];
                }

                data.resize(t.data.size());
                std::copy(t.data.begin(), t.data.end(), data.begin());
            } else {
                for (uint_t i = 0; i < Dim-1; ++i) {
                    phypp_check(dims[i+1] == t.dims[i],
                        "push_back: incompatible dimensions (", dims, " vs ", t.dims, ")");
                }

                data.insert(data.end(), t.data.begin(), t.data.end());
                ++dims[0];
            }
        }

        template<typename T = Type, typename enable =
            typename std::enable_if<std::is_convertible<T,Type>::value>::type>
        void push_back(const vec<Dim-1,T*>& t) {
            static_assert(Dim > 1, "cannot call push_back(vec<D-1>) on monodimensional vectors");
            push_back(t.concretise());
        }

        void reserve(uint_t n) {
            data.reserve(n);
        }

        const vec<Dim,Type>& concretise() const {
            return *this;
        }

        template<typename T>
        uint_t to_idx_(T ui, meta::cte_t<false>) const {
            phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
                data.size(), ")");
            return ui;
        }

        template<typename T>
        uint_t to_idx_(T i, meta::cte_t<true>) const {
            if (i < 0) i += data.size();
            phypp_check(i >= 0, "operator[]: index out of bounds (", i+data.size(), " vs. ",
                data.size(), ")");
            uint_t ui(i);
            phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
                data.size(), ")");
            return ui;
        }

        template<std::size_t D, typename T>
        uint_t to_idx_(T ui, meta::cte_t<false>) const {
            phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
                dims[D], ")");
            return ui;
        }

        template<std::size_t D, typename T>
        uint_t to_idx_(T i, meta::cte_t<true>) const {
            if (i < 0) i += dims[D];
            phypp_check(i >= 0, "operator(): index out of bounds (", i+data.size(), " vs. ",
                dims[D], ")");
            uint_t ui(i);
            phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
                dims[D], ")");
            return ui;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        uint_t to_idx(T ix) const {
            return to_idx_(ix, meta::cte_t<std::is_signed<T>::value>());
        }

        template<std::size_t D, typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        uint_t to_idx(T ix) const {
            return to_idx_<D>(ix, meta::cte_t<std::is_signed<T>::value>());
        }

        uint_t pitch(uint_t i) const {
            uint_t p = 1;
            for (uint_t j = i+1; j < Dim; ++j) {
                p *= dims[j];
            }
            return p;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        Type& operator [] (T i) {
            return const_cast<Type&>(const_cast<const vec&>(*this)[i]);
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        const Type& operator [] (T i) const {
            return reinterpret_cast<const Type&>(data[to_idx(i)]);
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T>& i) {
            vec<1,Type*> v(impl::vec_ref_tag, *this);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = impl::ptr<Type>(data[to_idx(i.safe[j])]);
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T*>& i) {
            vec<1,Type*> v(impl::vec_ref_tag, *this);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = impl::ptr<Type>(data[to_idx(i.safe[j])]);
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (std::initializer_list<T> l) {
            return operator[](vec<1,T>{l});
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T>& i) const {
            vec<1,const Type*> v(impl::vec_ref_tag, *this);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = impl::ptr<Type>(data[to_idx(i.safe[j])]);
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
            vec<1,const Type*> v(impl::vec_ref_tag, *this);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = impl::ptr<Type>(data[to_idx(i.safe[j])]);
            }
            return v;
        }

        vec<1,Type*> operator [] (impl::range_impl::full_range_t rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::left_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::right_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::left_right_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (impl::range_impl::full_range_t rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::left_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::right_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::left_right_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) ->
            typename impl::vec_access::helper<true, false, Dim, Type, Args...>::type {
            static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return impl::vec_access::helper<true, false, Dim, Type, Args...>::access(*this, i...);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) const ->
            typename impl::vec_access::helper<true, true, Dim, Type, Args...>::type {
            static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return impl::vec_access::helper<true, true, Dim, Type, Args...>::access(*this, i...);
        }

        struct safe_proxy {
            vec& parent;

            safe_proxy(vec& p) : parent(p) {}
            safe_proxy(const vec&) = delete;
            safe_proxy(vec&&) = delete;
            safe_proxy& operator=(const vec&) = delete;
            safe_proxy& operator=(vec&&) = delete;

            Type& operator [] (uint_t i) {
                return const_cast<Type&>(const_cast<const safe_proxy&>(*this)[i]);
            }

            const Type& operator [] (uint_t i) const {
                return reinterpret_cast<const Type&>(parent.data[i]);
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,Type*> operator [] (const vec<1,T>& i) {
                vec<1,Type*> v(impl::vec_ref_tag, parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = impl::ptr<Type>(parent.data[i.safe[j]]);
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_integral<T>::value>::type>
            vec<1,Type*> operator [] (std::initializer_list<T> l) {
                return operator[](vec<1,T>{l});
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,Type*> operator [] (const vec<1,T*>& i) {
                vec<1,Type*> v(impl::vec_ref_tag, parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = impl::ptr<Type>(parent.data[i.safe[j]]);
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,const Type*> operator [] (const vec<1,T>& i) const {
                vec<1,const Type*> v(impl::vec_ref_tag, parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = impl::ptr<Type>(parent.data[i.safe[j]]);
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
                vec<1,const Type*> v(impl::vec_ref_tag, parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = impl::ptr<Type>(parent.data[i.safe[j]]);
                }
                return v;
            }

            vec<1,Type*> operator [] (impl::range_impl::full_range_t rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::left_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::right_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::left_right_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (impl::range_impl::full_range_t rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::left_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::right_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::left_right_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            template<typename ... Args>
            auto operator () (const Args& ... i) ->
                typename impl::vec_access::helper<false, false, Dim, Type, Args...>::type {
                static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                    "wrong number of indices for this vector");
                return impl::vec_access::helper<false, false, Dim, Type, Args...>::access(parent, i...);
            }

            template<typename ... Args>
            auto operator () (const Args& ... i) const ->
                typename impl::vec_access::helper<false, true, Dim, Type, Args...>::type {
                static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                    "wrong number of indices for this vector");
                return impl::vec_access::helper<false, true, Dim, Type, Args...>::access(parent, i...);
            }
        } safe;

        vec operator + () const {
            return *this;
        }

        vec operator - () const {
            vec v = *this;
            for (auto& t : v) {
                t = -t;
            }

            return v;
        }

        #define OPERATOR(op) \
            template<typename U> \
            vec& operator op (const vec<Dim,U>& u) { \
                phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                    "' (", dims, " vs ", u.dims, ")"); \
                if (u.view_same(*this)) { \
                    std::vector<dtype> t; t.resize(data.size()); \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        t[i] = u.safe[i]; \
                    } \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        data[i] op t[i]; \
                    } \
                } else { \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        data[i] op u.safe[i]; \
                    } \
                } \
                return *this; \
            } \
            template<typename U> \
            vec& operator op (U u) { \
                for (auto& v : data) { \
                    v op u; \
                } \
                return *this; \
            }

        OPERATOR(*=)
        OPERATOR(/=)
        OPERATOR(+=)
        OPERATOR(-=)

        #undef OPERATOR

        using iterator = typename impl::vec_iterator_type<vec,impl::default_iterator_policy<vec>>::iterator;
        using const_iterator = typename impl::vec_iterator_type<vec,impl::default_iterator_policy<vec>>::const_iterator;


        iterator begin() {
            return data.begin();
        }

        iterator end() {
            return data.end();
        }

        const_iterator begin() const {
            return data.begin();
        }

        const_iterator end() const {
            return data.end();
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::iterator begin(Property p = Property()) {
            return typename impl::vec_iterator_type<vec,Property>::iterator(
                data.begin(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::iterator end(Property p = Property()) {
            return typename impl::vec_iterator_type<vec,Property>::iterator(
                data.end(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::const_iterator begin(Property p = Property()) const {
            return typename impl::vec_iterator_type<vec,Property>::const_iterator(
                data.begin(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::const_iterator end(Property p = Property()) const {
            return typename impl::vec_iterator_type<vec,Property>::const_iterator(
                data.end(), p
            );
        }

        template<typename ... Args>
        impl::vec_access::strided_range<vec,impl::strided_iterator_policy<vec>> stride(const Args& ... i) {
            return impl::vec_access::strided_range<vec,impl::strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<const vec,impl::strided_iterator_policy<vec>> stride(const Args& ... i) const {
            return impl::vec_access::strided_range<const vec,impl::strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<vec,impl::raw_strided_iterator_policy<vec>> raw_stride(const Args& ... i) {
            return impl::vec_access::strided_range<vec,impl::raw_strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<const vec,impl::raw_strided_iterator_policy<vec>> raw_stride(const Args& ... i) const {
            return impl::vec_access::strided_range<const vec,impl::raw_strided_iterator_policy<vec>>(*this, i...);
        }
    };


    ////////////////////////////////////////////
    //              Vector view               //
    ////////////////////////////////////////////

    template<std::size_t Dim, typename Type>
    struct vec<Dim,Type*> {
        using rtype = meta::rtype_t<Type*>;
        using effective_type = vec<Dim,rtype>;
        using dtype = Type;
        using drtype = rtype;
        using vtype = std::vector<dtype*>;
        using dim_type = std::array<std::size_t, Dim>;
        struct comparator {
            bool operator() (const dtype* t1, const dtype* t2) {
                return *t1 < *t2;
            }
            template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
            bool operator() (const dtype* t1, const U& t2) {
                return *t1 < t2;
            }
            template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
            bool operator() (const U& t1, const dtype* t2) {
                return t1 < *t2;
            }

            constexpr bool operator() (const dtype& t1, const dtype& t2) const {
                return t1 < t2;
            }
            template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
            constexpr bool operator() (const dtype& t1, const U& t2) const {
                return t1 < t2;
            }
            template<typename U, typename enable = typename std::enable_if<!std::is_pointer<U>::value>::type>
            constexpr bool operator() (const U& t1, const dtype& t2) const {
                return t1 < t2;
            }
        };

        void*    parent = nullptr;
        vtype    data;
        dim_type dims = {{0}};

        vec() = delete;
        vec(const vec& v) : parent(v.parent), data(v.data), dims(v.dims), safe(*this) {}
        vec(impl::vec_nocopy_tag_t, const vec& v) : parent(v.parent), dims(v.dims), safe(*this) {}
        vec(vec&& v) : parent(v.parent), data(std::move(v.data)), dims(std::move(v.dims)), safe(*this) {}

        template<std::size_t D, typename T, typename enable =
            typename std::enable_if<std::is_same<meta::rtype_t<T>, rtype>::value>::type>
        explicit vec(impl::vec_ref_tag_t, const vec<D,T>& p) : safe(*this) {
            parent = static_cast<void*>(const_cast<vec<D,T>*>(&p));
        }

        explicit vec(impl::vec_ref_tag_t, void* p) : safe(*this) {
            parent = p;
        }

        vec& operator = (const Type& t) {
            for (uint_t i = 0; i < data.size(); ++i) {
                *data[i] = t;
            }
            return *this;
        }

        vec& operator = (meta::nested_initializer_list<Dim,meta::dtype_t<Type>> il) {
            impl::vec_ilist::helper<Dim, Type*>::fill(*this, il);
            return *this;
        }

        template<typename T>
        void assign_(const vec<Dim,T>& v) {
            if (view_same(v)) {
                // Make a copy to prevent aliasing
                std::vector<dtype> t(v.data.size());
                for (uint_t i : range(v)) {
                    t[i] = v.safe[i];
                }

                // Actual assignment
                for (uint_t i : range(v)) {
                    *data[i] = t[i];
                }
            } else {
                // No aliasing possible, assign directly
                for (uint_t i : range(v)) {
                    *data[i] = v.safe[i];
                }
            }
        }

        template<typename T>
        vec& operator = (const vec<Dim,T>& v) {
            static_assert(meta::vec_implicit_convertible<T,Type>::value, "could not assign vectors of "
                "non-implicitly-convertible types");
            phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
                v.dims, " to ", dims, ")");

            assign_(v);

            return *this;
        }

        vec& operator = (const vec& v) {
            phypp_check(data.size() == v.data.size(), "incompatible size in assignment (assigning ",
                v.dims, " to ", dims, ")");

            assign_(v);

            return *this;
        }

        template<std::size_t D, typename T>
        bool view_same(const vec<D,T>& v) const {
            return false;
        }

        template<std::size_t D>
        bool view_same(const vec<D,Type*>& v) const {
            return v.parent == parent;
        }

        template<std::size_t D>
        bool view_same(const vec<D,Type>& v) const {
            return static_cast<void*>(const_cast<vec<D,Type>*>(&v)) == parent;
        }

        bool empty() const {
            return data.empty();
        }

        uint_t size() const {
            return data.size();
        }

        void resize() {
            uint_t size = 1;
            for (uint_t i = 0; i < Dim; ++i) {
                size *= dims[i];
            }

            data.resize(size);
        }

        effective_type concretise() const {
            return *this;
        }

        const Type& back() const {
            static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
            return reinterpret_cast<const Type&>(*data.back());
        }

        Type& back() {
            static_assert(Dim == 1, "cannot call back() on multidimensional verctors");
            return const_cast<Type&>(const_cast<const vec&>(*this).back());
        }

        const Type& front() const {
            static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
            return reinterpret_cast<const Type&>(*data.front());
        }

        Type& front() {
            static_assert(Dim == 1, "cannot call front() on multidimensional verctors");
            return const_cast<Type&>(const_cast<const vec&>(*this).front());
        }

        template<typename T>
        T to_idx_(T ui, meta::cte_t<false>) const {
            phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
                data.size(), ")");
            return ui;
        }

        template<typename T>
        uint_t to_idx_(T i, meta::cte_t<true>) const {
            if (i < 0) i += data.size();
            phypp_check(i >= 0, "operator[]: index out of bounds (", i+data.size(), " vs. ",
                data.size(), ")");
            uint_t ui(i);
            phypp_check(ui < data.size(), "operator[]: index out of bounds (", ui, " vs. ",
                data.size(), ")");
            return ui;
        }

        template<std::size_t D, typename T>
        T to_idx_(T ui, meta::cte_t<false>) const {
            phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
                dims[D], ")");
            return ui;
        }

        template<std::size_t D, typename T>
        uint_t to_idx_(T i, meta::cte_t<true>) const {
            if (i < 0) i += dims[D];
            phypp_check(i >= 0, "operator(): index out of bounds (", i+data.size(), " vs. ",
                dims[D], ")");
            uint_t ui(i);
            phypp_check(ui < dims[D], "operator(): index out of bounds (", ui, " vs. ",
                dims[D], ")");
            return ui;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
            return to_idx_(ix, meta::cte_t<std::is_signed<T>::value>());
        }

        template<std::size_t D, typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        typename std::conditional<std::is_signed<T>::value, uint_t, T>::type to_idx(T ix) const {
            return to_idx_<D>(ix, meta::cte_t<std::is_signed<T>::value>());
        }

        uint_t pitch(uint_t i) const {
            uint_t p = 1;
            for (uint_t j = i+1; j < Dim; ++j) {
                p *= dims[j];
            }
            return p;
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        Type& operator [] (T i) {
            return const_cast<Type&>(const_cast<const vec&>(*this)[i]);
        }

        template<typename T, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        const Type& operator [] (T i) const {
            return *data[to_idx(i)];
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T>& i) {
            vec<1,Type*> v(impl::vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = data[to_idx(i.safe[j])];
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (std::initializer_list<T> l) {
            return operator[](vec<1,T>{l});
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,Type*> operator [] (const vec<1,T*>& i) {
            vec<1,Type*> v(impl::vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = data[to_idx(i.safe[j])];
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T>& i) const {
            vec<1,const Type*> v(impl::vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = data[to_idx(i.safe[j])];
            }
            return v;
        }

        template<typename T = uint_t, typename enable =
            typename std::enable_if<std::is_integral<T>::value>::type>
        vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
            vec<1,const Type*> v(impl::vec_ref_tag, parent);
            v.data.resize(i.data.size());
            v.dims[0] = i.data.size();
            for (uint_t j : range(i)) {
                v.data[j] = data[to_idx(i.safe[j])];
            }
            return v;
        }

        vec<1,Type*> operator [] (impl::range_impl::full_range_t rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::left_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::right_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,Type*> operator [] (const impl::range_impl::left_right_range_t& rng) {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (impl::range_impl::full_range_t rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::left_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::right_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        vec<1,const Type*> operator [] (const impl::range_impl::left_right_range_t& rng) const {
            return impl::vec_access::bracket_access(*this, rng);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) ->
            typename impl::vec_access::helper<true, false, Dim, Type*, Args...>::type {
            static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return impl::vec_access::helper<true, false, Dim, Type*, Args...>::access(*this, i...);
        }

        template<typename ... Args>
        auto operator () (const Args& ... i) const ->
            typename impl::vec_access::helper<true, true, Dim, Type*, Args...>::type {
            static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                "wrong number of indices for this vector");
            return impl::vec_access::helper<true, true, Dim, Type*, Args...>::access(*this, i...);
        }

        struct safe_proxy {
            vec& parent;

            safe_proxy(vec& p) : parent(p) {}
            safe_proxy(const vec&) = delete;
            safe_proxy(vec&&) = delete;
            safe_proxy& operator=(const vec&) = delete;
            safe_proxy& operator=(vec&&) = delete;

            Type& operator [] (uint_t i) {
                return const_cast<Type&>(const_cast<const safe_proxy&>(*this)[i]);
            }

            const Type& operator [] (uint_t i) const {
                return *parent.data[i];
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,Type*> operator [] (const vec<1,T>& i) {
                vec<1,Type*> v(impl::vec_ref_tag, parent.parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = parent.data[i.safe[j]];
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_integral<T>::value>::type>
            vec<1,Type*> operator [] (std::initializer_list<T> l) {
                return operator[](vec<1,T>{l});
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,Type*> operator [] (const vec<1,T*>& i) {
                vec<1,Type*> v(impl::vec_ref_tag, parent.parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = parent.data[i.safe[j]];
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,const Type*> operator [] (const vec<1,T>& i) const {
                vec<1,const Type*> v(impl::vec_ref_tag, parent.parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = parent.data[i.safe[j]];
                }
                return v;
            }

            template<typename T = uint_t, typename enable =
                typename std::enable_if<std::is_unsigned<T>::value>::type>
            vec<1,const Type*> operator [] (const vec<1,T*>& i) const {
                vec<1,const Type*> v(impl::vec_ref_tag, parent.parent);
                v.data.resize(i.data.size());
                v.dims[0] = i.data.size();
                for (uint_t j : range(i)) {
                    v.data[j] = parent.data[i.safe[j]];
                }
                return v;
            }

            vec<1,Type*> operator [] (impl::range_impl::full_range_t rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::left_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::right_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,Type*> operator [] (const impl::range_impl::left_right_range_t& rng) {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (impl::range_impl::full_range_t rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::left_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::right_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            vec<1,const Type*> operator [] (const impl::range_impl::left_right_range_t& rng) const {
                return impl::vec_access::bracket_access(parent, rng);
            }

            template<typename ... Args>
            auto operator () (const Args& ... i) ->
                typename impl::vec_access::helper<false, false, Dim, Type*, Args...>::type {
                static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                    "wrong number of indices for this vector");
                return impl::vec_access::helper<false, false, Dim, Type*, Args...>::access(parent, i...);
            }

            template<typename ... Args>
            auto operator () (const Args& ... i) const ->
                typename impl::vec_access::helper<false, true, Dim, Type*, Args...>::type {
                static_assert(impl::vec_access::accessed_dim<Args...>::value == Dim,
                    "wrong number of indices for this vector");
                return impl::vec_access::helper<false, true, Dim, Type*, Args...>::access(parent, i...);
            }
        } safe;

        effective_type operator + () const {
            return *this;
        }

        effective_type operator - () const {
            effective_type v = *this;
            for (auto& t : v) {
                t = -t;
            }

            return v;
        }

        #define OPERATOR(op) \
            template<typename U> \
            vec& operator op (const vec<Dim,U>& u) { \
                phypp_check(dims == u.dims, "incompatible dimensions in operator '" #op \
                    "' (", dims, " vs ", u.dims, ")"); \
                if (u.view_same(*this)) { \
                    std::vector<dtype> t; t.resize(data.size()); \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        t[i] = u.safe[i]; \
                    } \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        *data[i] op t[i]; \
                    } \
                } else { \
                    for (uint_t i = 0; i < data.size(); ++i) { \
                        *data[i] op u.safe[i]; \
                    } \
                } \
                return *this; \
            } \
            \
            template<typename U> \
            vec& operator op (U u) { \
                for (auto& v : data) { \
                    *v op u; \
                } \
                return *this; \
            }

        OPERATOR(*=)
        OPERATOR(/=)
        OPERATOR(%=)
        OPERATOR(+=)
        OPERATOR(-=)

        #undef OPERATOR

        using iterator = typename impl::vec_iterator_type<vec,impl::default_iterator_policy<vec>>::iterator;
        using const_iterator = typename impl::vec_iterator_type<vec,impl::default_iterator_policy<vec>>::const_iterator;

        iterator begin() {
            return data.begin();
        }

        iterator end() {
            return data.end();
        }

        const_iterator begin() const {
            return data.begin();
        }

        const_iterator end() const {
            return data.end();
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::iterator begin(Property p = Property()) {
            return typename impl::vec_iterator_type<vec,Property>::iterator(
                data.begin(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::iterator end(Property p = Property()) {
            return typename impl::vec_iterator_type<vec,Property>::iterator(
                data.end(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::const_iterator begin(Property p = Property()) const {
            return typename impl::vec_iterator_type<vec,Property>::const_iterator(
                data.begin(), p
            );
        }

        template<typename Property>
        typename impl::vec_iterator_type<vec,Property>::const_iterator end(Property p = Property()) const {
            return typename impl::vec_iterator_type<vec,Property>::const_iterator(
                data.end(), p
            );
        }

        template<typename ... Args>
        impl::vec_access::strided_range<vec,impl::strided_iterator_policy<vec>> stride(const Args& ... i) {
            return impl::vec_access::strided_range<vec,impl::strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<const vec,impl::strided_iterator_policy<vec>> stride(const Args& ... i) const {
            return impl::vec_access::strided_range<const vec,impl::strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<vec,impl::raw_strided_iterator_policy<vec>> raw_stride(const Args& ... i) {
            return impl::vec_access::strided_range<vec,impl::raw_strided_iterator_policy<vec>>(*this, i...);
        }

        template<typename ... Args>
        impl::vec_access::strided_range<const vec,impl::raw_strided_iterator_policy<vec>> raw_stride(const Args& ... i) const {
            return impl::vec_access::strided_range<const vec,impl::raw_strided_iterator_policy<vec>>(*this, i...);
        }
    };


    ////////////////////////////////////////////
    //            Type shortcuts              //
    ////////////////////////////////////////////

    // Shortcuts for the most used types: vec1f, vec2d, ... up to vec6.
    #define MAKE_TYPEDEFS(N) \
        using vec##N##f = vec<N, float>; \
        using vec##N##d = vec<N, double>; \
        using vec##N##i = vec<N, int_t>; \
        using vec##N##u = vec<N, uint_t>; \
        using vec##N##s = vec<N, std::string>; \
        using vec##N##c = vec<N, char>; \
        using vec##N##b = vec<N, bool>;

    MAKE_TYPEDEFS(1)
    MAKE_TYPEDEFS(2)
    MAKE_TYPEDEFS(3)
    MAKE_TYPEDEFS(4)
    MAKE_TYPEDEFS(5)
    MAKE_TYPEDEFS(6)

    #undef MAKE_TYPEDEFS
}

#define PHYPP_INCLUDING_CORE_VEC_BITS
#include "phypp/core/bits/operators.hpp"
#undef PHYPP_INCLUDING_CORE_VEC_BITS

#endif

