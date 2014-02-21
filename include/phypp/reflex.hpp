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
        member_t(std::string name_, const std::type_info& type_, data_t* parent_ = nullptr,
            void* value_ = nullptr, data_t* reflex_ = nullptr) :
            name(name_), type(type_), parent(parent_), value(value_), reflex(reflex_) {}

        std::string name;
        const std::type_info& type;
        data_t* parent = nullptr;
        void* value;
        data_t* reflex = nullptr;

        std::string full_name() const;
    };

    #define REFLEX_MEM_HEADER {'_', 'r', 'e', 'f', 'l', 'e', 'x', '_'}

    struct data_impl_t {
        data_impl_t() = default;
        data_impl_t(const std::type_info& type_, std::string name_,
            std::vector<member_t> members_) : type(type_), name(name_), members(members_) {}
        data_impl_t(data_impl_t&&) = default;
        data_impl_t(const data_impl_t& d) = delete;
        data_impl_t& operator = (data_impl_t&&) = default;
        data_impl_t& operator = (const data_impl_t&) = delete;

        const std::type_info& type;
        std::string name;
        std::vector<member_t> members;
    };

    struct data_t {
        template<typename T>
        void copy_members_(T&& d) {
            members.reserve(d.members.size());
            int_t diff = (char*)this - (const char*)&d;
            for (auto& m : d.members) {
                members.push_back({m.name, m.type, this, (void*)((char*)m.value + diff)});
                if (m.reflex) {
                    members.back().reflex = (data_t*)(void*)((char*)m.reflex + diff);
                    members.back().reflex->parent = this;
                }
            }
        }

        data_t() = default;
        data_t(data_impl_t&& d) : name(d.name), members(std::move(d.members)) {}

        data_t(data_t&& d) : name(std::move(d.name)) {
            copy_members_(std::move(d));
        }

        data_t(const data_t& d) : name(d.name) {
            copy_members_(d);
        }

        data_t& operator = (data_t&& d) { return *this; }
        data_t& operator = (const data_t& d) { return *this; }

        char header[8] = REFLEX_MEM_HEADER;
        const std::type_info& type = typeid(*this);
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
    std::string seek_name(const T& t, uint_t max_dist = 100000) {
        const char* c = reinterpret_cast<const char*>(&t);
        const char* oc = c;
        const char pattern[] = REFLEX_MEM_HEADER;

        // Probe the memory that follows the address of the object and look for a reflection header
        while (uint_t(c - oc) < max_dist) {
            while (*c != pattern[0]) ++c;
            const char* hstart = c;

            bool found = true;
            for (uint_t i = 1; i < sizeof(pattern); ++i) {
                if (*(++c) != pattern[i]) {
                    found = false;
                    break;
                }
            }

            if (!found) continue;

            // If one is found, go back to the corresponding data_t
            const reflex::data_t* data = reinterpret_cast<const reflex::data_t*>(hstart);
            const uint_t offset = &data->header[0] - reinterpret_cast<const char*>(data);
            data = reinterpret_cast<const reflex::data_t*>(hstart - offset);

            // Look for the members of this data_t, and see if one matches both address and type
            // (type is necessary because some objects, although different, share the same address,
            // for example a structure and its first member)
            for (auto& m : data->members) {
                if (m.type == typeid(t) && m.value == (const void*)&t) {
                    return data->full_name() + "." + m.name;
                }
            }

            // If no member matches, then this was not the right data_t, continue searching.
        }

        return "<?>";
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
        m.reflex = &d._reflex;
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
    void do_init_(type_list<Args...> tl, cte_t<0>, T* t, data_impl_t& data) {}

    template<typename ... Args, std::size_t N, typename T>
    void do_init_(type_list<Args...> tl, cte_t<N>, T* t, data_impl_t& data) {
        member_t& m = data.members[N-1];
        m.parent = &t->_reflex;
        auto& v = get_value<N-1>(tl, m);
        do_init_set_name_(m, v, &t->_reflex,
            cte_t<has_reflex<typename std::decay<decltype(v)>::type>::value>()
        );
        do_init_(tl, cte_t<N-1>(), t, data);
    }

    template<typename ... Args, typename T>
    data_impl_t do_init(type_list<Args...> tl, T* t, data_impl_t&& data) {
        data_impl_t d = std::move(data);
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

    #ifdef REFLECTION_STAGE
    template<typename T>
    struct struct_t {
        constify<T,data_t>& data;
        using member_types = empty_t::_reflex_types;
        static const std::size_t member_count = member_types::count;
    };
    #else
    template<typename T>
    struct struct_t {
        constify<T,data_t>& data;
        using member_types = typename T::_reflex_types;
        static const std::size_t member_count = member_types::count;
    };
    #endif

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
            return struct_t<T>{empty._reflex};
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
    auto wrap(T&& t) ->
        decltype(wrap_t<enabled<typename std::decay<T>::type>::value>::wrap(std::forward<T>(t))) {
        return wrap_t<enabled<typename std::decay<T>::type>::value>::wrap(std::forward<T>(t));
    }
}

#define MAKE_MEMBER(name) reflex::member_t{#name, typeid(name), nullptr, (void*)&name}

#define MEMBERS1(...) using _reflex_types = decltype(reflex::make_types_decay_(__VA_ARGS__))

#define MEMBERS2(name, ...) \
    reflex::data_t _reflex = reflex::do_init(_reflex_types(), this, { \
        typeid(*this), name, {__VA_ARGS__}})

#define NO_MEMBER(name) \
    using _reflex_types = decltype(reflex::make_types_()); \
    reflex::data_t _reflex = reflex::data_impl_t{typeid(*this), name, {}}

#endif
