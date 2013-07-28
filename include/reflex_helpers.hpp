#ifndef REFLEX_HELPERS_HPP
#define REFLEX_HELPERS_HPP

#include "reflex.hpp"
#include "print.hpp"
#include "string.hpp"

namespace reflex {
    template<typename T, typename F, std::size_t N>
    void foreach_member_(reflex::struct_t<T> d, F& f, cte_t<N>, cte_t<N>) {}

    template<typename T, typename F, std::size_t I, std::size_t N>
    void foreach_member_(reflex::struct_t<T> d, F& f, cte_t<I>, cte_t<N>) {
        using type = typename reflex::struct_t<T>::member_types::template get<I>;
        auto& m = d.data.members[I];
        f(m, reflex::wrap(*(constify<T,type>*)d.data.members[I].value));
        foreach_member_(d, f, cte_t<I+1>(), cte_t<N>());
    }

    template<typename T, typename F>
    void foreach_member(reflex::struct_t<T> d, F&& tf) {
        F& f = tf;
        foreach_member_(d, f, cte_t<0>(), cte_t<reflex::struct_t<T>::member_count>());
    }

    template<typename O, typename T>
    O& operator << (O& o, reflex::struct_t<T> d) {
        if (reflex::struct_t<T>::member_count == 0) {
            o << "{}";
            return o;
        }

        struct {
            O& o;
            bool first = true;

            template<typename P>
            void operator() (const member_t& m, const P& t) {
                if (!this->first) this->o << ", ";
                this->o << m.name << "=" << t;
                this->first = false;
            }
        } do_print{o};

        o << "{ ";
        reflex::foreach_member(d, do_print);
        o << " }";
        return o;
    }

    template<typename T>
    std::string type_name_of(T&& t) {
        return typeid(typename std::decay<T>::type).name();
    }
}
    
template<typename T, typename U>
void merge_elements(T& t, const U& u);

template<typename T, typename U>
void merge_elements_(T& t, const U& u) {
    t = u;
}

template<typename T, typename U>
void merge_elements_(reflex::struct_t<T> t, reflex::struct_t<U> u) {
    struct {
        reflex::struct_t<T> t;

        template<typename M>
        void operator () (const reflex::member_t& m, const M& v) {
            struct {
                const reflex::member_t& m;
                const M& v;

                void operator () (reflex::member_t& n, M& p) {
                    if (this->m.name == n.name) {
                        merge_elements(p, this->v);
                    }
                }

                template<typename P, typename enable =
                    typename std::enable_if<reflex::is_struct<M>::value>::type>
                void operator () (reflex::member_t& n, reflex::struct_t<P> p) {
                    merge_elements_(p, this->v);
                }

                template<typename P>
                void operator () (reflex::member_t& n, P&& p) {
                    phyu_check(this->m.name != n.name, "incompatible types in merging '",
                        this->m.full_name(), "' into '", n.full_name(), "'"
                    );
                }
            } do_run{m, v};

            reflex::foreach_member(this->t, do_run);
        }
    } do_merge{t};

    reflex::foreach_member(u, do_merge);
}

template<typename T, typename U>
void merge_elements(T& t, const U& u) {
    merge_elements_(reflex::wrap(t), reflex::wrap(u));
}

#endif
