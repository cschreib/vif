#ifndef REFLEX_HELPERS_HPP
#define REFLEX_HELPERS_HPP

#include "reflex.hpp"
#include "print.hpp"
#include "string.hpp"

namespace reflex {
    namespace impl {
        template<typename T>
        struct do_print {
            explicit do_print(T& t) : o(t) {}

            T& o;
            bool first = true;

            template<typename P>
            void operator() (const member_t& m, const P& t) {
                if (!this->first) this->o << ", ";
                this->o << m.name << "=" << t;
                this->first = false;
            }
        };
    }

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

        o << "{ ";
        reflex::foreach_member(d, impl::do_print<O>(o));
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
void merge_elements_(reflex::struct_t<T> t, reflex::struct_t<U> u);

namespace reflex {
    namespace impl {
        template<typename M, bool Struct>
        struct do_run_;

        template<typename M>
        struct do_run_<M, false> {
            explicit do_run_(const reflex::member_t& rm, const M& tv) : m(rm), v(tv) {}

            const reflex::member_t& m;
            const M& v;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    merge_elements(p, this->v);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename M>
        struct do_run_<M, true> {
            explicit do_run_(const reflex::member_t& rm, const M& tv) : m(rm), v(tv) {}

            const reflex::member_t& m;
            const M& v;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    merge_elements(p, this->v);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, reflex::struct_t<P> p) {
                merge_elements_(p, this->v);
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename M>
        using do_run = do_run_<M, reflex::is_struct<M>::value>;

        template<typename T>
        struct do_merge {
            explicit do_merge(reflex::struct_t<T> s) : t(s) {}

            reflex::struct_t<T> t;

            template<typename M>
            void operator () (const reflex::member_t& m, const M& v) {
                reflex::foreach_member(this->t, do_run<M>(m, v));
            }
        };
    }
}

template<typename T, typename U>
void merge_elements_(reflex::struct_t<T> t, reflex::struct_t<U> u) {
    reflex::foreach_member(u, reflex::impl::do_merge<T>(t));
}

template<typename T, typename U>
void merge_elements(T& t, const U& u) {
    merge_elements_(reflex::wrap(t), reflex::wrap(u));
}

template<typename T, typename U>
void merge_elements(T& t, const U& u, const vec1u& ids);

template<typename T, typename U>
void merge_elements_(T& t, const U& u, const vec1u& ids) {
    t = u;
}

template<typename T, typename U>
void merge_elements_(vec_t<1,T>& t, const vec_t<1,U>& u, const vec1u& ids) {
    t = u[ids];
}

template<typename T, typename U>
void merge_elements_(vec_t<2,T>& t, const vec_t<2,U>& u, const vec1u& ids) {
    t = u(ids,_);
}

template<typename T, typename U>
void merge_elements_(reflex::struct_t<T> t, reflex::struct_t<U> u, const vec1u& ids);

namespace reflex {
    namespace impl {
        template<typename M, bool Struct>
        struct do_run_ids_;

        template<typename M>
        struct do_run_ids_<M, false> {
            do_run_ids_(const reflex::member_t& tm, const M& tv, const vec1u& i) :
                m(tm), v(tv), ids(i) {}

            const reflex::member_t& m;
            const M& v;
            const vec1u& ids;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    merge_elements(p, this->v, this->ids);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename M>
        struct do_run_ids_<M, true> {
            do_run_ids_(const reflex::member_t& tm, const M& tv, const vec1u& i) :
                m(tm), v(tv), ids(i) {}

            const reflex::member_t& m;
            const M& v;
            const vec1u& ids;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    merge_elements(p, this->v, this->ids);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, reflex::struct_t<P> p) {
                merge_elements_(p, this->v, this->ids);
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename M>
        using do_run_ids = do_run_ids_<M, reflex::is_struct<M>::value>;

        template<typename T>
        struct do_merge_ids {
            do_merge_ids(reflex::struct_t<T> tt, const vec1u& i) : t(tt), ids(i) {}

            reflex::struct_t<T> t;
            const vec1u& ids;

            template<typename M>
            void operator () (const reflex::member_t& m, const M& v) {
                reflex::foreach_member(this->t, do_run_ids<M>(m, v, this->ids));
            }
        };
    }
}

template<typename T, typename U>
void merge_elements_(reflex::struct_t<T> t, reflex::struct_t<U> u, const vec1u& ids) {
    reflex::foreach_member(u, reflex::impl::do_merge_ids<T>(t, ids));
}

template<typename T, typename U>
void merge_elements(T& t, const U& u, const vec1u& ids) {
    merge_elements_(reflex::wrap(t), reflex::wrap(u), ids);
}

#endif
