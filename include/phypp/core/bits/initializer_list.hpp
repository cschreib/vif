#ifndef PHYPP_INCLUDING_CORE_VEC_BITS
#error this file is not meant to be included separately, include "phypp/core/vec.hpp" instead
#endif

namespace phypp {
namespace impl {
namespace vec_ilist {
    template<std::size_t Dim, typename Type>
    struct helper {
        using dtype = meta::dtype_t<Type>;

        static void set_dim_(vec<Dim,Type>& v, meta::nested_initializer_list<1,dtype> il, meta::cte_t<Dim-1>) {
            if (il.size() == 0) {
                for (uint_t i = 0; i < Dim; ++i) {
                    v.dims[i] = 0;
                }
            } else {
                v.dims[Dim-1] = il.size();
            }
        }

        template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
        static void set_dim_(vec<Dim,Type>& v, meta::nested_initializer_list<Dim-N,dtype> il, meta::cte_t<N>) {
            if (il.size() == 0) {
                for (uint_t i = 0; i < Dim; ++i) {
                    v.dims[i] = 0;
                }
            } else {
                v.dims[N] = il.size();
                set_dim_(v, *il.begin(), meta::cte_t<N+1>());
            }
        }

        static void fill_(vec<Dim,Type>& v, meta::nested_initializer_list<1,dtype> il, uint_t& idx, meta::cte_t<Dim-1>) {
            phypp_check(il.size() == v.dims[Dim-1], "heterogeneous intializer lists are not allowed");
            for (auto& t : il) {
                v.data[idx] = t;
                ++idx;
            }
        }

        template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
        static void fill_(vec<Dim,Type>& v, meta::nested_initializer_list<Dim-N,dtype> il, uint_t& idx, meta::cte_t<N>) {
            phypp_check(il.size() == v.dims[N], "heterogeneous intializer lists are not allowed");
            for (auto& t : il) {
                fill_(v, t, idx, meta::cte_t<N+1>());
            }
        }

        static void fill(vec<Dim,Type>& v, meta::nested_initializer_list<Dim,dtype> il) {
            set_dim_(v, il, meta::cte_t<0>());
            v.resize();
            if (!v.empty()) {
                uint_t idx = 0;
                fill_(v, il, idx, meta::cte_t<0>());
            }
        }
    };

    template<std::size_t Dim, typename Type>
    struct helper<Dim, Type*> {
        using dtype = meta::dtype_t<Type>;

        static void fill_(vec<Dim,Type*>& v, meta::nested_initializer_list<1,dtype> il, uint_t& idx, meta::cte_t<Dim-1>) {
            phypp_check(il.size() == v.dims[Dim-1], "heterogeneous intializer lists are not allowed");
            for (auto& t : il) {
                *v.data[idx] = t;
                ++idx;
            }
        }

        template<std::size_t N, typename enable = typename std::enable_if<N != Dim-1>::type>
        static void fill_(vec<Dim,Type*>& v, meta::nested_initializer_list<Dim-N,dtype> il, uint_t& idx, meta::cte_t<N>) {
            phypp_check(il.size() == v.dims[N], "heterogeneous intializer lists are not allowed");
            for (auto& t : il) {
                fill_(v, t, idx, meta::cte_t<N+1>());
            }
        }

        static void fill(vec<Dim,Type*>& v, meta::nested_initializer_list<Dim,dtype> il) {
            uint_t idx = 0;
            fill_(v, il, idx, meta::cte_t<0>());
        }
    };
}
}
}
