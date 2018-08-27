#ifndef VIF_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "vif/utilty/generic.hpp" instead
#endif

namespace vif {
    // Get multi-dim IDs from a flat ID
    template<std::size_t D>
    vec1u mult_ids(const std::array<uint_t,D>& dims, uint_t i) {
        vec1u r(D);

        for (uint_t j : range(D)) {
            r.safe[D-1-j] = i % dims[D-1-j];
            i /= dims[D-1-j];
        }

        return r;
    }

    template<std::size_t D>
    vec2u mult_ids(const std::array<uint_t,D>& dims, vec1u i) {
        vec2u r(D,i.size());

        for (uint_t j : range(D)) {
            r.safe(D-1-j,_) = i % dims[D-1-j];
            i /= dims[D-1-j];
        }

        return r;
    }

    template<std::size_t D, typename T>
    vec1u mult_ids(const vec<D,T>& v, uint_t i) {
        return mult_ids(v.dims, i);
    }

    template<std::size_t D, typename T>
    vec2u mult_ids(const vec<D,T>& v, vec1u i) {
        return mult_ids(v.dims, i);
    }

    // Get flat ID from multi-dim IDs
    namespace impl {
        template<std::size_t D>
        uint_t flat_id_(const std::array<uint_t,D>&, uint_t ret, meta::cte_t<D>) {
            return ret;
        }

        template<std::size_t D, std::size_t I, typename U, typename ... Args>
        uint_t flat_id_(const std::array<uint_t,D>& dims, uint_t ret,
            meta::cte_t<I>, U i, Args&& ... args) {

            uint_t di = impl::vec_access::pitch(dims, I)*
                impl::vec_access::template to_idx<I>(dims, i);
            return flat_id_(dims, ret + di, meta::cte_t<I+1>{},
                std::forward<Args>(args)...);
        }

        template<std::size_t D, std::size_t I, typename TI>
        uint_t flat_id_(const std::array<uint_t,D>& dims, uint_t ret,
            meta::cte_t<I>, const vec<1,TI>& ids) {

            uint_t di = impl::vec_access::pitch(dims, I)*
                impl::vec_access::template to_idx<I>(dims, ids[I]);
            return flat_id_(dims, ret + di, meta::cte_t<I+1>{}, ids);
        }
    }

    template<std::size_t D, typename ... Args>
    uint_t flat_id(const std::array<uint_t,D>& dims, Args&& ... args) {
        static_assert(sizeof...(Args) == D, "wrong number of IDs provided");
        return impl::flat_id_(dims, 0, meta::cte_t<0>{}, std::forward<Args>(args)...);
    }

    template<std::size_t D, typename T, typename ... Args>
    uint_t flat_id(const vec<D,T>& v, Args&& ... args) {
        return flat_id(v.dims, std::forward<Args>(args)...);
    }

    template<std::size_t D, typename TI>
    uint_t flat_id(const std::array<uint_t,D>& dims, const vec<1,TI>& ids) {
        vif_check(ids.size() == D, "wrong number of IDs provided");
        return impl::flat_id_(dims, meta::cte_t<0>{}, ids);
    }

    template<std::size_t D, typename T, typename TI>
    uint_t flat_id(const vec<D,T>& v, const vec<1,TI>& ids) {
        return flat_id(v.dims, ids);
    }

    // Increment a list of indices, step by step
    // This is an alternative to a recursive loop.
    inline void increment_index_list(vec1u& ids, const vec1u& n) {
        uint_t i = ids.size();
        do {
            --i;
            ++(ids.safe[i]);
            if (ids.safe[i] == n.safe[i]) {
                ids.safe[i] = 0;
            } else {
                break;
            }
        } while (i != 0);
    }

    inline void increment_index_list(vec1u& ids, const uint_t& n) {
        uint_t i = ids.size();
        do {
            --i;
            ++(ids.safe[i]);
            if (ids.safe[i] == n) {
                ids.safe[i] = 0;
            } else {
                break;
            }
        } while (i != 0);
    }
}
