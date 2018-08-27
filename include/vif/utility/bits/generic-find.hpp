#ifndef VIF_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "vif/utilty/generic.hpp" instead
#endif

namespace vif {
    // Return the indices of the vector where the value is 'true'.
    template<std::size_t Dim, typename Type, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<Type>,bool>::value>::type>
    vec1u where(const vec<Dim,Type>& v) {
        vec1u ids;
        ids.data.reserve(v.size());
        for (uint_t i : range(v)) {
            if (v.safe[i]) {
                ids.data.push_back(i);
            }
        }

        ids.dims[0] = ids.data.size();
        return ids;
    }

    // Returns the position of the first value in the vector that is 'true'
    // Returns 'npos' if no value is found
    template<std::size_t Dim, typename Type, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<Type>,bool>::value>::type>
    uint_t where_first(const vec<Dim,Type>& v) {
        if (v.empty()) return npos;

        for (uint_t i : range(v)) {
            if (v.safe[i]) {
                return i;
            }
        }

        return npos;
    }

    // Returns the position of the last value in the vector that is 'true'
    // Returns 'npos' if no value is found
    template<std::size_t Dim, typename Type, typename enable =
        typename std::enable_if<std::is_same<meta::rtype_t<Type>,bool>::value>::type>
    uint_t where_last(const vec<Dim,Type>& v) {
        if (v.empty()) return npos;

        for (uint_t i : range(v)) {
            if (v.safe[v.size()-1-i]) {
                return v.size()-1-i;
            }
        }

        return npos;
    }

    // Build the complement of a set of indices, i.e. return the indices that are not part of the
    // provided set.
    template<std::size_t Dim, typename Type>
    vec1u complement(const vec<Dim,Type>& v, const vec1u& ids) {
        if (ids.size() == v.size()) return vec1u();
        vif_check(ids.size() < v.size(), "incompatible size of ids (", ids.size(),
            " vs. ", v.size(), ")");

        vec1b sel(v.size());
        for (uint_t i : ids) {
            sel[i] = true;
        }

        vec1u res; res.reserve(v.size() - ids.size());
        for (uint_t i : range(v)) {
            if (!sel.safe[i]) {
                res.push_back(i);
            }
        }

        return res;
    }

    // In a sorted vector, return the first indices of each non unique sequence, effectively returning
    // indices to all the unique values in the vector.
    // By construction, the returned indices point to sorted values in the original vector.
    template<std::size_t Dim, typename Type>
    vec1u unique_ids_sorted(const vec<Dim,Type>& v) {
        vec1u r;
        if (v.empty()) return r;
        r.reserve(v.size()/4);

        meta::rtype_t<Type> last = v.safe[0];
        r.push_back(0);
        for (uint_t i : range(1, v.size())) {
            if (v.safe[i] != last) {
                r.push_back(i);
                last = v.safe[i];
            }
        }

        r.data.shrink_to_fit();
        return r;
    }

    // In a vector, return indices to all the unique values.
    // The returned indices are are sorted by value.
    template<std::size_t Dim, typename Type>
    vec1u unique_ids(const vec<Dim,Type>& v) {
        vec1u r;
        if (v.empty()) return r;
        r.reserve(v.size()/4);

        vec1u sid = sort(v);
        meta::rtype_t<Type> last = v.safe[sid.safe[0]];
        r.push_back(sid[0]);
        for (uint_t ti : range(1, sid.size())) {
            uint_t i = sid.safe[ti];
            if (v.safe[i] != last) {
                r.push_back(i);
                last = v.safe[i];
            }
        }

        r.data.shrink_to_fit();
        return r;
    }

    // In a vector, return indices to all the unique values. This version takes a second
    // argument with indices that sort the input vector, so that the function does not have
    // to compute the sort.
    // The returned indices are positions in the original (unsorted) vector, and are sorted by
    // value.
    template<std::size_t Dim, typename Type>
    vec1u unique_ids(const vec<Dim,Type>& v, const vec1u& sid) {
        vec1u r;
        if (sid.empty()) return r;
        r.reserve(v.size()/4);

        meta::rtype_t<Type> last = v[sid.safe[0]];
        r.push_back(sid.safe[0]);
        for (uint_t ti : range(1, sid.size())) {
            uint_t i = sid.safe[ti];
            if (v[i] != last) {
                r.push_back(i);
                last = v.safe[i];
            }
        }

        r.data.shrink_to_fit();
        return r;
    }

    // In a sorted vector, return all the unique values.
    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values_sorted(const vec<Dim,Type>& v) {
        return v.safe[unique_ids_sorted(v)];
    }

    // In a vector, return all the unique values.
    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values(const vec<Dim,Type>& v) {
        return v.safe[unique_ids(v)];
    }

    // In a vector, return all the unique values. This version takes a second
    // argument with indices that sort the input vector, so that the function does not have
    // to compute the sort.
    template<std::size_t Dim, typename Type>
    vec<1,meta::rtype_t<Type>> unique_values(const vec<Dim,Type>& v, const vec1u& sid) {
        return v.safe[unique_ids(v, sid)];
    }

    // For each value of the first vector, return 'true' if it is equal to any of the values of the
    // second vector, and 'false' otherwise.
    template<typename Type1, std::size_t Dim2, typename Type2 = Type1>
    bool is_any_of(const Type1& v1, const vec<Dim2,Type2>& v2) {
        for (uint_t j : range(v2)) {
            if (v1 == v2.safe[j]) {
                return true;
            }
        }

        return false;
    }

    // For each value of the first vector, return 'true' if it is equal to any of the values of the
    // second vector, and 'false' otherwise.
    template<std::size_t Dim1, typename Type1, std::size_t Dim2 = Dim1, typename Type2 = Type1>
    vec<Dim1,bool> is_any_of(const vec<Dim1,Type1>& v1, const vec<Dim2,Type2>& v2) {
        vec<Dim1,bool> r(v1.dims);
        for (uint_t i : range(v1))
        for (uint_t j : range(v2)) {
            if (v1.safe[i] == v2.safe[j]) {
                r.safe[i] = true;
                break;
            }
        }

        return r;
    }

    // Compare the two provided vectors and push indices where the two match into 'id1' and 'id2'.
    // If both 'v1' and 'v2' are only composed of unique values, this function is symmetric. Else,
    // only the index of the first value of 'v2' that matches that of 'v1' is stored.
    template<std::size_t D1 = 1, std::size_t D2 = 1, typename Type1, typename Type2 = Type1>
    void match(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2, vec1u& id1, vec1u& id2) {
        uint_t n1 = v1.size();
        uint_t n2 = v2.size();
        uint_t n = std::min(n1, n2);

        id1.clear();
        id2.clear();
        id1.data.reserve(n/2);
        id2.data.reserve(n/2);

        // Choose algorithm based on data size (value chosen by profiling)
        if (n < 64) {
            // Brute force
            for (uint_t i : range(n1))
            for (uint_t j : range(n2)) {
                if (v2.safe[j] == v1.safe[i]) {
                    id1.data.push_back(i);
                    id2.data.push_back(j);
                    break;
                }
            }

            id1.dims[0] = id1.data.size();
            id2.dims[0] = id2.data.size();
        } else {
            // Sort the two vectors and traverse
            vec1u ids1 = sort(v1);
            vec1u ids2 = sort(v2);

            uint_t i1 = 0, i2 = 0;
            while (i1 < n1 && i2 < n2) {
                auto tv1 = v1.safe[ids1.safe[i1]];
                auto tv2 = v2.safe[ids2.safe[i2]];
                if (tv1 < tv2) {
                    ++i1;
                } else if (tv1 > tv2) {
                    ++i2;
                } else {
                    id1.data.push_back(ids1.safe[i1]);
                    id2.data.push_back(ids2.safe[i2]);
                    ++i1; ++i2;
                }
            }

            id1.dims[0] = id1.data.size();
            id1.data.shrink_to_fit();
            id2.dims[0] = id2.data.size();
            id2.data.shrink_to_fit();
        }
    }

    // Returns a new vector containing only the values that are found in common in v1 and v2.
    // Values are returned sorted, duplicates are preserved. Assumes that input vectors
    // are sorted.
    template<std::size_t D1 = 1, std::size_t D2 = 1, typename Type1, typename Type2 = Type1>
    auto set_intersection_sorted(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2) ->
        vec<1,decltype(v1[0]*v2[0])> {

        vec<1,decltype(v1[0]*v2[0])> ret;

        std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
            std::back_inserter(ret.data), typename vec<D1,Type1>::comparator_less()
        );

        ret.dims[0] = ret.data.size();

        return ret;
    }

    namespace impl {
        template<std::size_t D, typename T>
        vec<D,T> make_sorted_(vec<D,T> v) {
            inplace_sort(v);
            return std::move(v);
        }

        template<std::size_t D, typename T>
        vec<D,meta::rtype_t<T>> make_sorted_(vec<D,T*> v) {
            auto tv = v.concretize();
            inplace_sort(tv);
            return std::move(tv);
        }
    }

    // Same as set_intersection_sorted() but without the assumption of sorted input.
    template<std::size_t D1 = 1, std::size_t D2 = 1, typename Type1, typename Type2 = Type1>
    auto set_intersection(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2) ->
        vec<1,decltype(v1[0]*v2[0])> {

        return set_intersection_sorted(impl::make_sorted_(v1), impl::make_sorted_(v2));
    }

    // Returns a new vector containing the values from v1 and v2, with values found in both
    // vectors merged into one. Values are returned sorted, duplicates are preserved.
    // Assumes that input vectors are sorted.
    template<std::size_t D1 = 1, std::size_t D2 = 1, typename Type1, typename Type2 = Type1>
    auto set_union_sorted(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2) ->
        vec<1,decltype(v1[0]*v2[0])> {

        vec<1,decltype(v1[0]*v2[0])> ret;

        std::set_union(v1.begin(), v1.end(), v2.begin(), v2.end(),
            std::back_inserter(ret.data), typename vec<D1,Type1>::comparator_less()
        );

        ret.dims[0] = ret.data.size();

        return ret;
    }

    // Same as set_intersection_sorted() but without the assumption of sorted input.
    template<std::size_t D1 = 1, std::size_t D2 = 1, typename Type1, typename Type2 = Type1>
    auto set_union(const vec<D1,Type1>& v1, const vec<D2,Type2>& v2) ->
        vec<1,decltype(v1[0]*v2[0])> {

        return set_union_sorted(impl::make_sorted_(v1), impl::make_sorted_(v2));
    }

    // Returns the position of the last value in the array that is less than or equal to 'x'.
    // Returns 'npos' if no value satisfy this criterium.
    // Note: assumes that 'v' is sorted and does not contain NaN values.
    template<typename T, std::size_t Dim, typename Type>
    uint_t lower_bound(const vec<Dim,Type>& v, T x) {
        if (v.empty()) return npos;

        auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
            typename vec<Dim,Type>::comparator_less());

        if (iter == v.data.begin()) {
            return npos;
        } else {
            return iter - v.data.begin() - 1;
        }
    }

    // Returns the position of the first value in the array that is greater than 'x'.
    // Returns 'npos' if no value satisfy this criterium.
    // Note: assumes that 'v' is sorted and does not contain NaN values.
    template<typename T, std::size_t Dim, typename Type>
    uint_t upper_bound(const vec<Dim,Type>& v, T x) {
        if (v.empty()) return npos;

        auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
            typename vec<Dim,Type>::comparator_less());

        if (iter == v.data.end()) {
            return npos;
        } else {
            return iter - v.data.begin();
        }
    }

    // Return the position of the last value in 'v' that is less than or equal to 'x' and
    // the position of the first value in 'v' that is greater than 'x'.
    // Note: assumes that 'v' is sorted and does not contain NaN values.
    template<typename T, std::size_t Dim, typename Type>
    std::array<uint_t,2> bounds(const vec<Dim,Type>& v, T x) {
        std::array<uint_t,2> res;

        if (v.empty()) {
            res[0] = npos;
            res[1] = npos;
        } else {
            auto iter = std::upper_bound(v.data.begin(), v.data.end(), x,
                typename vec<Dim,Type>::comparator_less());

            if (iter == v.data.end()) {
                res[0] = v.size() - 1;
                res[1] = npos;
            } else if (iter == v.data.begin()) {
                res[0] = npos;
                res[1] = 0;
            } else {
                res[1] = iter - v.data.begin();
                res[0] = res[1] - 1;
            }
        }

        return res;
    }

    // Return the position of the last value in 'v' that is less than or equal to 'x1' and
    // the position of the first value in 'v' that is greater than 'x2'.
    // Note: assumes that:
    //  1) 'v' is sorted and does not contain NaN values,
    //  2) 'x2' is greater than or equal to 'x1'.
    template<typename T, typename U, std::size_t Dim, typename Type>
    std::array<uint_t,2> bounds(const vec<Dim,Type>& v, T x1, U x2) {
        std::array<uint_t,2> res;

        if (v.empty()) {
            res[0] = npos;
            res[1] = npos;
        } else {
            auto iter = std::upper_bound(v.data.begin(), v.data.end(), x1,
                typename vec<Dim,Type>::comparator_less());

            if (iter == v.data.begin()) {
                res[0] = npos;
            } else {
                res[0] = iter - v.data.begin() - 1;
            }

            iter = std::upper_bound(iter, v.data.end(), x2,
                typename vec<1,Type>::comparator_less());

            if (iter == v.data.end()) {
                res[1] = npos;
            } else {
                res[1] = iter - v.data.begin();
            }
        }

        return res;
    }

    // Return the indices of all the first and last values in the array that are equal to 'x'.
    // Note: assumes that 'v' is sorted and does not contain NaN values.
    template<typename T, std::size_t Dim, typename Type>
    std::array<uint_t,2> equal_range(const vec<Dim,Type>& v, T x) {
        std::array<uint_t,2> res;
        auto tres = std::equal_range(v.data.begin(), v.data.end(), x,
            typename vec<Dim,Type>::comparator_less());

        if (tres.first != v.data.end()) {
            res[0] = tres.first - v.data.begin();
            res[1] = tres.second - v.data.begin();
        } else {
            res[0] = res[1] = npos;
        }

        return res;
    }

    // From a starting position on a 2D boolean map, find the closest point that is 'true'.
    inline bool astar_find(const vec2b& map, uint_t& x, uint_t& y) {
        vif_check(!map.empty(), "this algorithm requires a non empty 2D vector");

        if (x >= map.dims[0]) x = map.dims[0]-1;
        if (y >= map.dims[1]) y = map.dims[1]-1;

        if (map.safe(x,y)) return true;

        using vec_pair = vec<1,std::pair<uint_t,uint_t>>;

        vec_pair open;
        open.push_back(std::make_pair(x,y));

        vec2b visit(map.dims);
        visit.safe(x,y) = true;

        while (!open.empty()) {
            vec_pair old_open = std::move(open);

            for (auto& p : old_open) {
                int_t ox = p.first;
                int_t oy = p.second;

                for (uint_t d : range(4)) {
                    int_t tnx, tny;
                    if (d == 0) {
                        tnx = ox;   tny = oy+1;
                    } else if (d == 1) {
                        tnx = ox+1; tny = oy;
                    } else if (d == 2) {
                        tnx = ox;   tny = oy-1;
                    } else {
                        tnx = ox-1; tny = oy;
                    }

                    if (tnx < 0 || tny < 0) continue;

                    x = tnx;
                    y = tny;
                    if (x >= map.dims[0] || y >= map.dims[1] || visit.safe(x,y)) continue;

                    if (!map.safe(x,y)) {
                        open.push_back(std::make_pair(x,y));
                        visit.safe(x,y) = true;
                    } else {
                        return true;
                    }
                }
            }
        }

        return false;
    }
}
