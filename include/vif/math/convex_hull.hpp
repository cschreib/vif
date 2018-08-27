#ifndef VIF_MATH_CONVEX_HULL_HPP
#define VIF_MATH_CONVEX_HULL_HPP

#include "vif/core/vec.hpp"
#include "vif/core/error.hpp"
#include "vif/core/range.hpp"
#include "vif/utility/generic.hpp"
#include "vif/math/base.hpp"


namespace vif {
    // Structure holding the data points of a convex hull
    template <typename T = double>
    struct convex_hull {
        struct no_check {};

        static_assert(!std::is_pointer<T>::value,
            "convex_hull can only be created from plain types, pointers and vector views "
            "are not supported");

        // Default constructed hull is empty
        convex_hull() {}

        // Build a new convex hull and make sure it is valid
        convex_hull(vec<1,T> tx, vec<1,T> ty) : x(std::move(tx)), y(std::move(ty)) {
            validate();
        }

        // Build a hull without security checks
        // example: convex_hull(x, y, no_check{});
        template <typename U = double>
        convex_hull(vec<1,T> tx, vec<1,T> ty, typename convex_hull<U>::no_check) :
            x(std::move(tx)), y(std::move(ty)), validated(true) {}

        // Number of vertices in the hull
        uint_t size() const {
            return x.size();
        }

        // Check the validity of this convex hull
        // Will only be done once unless 'force' is set to 'true'.
        void validate(bool force = false) const {
            if (!validated || force) {
                vif_check(x.dims == y.dims, "incompatible dimensions between X and Y "
                    "(", x.dims, " vs. ", y.dims, ")");
                vif_check(x.size() >= 3, "a hull must have at least 3 elements "
                    "(got ", x.size(), ")");

                const auto eps = std::numeric_limits<T>::epsilon();
                const uint_t hend = size()-1;
                closed = (abs(x.safe[0] - x.safe[hend]) < eps &&
                     abs(y.safe[0] - y.safe[hend]) < eps);

                orient = ((x.safe[1] - x.safe[0])*(y.safe[2] - y.safe[1]) -
                         (y.safe[1] - y.safe[0])*(x.safe[2] - x.safe[1]) > 0) ? 1 : -1;

                validated = true;
            }
        }

        // Create the (ux, uy, nx, ny) vectors if they do not exist yet
        void build_vectors() const {
            if (!vectors_built) {
                ux.resize(size()-1); uy.resize(size()-1);
                nx.resize(size()-1); ny.resize(size()-1);
                length.resize(size()-1);

                for (uint_t i : range(size()-1)) {
                    // Get unit vector
                    ux.safe[i] = x.safe[i+1] - x.safe[i];
                    uy.safe[i] = y.safe[i+1] - y.safe[i];
                    // Get perpendicular vector, pointing inside the hull
                    nx.safe[i] = orient == 1 ? y.safe[i+1] - y.safe[i] : y.safe[i] - y.safe[i+1];
                    ny.safe[i] = orient == 1 ? x.safe[i] - x.safe[i+1] : x.safe[i+1] - x.safe[i];
                    // Normalize
                    auto l = sqrt(sqr(nx.safe[i]) + sqr(ny.safe[i]));
                    length.safe[i] = l;
                    ux.safe[i] /= l; uy.safe[i] /= l; nx.safe[i] /= l; ny.safe[i] /= l;
                }

                vectors_built = true;
            }
        }

        // If manual adjustments are made on the vertices, call this function
        void mark_as_dirty() {
            vectors_built = false;
            validated = false;
        }

        vec<1,T> x, y; // vertices

        // Cached computations
        mutable bool validated = false; // flag if the hull has been validated or not
        mutable bool closed = false;    // is first vertex = last vertex ?
        mutable int_t orient = 0;       // orientation of the hull, 1: ccw, -1: cw

        mutable bool vectors_built = false; // flag if segment data has been computed
        mutable vec<1,T> ux, uy;            // unit vector
        mutable vec<1,T> nx, ny;            // normal vector
        mutable vec<1,T> length;            // vector length
    };

    // Build the convex hull of a set of points, returning the indices of the points that form the hull
    // in counter-clockwise order.
    // Uses the monotone chain algorithm, taken from:
    // http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C.2B.2B
    template<typename T>
    convex_hull<meta::rtype_t<T>> build_convex_hull(const vec<1,T>& x, const vec<1,T>& y) {
        vif_check(x.dims == y.dims, "incompatible dimensions between X and Y "
            "(", x.dims, " vs. ", y.dims, ")");

        uint_t npt = x.size();
        vec1u res = uintarr(2*npt);
        vec1u ids = uindgen(npt);
        std::sort(ids.data.begin(), ids.data.end(), [&](uint_t i, uint_t j) {
            if (x.safe[i] < x.safe[j]) return true;
            if (x.safe[i] > x.safe[j]) return false;
            if (y.safe[i] < y.safe[j]) return true;
            return false;
        });

        auto cross = [&](uint_t i, uint_t j, uint_t k) {
            return (x.safe[j] - x.safe[i])*(y.safe[k] - y.safe[i])
                 - (y.safe[j] - y.safe[i])*(x.safe[k] - x.safe[i]);
        };

        uint_t k = 0;
        for (uint_t i = 0; i < npt; ++i) {
            while (k >= 2 && cross(res.safe[k-2], res.safe[k-1], ids.safe[i]) <= 0) --k;
            res.safe[k] = ids.safe[i];
            ++k;
        }

        uint_t t = k+1;
        for (uint_t i = npt-2; i != npos; --i) {
            while (k >= t && cross(res.safe[k-2], res.safe[k-1], ids.safe[i]) <= 0) --k;
            res.safe[k] = ids.safe[i];
            ++k;
        }

        res.data.resize(k);
        res.dims[0] = k;

        convex_hull<meta::rtype_t<T>> hull(x.safe[res], y.safe[res], convex_hull<>::no_check{});
        hull.orient = 1; // counter-clockwise by construction
        hull.closed = true; // closed by construction

        return hull;
    }

    template<typename TX, typename TY, typename H>
    bool in_convex_hull(const TX& x, const TY& y, const convex_hull<H>& hull) {
        hull.validate();
        vif_check(hull.closed, "the provided hull must be closed");

        // Find out if the hull is built counter-clockwise or not
        for (uint_t i : range(hull.size()-1)) {
            auto cross = (hull.x.safe[i+1] - hull.x.safe[i])*(y - hull.y.safe[i]) -
                         (hull.y.safe[i+1] - hull.y.safe[i])*(x - hull.x.safe[i]);
            if (cross*hull.orient < 0) return false;
        }

        return true;
    }

    template<std::size_t Dim, typename TX, typename TY, typename H>
    vec<Dim,bool> in_convex_hull(const vec<Dim,TX>& x, const vec<Dim,TY>& y,
        const convex_hull<H>& hull) {

        vec<Dim,bool> res = replicate(true, x.dims);

        hull.validate();
        vif_check(hull.closed, "the provided hull must be closed");
        vif_check(x.dims == y.dims, "incompatible dimensions between X and Y "
            "(", x.dims, " vs. ", y.dims, ")");

        for (uint_t i : range(hull.size()-1)) {
            for (uint_t p : range(x)) {
                if (!res.safe[p]) continue;

                auto cross = (hull.x.safe[i+1] - hull.x.safe[i])*(y.safe[p] - hull.y.safe[i]) -
                             (hull.y.safe[i+1] - hull.y.safe[i])*(x.safe[p] - hull.x.safe[i]);

                if (cross*hull.orient < 0) res.safe[p] = false;
            }
        }

        return res;
    }

    // Compute the signed distance of a set of points with respect to the provided convex
    // hull. Positive distances mean that the point lies inside the hull.
    template<typename TX, typename TY, typename H, typename enable =
        typename std::enable_if<!meta::is_vec<TX>::value && !meta::is_vec<TY>::value>::type>
    auto convex_hull_distance(const TX& x, const TY& y, const convex_hull<H>& hull)
        -> decltype(sqrt(x*y)) {

        hull.validate();
        vif_check(hull.closed, "the provided hull must be closed");
        hull.build_vectors();

        decltype(sqrt(x*y)) res = finf;
        bool inhull = true;

        for (uint_t i : range(hull.size()-1)) {
            // Compute signed distance to current hull face line
            auto dx = x - hull.x.safe[i];
            auto dy = y - hull.y.safe[i];
            auto d  = dx*hull.nx.safe[i] + dy*hull.ny.safe[i];

            // Check if the point is inside the hull or not
            if (d > 0) {
                inhull = false;
            }

            d = abs(d);
            if (res > d) {
                // Find if the projection of the point on the face lies on the segment
                auto proj = dx*hull.ux.safe[i] + dy*hull.uy.safe[i];
                if (proj < 0) {
                    // Projection lies before starting point
                    d = sqrt(sqr(dx) + sqr(dy));
                } else if (proj > hull.length.safe[i]) {
                    // Projection lies after ending point
                    d = sqrt(sqr(x - hull.x.safe[i+1]) + sqr(y - hull.y.safe[i+1]));
                }

                // Keep this distance if it is minimal
                if (res > d) res = d;
            }
        }

        // Give negative distance to points outside the hull
        if (!inhull) {
            res = -res;
        }

        return res;
    }

    // Compute the signed distance of a set of points with respect to the provided convex
    // hull. Positive distances mean that the point lies inside the hull.
    template<std::size_t Dim, typename TX, typename TY, typename H>
    auto convex_hull_distance(const vec<Dim,TX>& x, const vec<Dim,TY>& y,
        const convex_hull<H>& hull) -> vec<Dim, decltype(sqrt(x[0]*y[0]))> {

        vif_check(x.dims == y.dims, "incompatible dimensions between X and Y "
            "(", x.dims, " vs. ", y.dims, ")");
        hull.validate();
        vif_check(hull.closed, "the provided must be closed");
        hull.build_vectors();

        vec<Dim,decltype(sqrt(x[0]*y[0]))> res = replicate(finf, x.dims);
        vec<Dim,bool> inhull = replicate(true, x.dims);

        for (uint_t i : range(hull.size()-1)) {
            for (uint_t p : range(x)) {
                // Compute signed distance to current hull face line
                auto dx = x.safe[p] - hull.x.safe[i];
                auto dy = y.safe[p] - hull.y.safe[i];
                auto d  = dx*hull.nx.safe[i] + dy*hull.ny.safe[i];

                // Check if the point is inside the hull or not
                if (d > 0) {
                    inhull.safe[p] = false;
                }

                d = abs(d);
                if (res.safe[p] > d) {
                    // Find if the projection of the point on the face lies on the segment
                    auto proj = dx*hull.ux.safe[i] + dy*hull.uy.safe[i];
                    if (proj < 0) {
                        // Projection lies before starting point
                        d = sqrt(sqr(dx) + sqr(dy));
                    } else if (proj > hull.length.safe[i]) {
                        // Projection lies after ending point
                        d = sqrt(sqr(x.safe[p] - hull.x.safe[i+1]) + sqr(y.safe[p] - hull.y.safe[i+1]));
                    }

                    // Keep this distance if it is minimal
                    if (res.safe[p] > d) res.safe[p] = d;
                }
            }
        }

        // Give negative distance to points outside the hull
        for (uint_t p : range(inhull)) {
            if (!inhull.safe[p]) {
                res.safe[p] *= -1.0;
            }
        }

        return res;
    }
}

#endif
