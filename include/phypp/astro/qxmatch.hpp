#ifndef PHYPP_ASTRO_QXMATCH_HPP
#define PHYPP_ASTRO_QXMATCH_HPP

#include "phypp/astro/astro.hpp"

namespace phypp {
namespace astro {
    struct qxmatch_res {
        vec2u id;
        vec2d d;
        vec1u rid;
        vec1d rd;

        // Reflection data
        MEMBERS1(id, d, rid, rd);
        MEMBERS2("qxmatch_res", MAKE_MEMBER(id), MAKE_MEMBER(d), MAKE_MEMBER(rid), MAKE_MEMBER(rd));
    };

    void qxmatch_save(const std::string& file, const qxmatch_res& r) {
        fits::write_table(file, ftable(r.id, r.d, r.rid, r.rd));
    }

    qxmatch_res qxmatch_restore(const std::string& file) {
        qxmatch_res r;
        fits::read_table(file, ftable(r.id, r.d, r.rid, r.rd));
        return r;
    }

    struct qxmatch_params {
        uint_t thread = 1u;
        uint_t nth = 1u;
        bool verbose = false;
        bool self = false;
        bool brute_force = false;
        bool no_mirror = false;
        bool linear = false;
    };
}

namespace impl {
    namespace qxmatch_impl {
        struct depth_cache {
            struct depth_t {
                vec1i bx, by;
                double max_dist;
            };

            std::vector<depth_t> depths;
            vec2b                visited;
            double               csize;
            std::size_t          max_depth;

            depth_cache(uint_t nx, uint_t ny, double cs) : csize(cs),
                max_depth(std::max(nx, ny)) {

                visited = vec2b(nx, ny);
                depths.reserve(20);

                // First depth is trivial
                depths.push_back(depth_t{});
                auto& depth = depths.back();
                depth.max_dist = csize/2.0;
                depth.bx = {0};
                depth.by = {0};
                visited(0,0) = true;

                // Generate a few in advance
                for (uint_t d : range(std::min(uint_t{9}, std::max(nx, ny)-1))) {
                    grow();
                }
            }

            const depth_t& operator[] (uint_t i) {
                if (i >= depths.size()) {
                    for (uint_t d : range(i+1 - depths.size())) {
                        grow();
                    }
                }

                return depths[i];
            }

            std::size_t size() const {
                return max_depth;
            }

            void grow() {
                phypp_check(depths.size() < visited.dims[0] ||
                    depths.size() < visited.dims[1], "cannot grow past the size of the "
                    "bucket field (", visited.dims, "; trying to reach ", depths.size(), ")");

                depths.push_back(depth_t{});
                auto& depth = depths.back();
                uint_t d = depths.size();

                // Look further by one cell size
                depth.max_dist = csize*(d + 0.5);

                // See which new buckets are reached.
                // We only need to do the maths for one quadrant, the other 3
                // can be deduced by symmetry.
                for (uint_t x : range(std::min(d+1, visited.dims[0])))
                for (uint_t y : range(std::min(d+1, visited.dims[1]))) {
                    double dist2 = sqr(csize)*min(vec1d{
                        sqr(x-0.5) + sqr(y-0.5), sqr(x+0.5) + sqr(y-0.5),
                        sqr(x+0.5) + sqr(y+0.5), sqr(x-0.5) + sqr(y+0.5)
                    });

                    if (dist2 <= sqr(depth.max_dist) && !visited(x, y)) {
                        depth.bx.push_back(x);
                        depth.by.push_back(y);
                        visited(x, y) = true;
                    }
                }

                // We have (+x,+y), fill the other 3 quadrants:
                // (-x,+y), (-x,-y), (-x,+y)
                vec1u idnew = uindgen(depth.bx.size());
                append(depth.bx, -depth.bx[idnew]);
                append(depth.by, depth.by[idnew]);
                append(depth.bx, -depth.bx[idnew]);
                append(depth.by, -depth.by[idnew]);
                append(depth.bx, depth.bx[idnew]);
                append(depth.by, -depth.by[idnew]);
            }
        };
    }
}

namespace astro {
    template<typename TypeR1, typename TypeD1, typename TypeR2, typename TypeD2>
    qxmatch_res qxmatch(const vec<1,TypeR1>& ra1, const vec<1,TypeD1>& dec1,
        const vec<1,TypeR2>& ra2, const vec<1,TypeD2>& dec2,
        qxmatch_params params = qxmatch_params{}) {

        qxmatch_res res;

        phypp_check(ra1.dims == dec1.dims, "first RA and Dec dimensions do not match (",
            ra1.dims, " vs ", dec1.dims, ")");
        phypp_check(ra2.dims == dec2.dims, "second RA and Dec dimensions do not match (",
            ra2.dims, " vs ", dec2.dims, ")");

        // Make sure we are dealing with finite coordinates...
        phypp_check(count(!is_finite(ra1) || !is_finite(dec1)) == 0,
            "first RA and Dec coordinates contain invalid values (infinite or NaN)");
        phypp_check(count(!is_finite(ra2) || !is_finite(dec2)) == 0,
            "second RA and Dec coordinates contain invalid values (infinite or NaN)");

        // Make sure that the requested number of matches is valid
        uint_t nth = clamp(params.nth, 1u, npos);

        const uint_t n1 = ra1.size();
        const uint_t n2 = ra2.size();

        res.id = replicate(npos, nth, n1);
        res.d  = replicate(dinf, nth, n1);

        if (!params.no_mirror) {
            res.rid = replicate(npos, n2);
            res.rd  = replicate(dinf, n2);
        }

        if (ra1.empty() || ra2.empty()) {
            return res;
        }

        // Convert input coordinates into proper units, if required
        vec1d dra1, ddec1, dcdec1, dra2, ddec2, dcdec2;
        if (params.linear) {
            // Not used
        } else {
            const double d2r = dpi/180.0;
            dra1  = ra1*d2r;
            ddec1 = dec1*d2r;
            dcdec1 = cos(ddec1);
            dra2  = ra2*d2r;
            ddec2 = dec2*d2r;
            dcdec2 = cos(ddec2);
        }

        auto true_distance = [&](double x1, double y1, double x2, double y2) {
            // Function to provide the distance between two points.
            if (params.linear) {
                return sqrt(sqr(x2 - x1) + sqr(y2 - y1));
            } else {
                return angdist(x1, y1, x2, y2);
            }
        };

        auto distance_proxy = [&](uint_t i, uint_t j) {
            // Function to provide a distance *indicator*.
            // Note that this is not the 'true' distance, but this it is sufficient
            // to find the nearest neighbors as long as this indicator is related to the
            // true distance through a monotonous function (that does not change ordering).
            if (params.linear) {
                return sqr(ra2.safe[j] - ra1.safe[i]) + sqr(dec2.safe[j] - dec1.safe[i]);
            } else {
                double sra = sin(0.5*(dra2.safe[j] - dra1.safe[i]));
                double sde = sin(0.5*(ddec2.safe[j] - ddec1.safe[i]));
                return sde*sde + sra*sra*dcdec2.safe[j]*dcdec1.safe[i];
            }
        };

        auto true_to_proxy = vectorize_lambda([&](double dist) {
            // Function to convert a "true" distance into the "proxy" distance
            // used by "distance_proxy".
            if (params.linear) {
                return sqr(dist);
            } else {
                return sqr(sin(dist/(3600.0*(180.0/dpi)*2)));
            }
        });
        auto proxy_to_true = vectorize_lambda([&](double dist) {
            // Function to convert a "proxy" distance into the "true" distance.
            if (params.linear) {
                return sqrt(dist);
            } else {
                return 3600.0*(180.0/dpi)*2*asin(sqrt(dist));
            }
        });

        if (!params.brute_force && ra2.size() < 10) {
            // The bucket algorithm was requested, but there are too few objects to
            // cross match to. The brute force algorithm is likely to be the fastest.
            params.brute_force = false;
        }

        if (!params.brute_force) {
            // Get bounds of the fields
            vec1d rra1  = {min(ra1),  max(ra1)};
            vec1d rra2  = {min(ra2),  max(ra2)};
            vec1d rdec1 = {min(dec1), max(dec1)};
            vec1d rdec2 = {min(dec2), max(dec2)};

            vec1d rra  = {std::min(rra1[0],  rra2[0]),  std::max(rra1[1],  rra2[1])};
            vec1d rdec = {std::min(rdec1[0], rdec2[0]), std::max(rdec1[1], rdec2[1])};

            // Choose a bucket size (arcsec)
            const double overgrowth = 10.0;
            uint_t nc2 = ceil(0.5*sqrt(dpi*ra2.size()/nth/overgrowth));
            vec1d hx = vec1d{rra2[0],  rra2[0],  rra2[1],  rra2[1]};
            vec1d hy = vec1d{rdec2[0], rdec2[1], rdec2[0], rdec2[1]};

            double cell_size = params.linear ? sqrt((rra2[1] - rra2[0])*(rdec2[1] - rdec2[0])) :
                std::max(3600.0*sqrt(field_area_hull(hx, hy))/nc2, 1.0);

            // Be careful that RA and Dec are spherical coordinates
            double dra, ddec;
            if (params.linear) {
                dra = cell_size;
                ddec = cell_size;
            } else {
                dra = cell_size/abs(cos(mean(rdec2)*dpi/180.0))/3600.0;
                ddec = cell_size/3600.0;
            }

            // Add some padding to prevent border issues
            rra[0]  -= dra;
            rra[1]  += dra;
            rdec[0] -= ddec;
            rdec[1] += ddec;

            // Final number of buckets
            uint_t nra  = (rra[1]  - rra[0])/dra;
            uint_t ndec = (rdec[1] - rdec[0])/ddec;

            // Build the buckets
            struct bucket_t {
                std::vector<uint_t> ids1;
                std::vector<uint_t> ids2;
            };

            vec<2,bucket_t> buckets(nra, ndec);

            // Fill the buckets
            vec1u idx1 = floor((ra1 - rra[0])/dra);
            vec1u idy1 = floor((dec1 - rdec[0])/ddec);
            for (uint_t i : range(ra1)) {
                buckets(idx1[i],idy1[i]).ids1.push_back(i);
            }

            vec1u idx2 = floor((ra2 - rra[0])/dra);
            vec1u idy2 = floor((dec2 - rdec[0])/ddec);
            for (uint_t j : range(ra2)) {
                buckets(idx2[j],idy2[j]).ids2.push_back(j);
            }

            // Precompute generic bucket geometry
            impl::qxmatch_impl::depth_cache depths(nra, ndec, cell_size);

            auto work1 = [&] (uint_t i, impl::qxmatch_impl::depth_cache& tdepths, qxmatch_res& tres) {
                int_t x0 = idx1[i];
                int_t y0 = idy1[i];

                // Compute the distance of this source from the bucket center
                double cell_dist = true_distance(ra1[i], dec1[i],
                    rra[0] + (x0+0.5)*dra, rdec[0] + (y0+0.5)*ddec);

                // Scan the buckets around this source until all neighbors are found
                uint_t d = 0;
                double reached_distance = 0.0;

                do {
                    auto& depth = tdepths[d];
                    for (uint_t b : range(depth.bx)) {
                        int_t x = x0+depth.bx.safe[b];
                        int_t y = y0+depth.by.safe[b];
                        if (x < 0 || uint_t(x) >= nra || y < 0 || uint_t(y) >= ndec) continue;

                        auto& bucket = buckets(uint_t(x),uint_t(y));
                        for (uint_t j : bucket.ids2) {
                            if (params.self && i == j) continue;

                            auto sd = distance_proxy(i, j);

                            // Compare this new distance to the largest one that is in the
                            // Nth nearest neighbor list. If it is lower than that, we
                            // insert it in the list, removing the old one, and sort the
                            // whole thing so that the largest distance goes as the end of
                            // the list.
                            if (sd < tres.d.safe(nth-1,i)) {
                                tres.id.safe(nth-1,i) = j;
                                tres.d.safe(nth-1,i) = sd;
                                uint_t k = nth-2;
                                while (k != npos && tres.d.safe(k,i) > tres.d.safe(k+1,i)) {
                                    std::swap(tres.d.safe(k,i), tres.d.safe(k+1,i));
                                    std::swap(tres.id.safe(k,i), tres.id.safe(k+1,i));
                                    --k;
                                }
                            }
                        }
                    }

                    reached_distance = true_to_proxy(
                        std::max(0.0, tdepths[d].max_dist - 1.1*cell_dist - 0.1*cell_size));

                    ++d;
                } while (tres.d.safe(nth-1,i) > reached_distance && d < tdepths.size());
            };

            auto work2 = [&] (uint_t j, impl::qxmatch_impl::depth_cache& tdepths, qxmatch_res& tres) {
                int_t x0 = idx2[j];
                int_t y0 = idy2[j];

                // Compute the distance of this source from the bucket center
                double cell_dist = true_distance(ra2[j], dec2[j],
                    rra[0] + (x0+0.5)*dra, rdec[0] + (y0+0.5)*ddec);

                // Scan the buckets around this source until all neighbors are found
                uint_t d = 0;
                double reached_distance = 0.0;

                do {
                    auto& depth = tdepths[d];
                    for (uint_t b : range(depth.bx)) {
                        int_t x = x0+depth.bx[b];
                        int_t y = y0+depth.by[b];

                        if (x < 0 || uint_t(x) >= nra || y < 0 || uint_t(y) >= ndec) continue;

                        auto& bucket = buckets(uint_t(x),uint_t(y));
                        for (uint_t i : bucket.ids1) {
                            double sd = distance_proxy(i, j);

                            // Just keep the nearest match
                            if (sd < tres.rd.safe[j]) {
                                tres.rd.safe[j] = sd;
                                tres.rid.safe[j] = i;
                            }
                        }
                    }

                    reached_distance = true_to_proxy(
                        std::max(0.0, tdepths[d].max_dist - 1.1*cell_dist - 0.1*cell_size));

                    ++d;
                } while (tres.rd[j] > reached_distance && d < tdepths.size());
            };

            if (params.thread <= 1) {
                if (n2 < nth) {
                    // We asked more neighbors than there are sources in the second catalog...
                    // Lower 'nth' to prevent this from blocking the algorithm.
                    nth = n2;
                }

                // When using a single thread, all the work is done in the main thread
                auto p = progress_start(n1+(!params.self && !params.no_mirror ? n2 : 0));
                for (uint_t i : range(ra1)) {
                    work1(i, depths, res);
                    if (params.verbose) progress(p, 31);
                }
                if (!params.self && !params.no_mirror) {
                    for (uint_t j : range(ra2)) {
                        work2(j, depths, res);
                        if (params.verbose) progress(p, 31);
                    }
                }
            } else {
                // When using more than one thread, the work load is evenly separated between
                // all the available threads, such that they should more or less all end at
                // the same time.
                std::atomic<uint_t> iter(0);

                // Create the thread pool and launch the threads
                std::vector<qxmatch_res> vres(params.thread);
                for (auto& r : vres) {
                    r.id = replicate(npos, nth, n1);
                    r.d  = dblarr(nth, n1)+dinf;
                    r.rid = replicate(npos, n2);
                    r.rd  = dblarr(n2)+dinf;
                }

                if (n2 < nth) {
                    // We asked more neighbors than there are sources in the second catalog...
                    // Lower 'nth' to prevent this from blocking the algorithm.
                    nth = n2;
                }

                vec1u tbeg1(params.thread);
                vec1u tend1(params.thread);
                vec1u tbeg2(params.thread);
                vec1u tend2(params.thread);
                auto pool = thread::pool(params.thread);
                uint_t total1 = 0;
                uint_t total2 = 0;
                uint_t assigned1 = floor(n1/float(params.thread));
                uint_t assigned2 = floor(n2/float(params.thread));
                for (uint_t t = 0; t < params.thread; ++t) {
                    if (t == params.thread-1) {
                        assigned1 = n1 - total1;
                        assigned2 = n2 - total2;
                    }

                    tbeg1[t] = total1;
                    tend1[t] = total1+assigned1;
                    tbeg2[t] = total2;
                    tend2[t] = total2+assigned2;

                    pool[t].start([&iter, &work1, &work2, &vres,
                        t, tbeg1, tbeg2, tend1, tend2, depths, params]() mutable {
                        for (uint_t i = tbeg1[t]; i < tend1[t]; ++i) {
                            work1(i, depths, vres[t]);
                            ++iter;
                        }

                        if (!params.self && !params.no_mirror) {
                            for (uint_t j = tbeg2[t]; j < tend2[t]; ++j) {
                                work2(j, depths, vres[t]);
                                ++iter;
                            }
                        }
                    });

                    total1 += assigned1;
                    total2 += assigned2;
                }

                // Wait for the computation to finish.
                // Here the main thread does nothing except sleeping, occasionally waking
                // up once in a while to update the progress bar if any.
                uint_t niter = n1+(!params.self && !params.no_mirror ? n2 : 0);
                auto p = progress_start(niter);
                while (iter < niter) {
                    thread::sleep_for(0.2);
                    if (params.verbose) print_progress(p, iter);
                }

                // By now, all threads should have ended their tasks.
                // We must ask them to terminate nicely.
                for (auto& t : pool) {
                    t.join();
                }

                // Merge back the results of each thread
                for (uint_t t = 0; t < params.thread; ++t) {
                    auto ids = tbeg1[t]-_-(tend1[t]-1);
                    res.id.safe(_,ids) = vres[t].id.safe(_,ids);
                    res.d.safe(_,ids) = vres[t].d.safe(_,ids);

                    if (!params.no_mirror) {
                        ids = tbeg2[t]-_-(tend2[t]-1);
                        res.rid.safe[ids] = vres[t].rid.safe[ids];
                        res.rd.safe[ids] = vres[t].rd.safe[ids];
                    }
                }
            }
        } else {
            auto work = [&, nth] (uint_t i, uint_t j, qxmatch_res& tres) {
                double sd = distance_proxy(i, j);

                // We compare this new distance to the largest one that is in the Nth
                // nearest neighbor list. If it is lower than that, we insert it in the list,
                // removing the old one, and sort the whole thing so that the largest distance
                // goes as the end of the list.
                if (sd < tres.d.safe(nth-1,i)) {
                    tres.id.safe(nth-1,i) = j;
                    tres.d.safe(nth-1,i) = sd;
                    uint_t k = nth-2;
                    while (k != npos && tres.d.safe(k,i) > tres.d.safe(k+1,i)) {
                        std::swap(tres.d.safe(k,i), tres.d.safe(k+1,i));
                        std::swap(tres.id.safe(k,i), tres.id.safe(k+1,i));
                        --k;
                    }
                }

                // Take care of reverse search
                if (!params.no_mirror && sd < tres.rd.safe[j]) {
                    tres.rid.safe[j] = i;
                    tres.rd.safe[j] = sd;
                }
            };

            if (params.thread <= 1) {
                // When using a single thread, all the work is done in the main thread
                auto p = progress_start(n1);
                for (uint_t i = 0; i < n1; ++i) {
                    for (uint_t j = 0; j < n2; ++j) {
                        if (params.self && i == j) continue;
                        work(i,j,res);
                    }

                    if (params.verbose) progress(p);
                }
            } else {
                // When using more than one thread, the work load is evenly separated between
                // all the available threads, such that they should more or less all end at
                // the same time.
                std::atomic<uint_t> iter(0);

                // Create the thread pool and launch the threads
                std::vector<qxmatch_res> vres(params.thread);
                for (auto& r : vres) {
                    r.id = replicate(npos, nth, n1);
                    r.d  = replicate(dinf, nth, n1);
                    if (!params.no_mirror) {
                        r.rid = replicate(npos, n2);
                        r.rd  = replicate(dinf, n2);
                    }
                }

                vec1u tbeg(params.thread);
                vec1u tend(params.thread);
                auto pool = thread::pool(params.thread);
                uint_t total = 0;
                uint_t assigned = floor(n1/float(params.thread));
                for (uint_t t = 0; t < params.thread; ++t) {
                    if (t == params.thread-1) {
                        assigned = n1 - total;
                    }

                    tbeg[t] = total;
                    tend[t] = total+assigned;

                    pool[t].start([&iter, &work, &vres, t, tbeg, tend, params, n2]() {
                        for (uint_t i = tbeg[t]; i < tend[t]; ++i) {
                            for (uint_t j = 0; j < n2; ++j) {
                                if (params.self && i == j) continue;
                                work(i,j,vres[t]);
                            }

                            ++iter;
                        }
                    });

                    total += assigned;
                }

                // Wait for the computation to finish.
                // Here the main thread does nothing except sleeping, occasionally waking up
                // once in a while to update the progress bar if any.
                auto p = progress_start(n1);
                while (iter < n1) {
                    thread::sleep_for(0.2);
                    if (params.verbose) print_progress(p, iter);
                }

                // By now, all threads should have ended their tasks.
                // We must ask them to terminate nicely.
                for (auto& t : pool) {
                    t.join();
                }

                // Merge back the results of each thread
                for (uint_t t = 0; t < params.thread; ++t) {
                    auto ids = tbeg[t]-_-(tend[t]-1);
                    res.id.safe(_,ids) = vres[t].id.safe(_,ids);
                    res.d.safe(_,ids) = vres[t].d.safe(_,ids);

                    if (!params.no_mirror) {
                        if (t == 0) {
                            res.rid = vres[t].rid;
                            res.rd = vres[t].rd;
                        } else {
                            for (uint_t j = 0; j < n2; ++j) {
                                if (res.rd.safe[j] >= vres[t].rd.safe[j]) {
                                    res.rid.safe[j] = vres[t].rid.safe[j];
                                    res.rd.safe[j] = vres[t].rd.safe[j];
                                }
                            }
                        }
                    }
                }
            }
        }

        // Convert the distance estimator to a real distance
        res.d = proxy_to_true(res.d);
        if (!params.no_mirror) {
            res.rd = proxy_to_true(res.rd);
        }

        return res;
    }

    template<typename C1, typename C2,
        typename enable = typename std::enable_if<!meta::is_vec<C1>::value>::type>
    qxmatch_res qxmatch(const C1& cat1, const C2& cat2, qxmatch_params params = qxmatch_params{}) {
        return qxmatch(cat1.ra, cat1.dec, cat2.ra, cat2.dec, params);
    }

    template<typename C1, typename enable = typename std::enable_if<!meta::is_vec<C1>::value>::type>
    qxmatch_res qxmatch(const C1& cat1, qxmatch_params params = qxmatch_params{}) {
        return qxmatch(cat1.ra, cat1.dec, params);
    }

    template<typename TypeR1, typename TypeD1>
    qxmatch_res qxmatch(const vec<1,TypeR1>& ra1, const vec<1,TypeD1>& dec1,
        qxmatch_params params = qxmatch_params{}) {
        params.self = true;
        return qxmatch(ra1, dec1, ra1, dec1, params);
    }

    struct id_pair {
        vec1u id1, id2;
        vec1u lost;
    };

    id_pair xmatch_clean_best(const qxmatch_res& r) {
        const uint_t ngal = dim(r.id)[1];

        id_pair c;
        c.id1.data.reserve(ngal);
        c.id2.data.reserve(ngal);
        c.lost.data.reserve(ngal/6);

        for (uint_t i = 0; i < ngal; ++i) {
            if (r.rid[r.id[i]] == i) {
                c.id1.data.push_back(i);
                c.id2.data.push_back(r.id[i]);
            } else {
                c.lost.data.push_back(i);
            }
        }

        c.id1.dims[0] = c.id1.data.size();
        c.id2.dims[0] = c.id2.data.size();
        c.lost.dims[0] = c.lost.data.size();

        return c;
    }

    void xmatch_check_lost(const id_pair& p) {
        if (!p.lost.empty()) {
            warning(p.lost.size(), " sources failed to cross match");
        }
    }

    void xmatch_save_lost(const id_pair& p, const std::string& save) {
        if (!p.lost.empty()) {
            warning(p.lost.size(), " sources failed to cross match");
            fits::write_table(save, ftable(p.lost));
        }
    }

    struct qdist_params {
        bool verbose = false;
    };

    template<typename TypeR, typename TypeD>
    vec2d qdist(const vec<1,TypeR>& ra, const vec<1,TypeD>& dec,
        qdist_params params = qdist_params{}) {
        vec2d ret(ra.size(), ra.size());

        phypp_check(ra.dims == dec.dims, "first RA and Dec dimensions do not match (",
            ra.dims, " vs ", dec.dims, ")");

        const double d2r = dpi/180.0;
        auto dra  = ra*d2r;
        auto ddec = dec*d2r;
        auto dcdec = cos(ddec);

        const uint_t n = ra.size();

        auto distance = [&](uint_t i, uint_t j) {
            double sra = sin(0.5*(dra.safe[j] - dra.safe[i]));
            double sde = sin(0.5*(ddec.safe[j] - ddec.safe[i]));
            double d = sde*sde + sra*sra*dcdec.safe[j]*dcdec.safe[i];
            return 3600.0*(180.0/dpi)*2*asin(sqrt(d));
        };

        // When using a single thread, all the work is done in the main thread
        auto p = progress_start(n*(n-1)/2);
        for (uint_t i : range(n)) {
            for (uint_t j : range(i+1, n)) {
                ret(j,i) = distance(i, j);

                if (params.verbose) progress(p, 113);
            }
        }

        return ret;
    }
}
}

#endif
