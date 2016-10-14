#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 3) {
        return 0;
    }

    std::string img_src_file = argv[1];
    std::string out_file = argv[2];

    std::string tpl;
    bool verbose = false;
    double aspix = dnan;
    double ratio = dnan;
    read_args(argc-2, argv+2, arg_list(verbose, name(tpl, "template"), aspix, ratio));

    // Read source image
    fits::input_image fimgs(img_src_file);
    vec2d imgs;
    fimgs.read(imgs);

    // Define input and output grids
    astro::wcs astros;  // input astrometry
    astro::wcs astrod;  // output astrometry
    vec1u dims;        // dimensions of output image
    fits::header hdrd; // FITS header of output image (with astrometry)

    bool has_wcs = fimgs.has_keyword("CTYPE1");
    if (has_wcs) {
        astros = astro::wcs(fimgs.read_header());
        if (!astros.is_valid()) {
            has_wcs = false;
        }
    }

    if (is_finite(ratio) && has_wcs) {
        // The image has WCS information, transform the size ratio into a pixel scale change
        double aspix_orig;
        if (!astro::get_pixel_size(astros, aspix_orig)) {
            error("missing or incorrect WCS data in source image");
            return 1;
        }

        aspix = aspix_orig/ratio;
        ratio = dnan;
    }

    if (is_finite(ratio)) {
        // Simple physical rescaling
        astro::make_wcs_header_params params;
        params.pixel_scale = 1.0;
        params.dims_x = imgs.dims[1];
        params.dims_y = imgs.dims[0];
        params.sky_ref_ra = 0.0;
        params.sky_ref_dec = 0.0;
        params.pixel_ref_x = imgs.dims[1]/2;
        params.pixel_ref_y = imgs.dims[0]/2;

        fits::header hdrs;
        astro::make_wcs_header(params, hdrs);
        astros = astro::wcs(hdrs);

        dims = max(ceil(fimgs.image_dims()*ratio), 1);
        params.pixel_scale = 1/ratio;
        astro::make_wcs_header(params, hdrd);
        astrod = astro::wcs(hdrd);
    } else {
        // Read astrometry of input image
        if (!has_wcs) {
            error("missing WCS data in source image");
            return 1;
        }

        if (!tpl.empty()) {
            // Read output grid from template image
            fits::input_image fimgd(tpl);
            if (!fimgd.has_keyword("CTYPE1")) {
                error("missing WCS data in template image");
                return 1;
            }

            dims = fimgd.image_dims();
            hdrd = fimgd.read_header();
            astrod = astro::wcs(hdrd);
            if (!astrod.is_valid()) {
                error("invalid WCS data in template image");
                return 1;
            }
        } else if (is_finite(aspix)) {
            // Make output grid from input grid, same center, just different pixel scale
            double aspix_orig;
            if (!astro::get_pixel_size(astros, aspix_orig)) {
                error("missing or incorrect WCS data in source image");
                return 1;
            }

            ratio = aspix/aspix_orig;
            dims = max(ceil(fimgs.image_dims()/ratio), 1);

            astro::make_wcs_header_params params;
            params.pixel_scale = aspix;
            params.dims_x = dims[1];
            params.dims_y = dims[0];
            astro::xy2ad(astros, imgs.dims[0]/2 + 1, imgs.dims[1]/2 + 1, params.sky_ref_ra, params.sky_ref_dec);
            params.pixel_ref_x = dims[1]/2;
            params.pixel_ref_y = dims[0]/2;

            astro::make_wcs_header(params, hdrd);
            astrod = astro::wcs(hdrd);
        } else {
            error("must specify output grid with ratio=..., aspix=..., or template=...");
            return 1;
        }
    }

    // Regridded image
    vec2d res(dims[0], dims[1]);

    // Convenience functions
    auto getmin = [](const vec1d& v) {
        double tmp = floor(min(v));
        return tmp > 0 ? uint_t(tmp) : uint_t(0);
    };
    auto getmax = [](const vec1d& v, uint_t n) {
        double tmp = ceil(max(v));
        return tmp > double(n) ? n : uint_t(tmp);
    };

    // Function to find if a point lies on the right side of a polygon's edge
    // orient: orientation of the polygon
    // cx1, cy1, cx2, cy2: X and Y coordinates of the two nodes of the edge
    // x, y: X and Y coordinates of the point to test
    // return: true if the point is on the right side, false otherwise
    auto in_poly_edge = [](int_t orient, double cx1, double cy1, double cx2, double cy2,
        double x, double y) {

        double cross = (cx2 - cx1)*(y - cy1) - (cy2 - cy1)*(x - cx1);
        return cross*orient < 0;
    };

    // Function to find the position and existence of the intersection of two lines
    // l1x1, l1y1, l1x2, l1y2: X and Y coordinates of two points defining the first line
    // l2x1, l2y1, l2x2, l2y2: X and Y coordinates of two points defining the second line
    // x, y: output X and Y coordinates of the intersection point (if any)
    // return: true if intersection exists, false otherwise
    auto segment_intersect = [](double l1x1, double l1y1, double l1x2, double l1y2,
                                double l2x1, double l2y1, double l2x2, double l2y2,
                                double& x, double& y) {
        // Find the intersection point
        // http://stackoverflow.com/a/1968345/1565581
        double s1x = l1x2 - l1x1;
        double s1y = l1y2 - l1y1;
        double s2x = l2x2 - l2x1;
        double s2y = l2y2 - l2y1;

        double det = s1x*s2y - s1y*s2x;
        if (abs(det) < 5*std::numeric_limits<double>::epsilon()) {
            // 'det' is zero: the lines are parallel, no intersection
            return false;
        }

        double s12x = l1x1 - l2x1;
        double s12y = l1y1 - l2y1;
        double t = (s2x*s12y - s2y*s12x)/det;
        x = l1x1 + t*s1x;
        y = l1y1 + t*s1y;

        return true;
    };

    // Precompute the projection of the new pixel grid on the old
    auto pg = progress_start(res.size());
    vec1d plx(res.dims[1]+1);
    vec1d ply(res.dims[1]+1);
    for (uint_t ix : range(res.dims[1]+1)) {
        double tra, tdec;
        astro::xy2ad(astrod, ix+0.5, 0.5, tra, tdec);
        astro::ad2xy(astros, tra, tdec, plx.safe[ix], ply.safe[ix]);
        plx.safe[ix] -= 1.0; ply.safe[ix] -= 1.0;
    }

    for (uint_t iy : range(res.dims[0])) {
        vec1d pux(res.dims[1]+1);
        vec1d puy(res.dims[1]+1);
        for (uint_t ix : range(res.dims[1]+1)) {
            double tra, tdec;
            astro::xy2ad(astrod, ix+0.5, iy+1.5, tra, tdec);
            astro::ad2xy(astros, tra, tdec, pux.safe[ix], puy.safe[ix]);
            pux.safe[ix] -= 1.0; puy.safe[ix] -= 1.0;
        }

        for (uint_t ix : range(res.dims[1])) {
            // Find projection of each pixel of the new grid on the original image
            // NB: assumes the astrometry is such that this projection is
            // reasonably approximated by a 4-edge polygon (i.e.: varying pixel scales,
            // pixel offsets and rotations are fine, but weird things may happen close
            // to the poles of the projection where things become non-linear)

            vec1d xps = {plx.safe[ix], plx.safe[ix+1], pux.safe[ix+1], pux.safe[ix]};
            vec1d yps = {ply.safe[ix], ply.safe[ix+1], puy.safe[ix+1], puy.safe[ix]};

            // Get bounds of this projection
            uint_t ymin = getmin(yps-0.5), ymax = getmax(yps+0.5, imgs.dims[0]-1);
            uint_t xmin = getmin(xps-0.5), xmax = getmax(xps+0.5, imgs.dims[1]-1);

            // Get polygon orientation
            int_t orient = ((xps.safe[1] - xps.safe[0])*(yps.safe[2] - yps.safe[1]) -
                (yps.safe[1] - yps.safe[0])*(xps.safe[2] - xps.safe[1]) > 0) ? 1 : -1;

            // Sum flux from original pixels that fall inside the projection
            double flx = 0.0;
            for (uint_t ipy = ymin; ipy <= ymax; ++ipy)
            for (uint_t ipx = xmin; ipx <= xmax; ++ipx) {
                // Construct the intersection polygon of the original pixel and the projection
                // https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
                vec1d cpx = {ipx-0.5, ipx+0.5, ipx+0.5, ipx-0.5};
                vec1d cpy = {ipy-0.5, ipy-0.5, ipy+0.5, ipy+0.5};

                uint_t c2 = xps.size()-1;
                for (uint_t c1 : range(xps)) {
                    if (cpx.empty()) break;

                    vec1d icpx = cpx, icpy = cpy;
                    cpx.clear(); cpy.clear();

                    // Find out if all the current polygon's points are "inside" the projection's current edge
                    vec1b in(icpx.dims);
                    for (uint_t i : range(icpx)) {
                        in.safe[i] = in_poly_edge(orient,
                            xps.safe[c1], yps.safe[c1], xps.safe[c2], yps.safe[c2],
                            icpx.safe[i], icpy.safe[i]
                        );
                    }

                    uint_t i2 = icpx.size()-1;
                    for (uint_t i1 : range(icpx)) {
                        if (in.safe[i2] != in.safe[i1]) {
                            // This edge [i2-i1] is intersected by the projection's current edge,
                            // find the intersection point and add it to the polygon
                            double tx, ty;
                            segment_intersect(
                                xps.safe[c1], yps.safe[c1], xps.safe[c2], yps.safe[c2],
                                icpx.safe[i1], icpy.safe[i1], icpx.safe[i2], icpy.safe[i2],
                                tx, ty
                            );

                            cpx.push_back(tx);
                            cpy.push_back(ty);
                        }

                        if (in.safe[i1]) {
                            // The point i1 is "inside" the projection's current edge, keep it for now
                            cpx.push_back(icpx.safe[i1]);
                            cpy.push_back(icpy.safe[i1]);
                        }

                        i2 = i1;
                    }

                    c2 = c1;
                }

                // No intersection, just discard that pixel
                if (cpx.size() < 3) continue;

                // Compute the area of this intersection (1: full coverage, 0: no coverage)
                double area = 0.0;
                uint_t i3 = cpx.size()-2;
                uint_t i2 = cpx.size()-1;
                for (uint_t i1 : range(cpx.size()-2)) {
                    area += 0.5*abs(
                        cpx.safe[i1]*(cpy.safe[i2]-cpy.safe[i3]) +
                        cpx.safe[i2]*(cpy.safe[i3]-cpy.safe[i1]) +
                        cpx.safe[i3]*(cpy.safe[i1]-cpy.safe[i2])
                    );

                    i3 = i2;
                    i2 = i1;
                }


                flx += imgs.safe(ipy,ipx)*area;
            }

            res.safe(iy,ix) = flx;

            if (verbose) progress(pg, 31);
        }

        plx = std::move(pux);
        ply = std::move(puy);
    }

    // Save regridded image
    file::mkdir(file::get_directory(out_file));
    fits::write(out_file, res, hdrd);

    return 0;
}
