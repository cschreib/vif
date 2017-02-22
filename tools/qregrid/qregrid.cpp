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
    double pixfrac = 1.0;
    bool approx = false;
    std::string weight;
    std::string method = "drizzle";
    bool conserve_flux = false;
    read_args(argc-2, argv+2, arg_list(
        verbose, name(tpl, "template"), aspix, ratio, method, conserve_flux, weight, pixfrac, approx
    ));

    // Forward options
    astro::regrid_interpolate_params iopts;
    astro::regrid_drizzle_params     dopts;
    iopts.verbose = verbose;
    dopts.verbose = verbose;

    if (method == "drizzle") {
        dopts.pixfrac = pixfrac;
        if (approx) {
            dopts.dest_pixfrac = true;
        }
    } else if (method == "nearest") {
        iopts.conserve_flux = conserve_flux;
        iopts.method = astro::interpolation_method::nearest;
    } else if (method == "linear") {
        iopts.conserve_flux = conserve_flux;
        iopts.method = astro::interpolation_method::linear;
    } else {
        error("unknown regridding method '", method, "'");
        return 1;
    }

    // Read source image
    fits::input_image fimgs(img_src_file);
    vec2d imgs;
    fimgs.read(imgs);

    // Define input and output grids
    astro::wcs astros;  // input astrometry
    astro::wcs astrod;  // output astrometry
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
        params.pixel_ref_x = imgs.dims[1]/2 + 1;
        params.pixel_ref_y = imgs.dims[0]/2 + 1;

        fits::header hdrs;
        astro::make_wcs_header(params, hdrs);
        astros = astro::wcs(hdrs);

        vec1u dims = max(ceil(fimgs.image_dims()*ratio), 1);
        params.pixel_scale = 1/ratio;
        params.dims_x = dims[1];
        params.dims_y = dims[0];
        params.pixel_ref_x = dims[1]/2 + 1;
        params.pixel_ref_y = dims[0]/2 + 1;
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
            vec1u dims = max(ceil(fimgs.image_dims()/ratio), 1);

            astro::make_wcs_header_params params;
            params.pixel_scale = aspix;
            params.dims_x = dims[1];
            params.dims_y = dims[0];
            astro::xy2ad(astros, imgs.dims[0]/2 + 1, imgs.dims[1]/2 + 1, params.sky_ref_ra, params.sky_ref_dec);
            params.pixel_ref_x = dims[1]/2 + 1;
            params.pixel_ref_y = dims[0]/2 + 1;

            astro::make_wcs_header(params, hdrd);
            astrod = astro::wcs(hdrd);
        } else {
            error("must specify output grid with ratio=..., aspix=..., or template=...");
            return 1;
        }
    }

    // Regrid
    vec2d wei, res;
    if (method == "drizzle") {
        res = astro::regrid_drizzle(imgs, astros, astrod, wei, dopts);
    } else {
        res = astro::regrid_interpolate(imgs, astros, astrod, iopts);
    }

    // Save regridded image
    file::mkdir(file::get_directory(out_file));
    fits::write(out_file, res, hdrd);

    if (!weight.empty()) {
        // Save weights
        file::mkdir(file::get_directory(weight));
        fits::write(weight, wei, hdrd);
    }

    return 0;
}
