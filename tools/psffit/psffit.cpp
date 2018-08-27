#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string psf_model;
    std::string serr;
    float radius = 0.0;
    float psf_frac = 1.0;
    float fconv = 1.0;
    bool silent = false;
    bool nobg = false;
    bool save = false;
    std::string resi;

    read_args(argc-3, argv+3, arg_list(
        name(psf_model, "psf"), radius, psf_frac, name(serr, "error"), fconv, silent,
        name(resi, "residual"), nobg, save
    ));

    vec2d img, err;
    fits::read(argv[1], img);
    if (!serr.empty()) {
        fits::read(serr, err);
    } else {
        if (!silent) warning("no error map provided, reported errors will be wrong");
        err = img*0 + 1;
    }

    float fx0, fy0;
    from_string(argv[2], fx0);
    from_string(argv[3], fy0);
    int_t x0 = round(fx0), y0 = round(fy0);

    vec2d psf;
    if (!make_psf({{img.dims[0], img.dims[1]}}, x0, y0, psf_model, psf)) {
        return 1;
    }

    vec1u idf;
    if (psf_frac != 1.0) {
        idf = where(is_finite(img) && is_finite(err) && is_finite(psf) &&
            psf > (1.0 - psf_frac)*max(psf[where(is_finite(psf))]));
    } else if (radius != 0.0) {
        vec2d rad = generate_img({{img.dims[0], img.dims[1]}}, [=](uint_t x, uint_t y) {
            return sqr(double(x) - x0) + sqr(double(y) - y0);
        });

        idf = where(is_finite(img) && is_finite(err) && is_finite(psf) && rad < sqr(radius));
    } else {
        idf = where(is_finite(img) && is_finite(err) && is_finite(psf));
    }

    if (nobg) {
        auto res = linfit(img[idf], err[idf], psf[idf]);

        if (!silent) print("flux: ", fconv*res.params[0], " +/- ", fconv*res.errors[0], " (flux unit)");

        if (!resi.empty()) {
            vec2d ires = img - res.params[0]*psf;
            fits::write(resi, ires);
        }

        if (save) {
            std::string outname = argv[1];
            outname = outname.substr(0, outname.size() - length(".fits")) + "_flx.fits";
            fits::write_table(outname, "flux", fconv*res.params[0], "flux_err", fconv*res.errors[0]);
        }
    } else {
        auto res = linfit(img[idf], err[idf], 1.0, psf[idf]);

        if (!silent) {
            print("flux: ", fconv*res.params[1], " +/- ", fconv*res.errors[1], " (flux unit)");
            print("background: ", res.params[0], " +/- ", res.errors[0], " (map unit)");
        }

        if (!resi.empty()) {
            if (resi == "1") {
                std::string outname = argv[1];
                resi = outname.substr(0, outname.size() - length(".fits")) + "_res.fits";
            }

            vec2d ires = img - res.params[1]*psf;
            fits::write(resi, ires);
        }

        if (save) {
            std::string outname = argv[1];
            outname = outname.substr(0, outname.size() - length(".fits")) + "_flx.fits";
            fits::write_table(outname,
                "flux", fconv*res.params[1], "flux_err", fconv*res.errors[1],
                "bg", fconv*res.params[0], "bg_err", fconv*res.errors[0]
            );
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;
    print("psffit v1.0");
    paragraph("usage: psffit img.fits x y psf=...");
}
