#include <phypp.hpp>

void print_help();

bool get_flux(int argc, char* argv[], const vec3d& cube);
bool get_logdisp(int argc, char* argv[], const vec3d& cube);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string file = argv[1];
    if (!file::exists(file)) {
        error("could not find '"+file+"'");
        return 1;
    }

    std::string op = argv[2];

    if (op == "flux") {
        vec3d cube; fits::read(file, cube);
        return get_flux(argc-2, argv+2, cube) ? 0 : 1;
    } else if (op == "logdisp") {
        vec3d cube; fits::read(file, cube);
        return get_logdisp(argc-2, argv+2, cube) ? 0 : 1;
    } else if (op == "batch") {
        return run_batch(argc-2, argv+2, file) ? 0 : 1;
    } else {
        error("unknown operation '"+op+"'");
        return 1;
    }
}

struct flux_extractor {
    bool mea = true;
    uint_t nbstrap = 200;
    uint_t tseed = 42;

    vec2d psf;
    vec2b ipsf;

    bool config(program_arguments& pa) {
        std::string tpsf;
        bool med = false;
        bool norm = false;
        double frac = 0.1;

        pa.read(arg_list(
            name(tpsf, "psf"), name(mea, "mean"), name(med, "median"), frac, norm, nbstrap,
            name(tseed, "seed")
        ));

        if (med) mea = false;

        if (tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        fits::read(tpsf, psf);

        if (norm) {
            psf /= max(fabs(psf));
        }

        ipsf = fabs(psf) > frac*max(fabs(psf));
        vec1u id = where(ipsf);
        if (id.empty()) {
            error("no pixel to fit (frac=", frac,")");
            return false;
        }

        return true;
    }

    vec2d img;

    void extract(const vec3d& cube, double& flux, double& bg) {
        if (mea) {
            img = qstack_mean(cube);
        } else {
            img = qstack_median(cube);
        }

        vec1u id = where(ipsf && finite(img));
        if (id.empty()) {
            flux = dnan;
            bg = dnan;
            return;
        }

        auto res = linfit(img[id], 1.0, 1.0, psf[id]);
        flux = res.params[1];
        bg = res.params[0];
    }

    void bootstrap(const vec3d& cube, double& flux, double& bg) {
        auto seed = make_seed(tseed);

        vec1d fluxes, bgs;
        uint_t nsrc = cube.dims[0];
        qstack_bootstrap(cube, nbstrap, nsrc/2, seed, [&](const vec3d& sub) {
            double f, b;
            extract(sub, f, b);
            fluxes.push_back(f);
            bgs.push_back(b);
        });

        fluxes = fluxes[where(finite(fluxes))];
        bgs = bgs[where(finite(bgs))];

        flux = stddev(fluxes)/sqrt(2);
        bg = stddev(bgs)/sqrt(2);
    }
};

struct logdisp_extractor {
    uint_t nbstrap = 200;
    uint_t tseed = 42;

    vec2d psf;
    vec2b ipsf, ifpsf;

    bool config(program_arguments& pa) {
        std::string tpsf;
        bool norm = false;
        double frac = 0.5;
        double ffrac = 0.1;

        pa.read(arg_list(
            name(tpsf, "psf"), frac, ffrac, norm, nbstrap, name(tseed, "seed")
        ));

        if (tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        fits::read(tpsf, psf);

        if (norm) {
            psf /= max(fabs(psf));
        }

        ipsf = fabs(psf) > frac*max(fabs(psf));
        vec1u id = where(ipsf);
        if (id.empty()) {
            error("no pixel to fit (frac=", frac,")");
            return false;
        }

        ifpsf = fabs(psf) > ffrac*max(fabs(psf));
        id = where(ifpsf);
        if (id.empty()) {
            error("no pixel to fit (frac=", ffrac,")");
            return false;
        }

        return true;
    }

    vec2d med, dsp;

    void extract(const vec3d& cube, double& disp, double& bg) {
        med.resize(cube.dims[1], cube.dims[2]);
        dsp.resize(cube.dims[1], cube.dims[2]);

        run_dim_idx(cube, 0, [&](uint_t i, vec1d& d) {
            med[i] = fast_median(d);
            for (uint_t j : range(d)) {
                d[j] = fabs(d[j] - med[i]);
            }
            dsp[i] = fast_median(d);
        });

        vec1u idd = where(ipsf && finite(dsp));
        if (idd.empty()) {
            disp = dnan;
            bg = dnan;
            return;
        }

        auto res = linfit(sqr(dsp[idd]), 1.0, 1.0, sqr(psf[idd]));
        disp = sqrt(res.params[1]);
        bg = 1.483*sqrt(res.params[0]);

        vec1u idf = where(ifpsf && finite(med));
        if (idf.empty()) {
            disp = dnan;
            bg = dnan;
            return;
        }

        res = linfit(med[idf], 1.0, 1.0, psf[idf]);
        disp /= res.params[1];

        disp = (1.171/disp)*(1.0 - sqrt(1.0 - sqr(disp/0.953)));
    }

    void bootstrap(const vec3d& cube, double& disp, double& bg) {
        auto seed = make_seed(tseed);

        vec1d disps, bgs;
        uint_t nsrc = cube.dims[0];
        qstack_bootstrap(cube, nbstrap, nsrc/2, seed, [&](const vec3d& sub) {
            double d, b;
            extract(sub, d, b);
            disps.push_back(d);
            bgs.push_back(b);
        });

        disps = disps[where(finite(disps))];
        bgs = bgs[where(finite(bgs))];

        disp = stddev(disps)/sqrt(2);
        bg = stddev(bgs)/sqrt(2);
    }
};

bool get_flux(int argc, char* argv[], const vec3d& cube) {
    bool bstrap = false;
    std::string out = "";

    flux_extractor ex;

    {
        program_arguments pa(argc, argv);
        pa.read(arg_list(bstrap, out));
        if (!ex.config(pa)) return false;
    }

    print("nsrc: ", cube.dims[0]);

    double flux, bg;
    ex.extract(cube, flux, bg);
    if (bstrap) {
        double flux_err, bg_err;
        ex.bootstrap(cube, flux_err, bg_err);
        print("flux: ", flux, " +/- ", flux_err);
        print("background: ", bg, " +/- ", bg_err);
    } else {
        print("flux: ", flux);
        print("background: ", bg);
    }

    if (!out.empty()) {
        fits::write(out+".fits", ex.img);
    }

    return true;
}

bool get_logdisp(int argc, char* argv[], const vec3d& cube) {
    bool bstrap = false;
    std::string out = "";

    logdisp_extractor ex;

    {
        program_arguments pa(argc, argv);
        pa.read(arg_list(bstrap, out));
        if (!ex.config(pa)) return false;
    }

    print("nsrc: ", cube.dims[0]);

    double disp, bg;
    ex.extract(cube, disp, bg);
    if (bstrap) {
        double disp_err, bg_err;
        ex.bootstrap(cube, disp_err, bg_err);
        print("logdisp: ", disp, " +/- ", disp_err);
        print("background: ", bg, " +/- ", bg_err);
    } else {
        print("logdisp: ", disp);
        print("background: ", bg);
    }

    if (!out.empty()) {
        fits::write(out+".fits", ex.dsp);
    }

    return true;
}

void print_help() {
    using namespace format;

    print("fluxcube v1.0");
    header("Usage: fluxcube cube.fits [options]");
    header("Available operations:");
    bullet("flux", "perform mean or median stacking then PSF fitting to get the flux");
    bullet("logdisp", "perform MAD stacking then PSF fitting to get the lognormal dispersion");

    print("");
    paragraph("Copyright (c) 2014 C. Schreiber (corentin.schreiber@cea.fr)");

    paragraph("This software is provided 'as-is', without any express or implied warranty. In no "
        "event will the authors be held liable for any damages arising from the use of this "
        "software.");

    paragraph("Permission is granted to anyone to use this software for any purpose, including "
        "commercial applications, and to alter it and redistribute it freely, subject to the "
        "following restrictions:");

    bullet("1", "The origin of this software must not be misrepresented; you must not claim that "
        "you wrote the original software. If you use this software in a product, an acknowledgment "
        "in the product documentation would be appreciated but is not required.");
    bullet("2", "Altered source versions must be plainly marked as such, and must not be "
        "misrepresented as being the original software.");
    bullet("3", "This notice may not be removed or altered from any source distribution.");

    print("");
}
