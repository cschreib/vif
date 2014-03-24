#include <phypp.hpp>

void print_help();

bool get_flux(int argc, char* argv[], const vec3d& cube);
bool get_logdisp(int argc, char* argv[], const vec3d& cube);
bool run_batch(int argc, char* argv[], const std::string& file);

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
    bool large_bg = false;
    bool residual = false;
    double frac = 0.1;

    vec2d psf_orig;
    vec2d psf;
    vec2b ipsf;

    bool update_masks_() {
        ipsf = fabs(psf) >= frac*max(fabs(psf));
        vec1u id = where(ipsf);
        if (id.empty()) {
            error("no pixel to fit (frac=", frac,")");
            return false;
        }

        return true;
    }

    bool config(program_arguments& pa, bool init = false) {
        std::string tpsf;
        bool med = false;
        bool norm = false;

        pa.read(arg_list(
            name(tpsf, "psf"), name(mea, "mean"), name(med, "median"), frac, norm, nbstrap,
            name(tseed, "seed"), large_bg, residual
        ));

        if (med) mea = false;

        if (!init && psf.empty() && tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        if (!tpsf.empty()) {
            fits::read(tpsf, psf);

            if (norm) {
                psf /= max(fabs(psf));
            }

            psf_orig = psf;

            if (!update_masks_()) {
                return false;
            }
        }

        return true;
    }

    vec2d img, res;

    void extract(const vec3d& cube, double& flux, double& bg) {
        if (psf.dims[0] != cube.dims[1] || psf.dims[1] != cube.dims[2]) {
            int_t hxsize = cube.dims[1]/2;
            int_t hysize = cube.dims[2]/2;
            vec1i mid = psf_orig.ids(max_id(psf_orig));
            psf = subregion(psf_orig, {mid[0]-hxsize, mid[1]-hysize, mid[0]+hxsize, mid[1]+hysize});

            if (!update_masks_()) {
                flux = dnan;
                bg = dnan;
                return;
            }
        }

        if (mea) {
            img = qstack_mean(cube);
        } else {
            img = qstack_median(cube);
        }

        if (large_bg) {
            vec1u id = where(finite(img));
            if (id.empty()) {
                flux = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            auto r = linfit(img[id], 1.0, 1.0, psf[id]*ipsf[id]);
            flux = r.params[1];
            bg = r.params[0];
        } else {
            vec1u id = where(ipsf && finite(img));
            if (id.empty()) {
                flux = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            auto r = linfit(img[id], 1.0, 1.0, psf[id]);
            flux = r.params[1];
            bg = r.params[0];
        }

        if (residual) {
            res = img - psf*flux;
        }
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
    bool central = false;
    bool residual = false;
    bool large_bg = false;
    vec1d raper;
    double raper_bg = dnan;
    double apcor = 1.0;
    bool auto_apcor = false;
    bool norm = false;
    double frac = 0.4;
    double ffrac = 0.1;
    double bfrac = dnan;

    vec2d psf_orig;
    vec2d psf;
    vec2b ipsf, ifpsf, opsf, mpsf;
    uint_t imax;

    bool update_masks_() {
        if (raper.empty()) {
            ipsf = fabs(psf) >= frac*max(fabs(psf));
            vec1u id = where(ipsf);
            if (id.empty()) {
                error("no pixel to fit (frac=", frac,")");
                return false;
            }

            opsf = fabs(psf) < bfrac*max(fabs(psf));
            id = where(opsf);
            if (id.empty()) {
                error("no pixel in background (bfrac=", bfrac,")");
                return false;
            }

            mpsf = (fabs(psf) > 0 && fabs(psf) < bfrac*max(fabs(psf))) || ipsf;
            id = where(mpsf);
            if (id.empty()) {
                error("no pixel to fit (bfrac=", bfrac,", frac=", frac, ")");
                return false;
            }

            imax = max_id(psf);
        } else {
            int_t x0 = psf.dims[0]/2, y0 = psf.dims[1]/2;
            ipsf = generate_img({psf.dims[0], psf.dims[1]}, [&](int_t x, int_t y) {
                return sqr(x-x0) + sqr(y-y0) < sqr(raper[0]);
            });
            if (raper.size() == 3) {
                opsf = generate_img({psf.dims[0], psf.dims[1]}, [&](int_t x, int_t y) {
                    double r = sqr(x-x0) + sqr(y-y0);
                    return r > sqr(raper[1]) && r < sqr(raper[2]);
                });
            } else if (raper.size() == 2) {
                opsf = generate_img({psf.dims[0], psf.dims[1]}, [&](int_t x, int_t y) {
                    return sqr(x-x0) + sqr(y-y0) > sqr(raper[1]);
                });
            } else {
                opsf = ipsf;
                opsf[_] = true;
            }

            if (auto_apcor) {
                apcor = 1.0/(total(sqr(psf[where(ipsf)]) - mean(sqr(psf[where(opsf)]))));
            }
        }

        ifpsf = fabs(psf) >= ffrac*max(fabs(psf));
        vec1u id = where(ifpsf);
        if (id.empty()) {
            error("no pixel to fit (ffrac=", ffrac,")");
            return false;
        }

        return true;
    }

    bool config(program_arguments& pa, bool init = false) {
        std::string tpsf;

        pa.read(arg_list(
            name(tpsf, "psf"), frac, ffrac, norm, nbstrap, name(tseed, "seed"), central, large_bg,
            residual, bfrac, name(raper, "aperture"), apcor
        ));

        if (apcor == 0.0) auto_apcor = true;

        if (!finite(bfrac)) bfrac = ffrac;

        if (!init && psf.empty() && tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        if (!tpsf.empty()) {
            fits::read(tpsf, psf);

            if (norm) {
                psf /= max(fabs(psf));
            }

            psf_orig = psf;

            if (!update_masks_()) {
                return false;
            }
        }

        return true;
    }

    vec2d med, dsp, res;

    void extract(const vec3d& cube, double& disp, double& bg) {
        if (psf.dims[0] != cube.dims[1] || psf.dims[1] != cube.dims[2]) {
            int_t hxsize = cube.dims[1]/2;
            int_t hysize = cube.dims[2]/2;
            vec1i mid = psf_orig.ids(max_id(psf_orig));
            psf = subregion(psf_orig, {mid[0]-hxsize, mid[1]-hysize, mid[0]+hxsize, mid[1]+hysize});

            if (!update_masks_()) {
                disp = dnan;
                bg = dnan;
                return;
            }
        }

        med.resize(cube.dims[1], cube.dims[2]);
        dsp.resize(cube.dims[1], cube.dims[2]);

        run_dim_idx(cube, 0, [&](uint_t i, vec1d& d) {
            med[i] = inplace_median(d);
            for (uint_t j : range(d)) {
                d[j] = fabs(d[j] - med[i]);
            }
            dsp[i] = inplace_median(d);
        });

        if (!raper.empty()) {
            vec1u idbg = where(opsf && finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel in aperture background");
                return;
            }

            bg = mean(sqr(dsp[idbg]));

            vec1u idd = where(ipsf && finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel in aperture");
                return;
            }

            disp = sqrt(apcor*total(sqr(dsp[idd]) - bg));
            bg = sqrt(bg);
        } else if (central) {
            vec1u idbg = where(opsf && finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            bg = median(dsp[idbg]);
            disp = sqrt(sqr(dsp[imax]) - sqr(bg))/psf[imax];
        } else if (large_bg) {
            vec1u idd = where(mpsf && finite(dsp));
            if (idd.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            auto r = linfit(sqr(dsp[idd]), 1.0, 1.0, sqr(psf[idd])*ipsf[idd]);
            disp = sqrt(r.params[1]);
            bg = 1.483*sqrt(r.params[0]);
        } else {
            vec1u idd = where(ipsf && finite(dsp));
            if (idd.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            auto r = linfit(sqr(dsp[idd]), 1.0, 1.0, sqr(psf[idd]));
            disp = sqrt(r.params[1]);
            bg = 1.483*sqrt(r.params[0]);
        }

        if (residual) {
            res = sqrt(sqr(dsp) - sqr(psf*disp));
        }

        double fmed;
        vec1u idf = where(ifpsf && finite(med));
        if (idf.empty()) {
            disp = dnan;
            bg = dnan;
            error("no pixel to fit");
            return;
        }

        auto r = linfit(med[idf], 1.0, 1.0, psf[idf]);
        fmed = r.params[1];

        disp /= fmed;
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
    bool residual = false;
    std::string out = "";

    flux_extractor ex;

    {
        program_arguments pa(argc, argv);
        pa.read(arg_list(bstrap, residual, out));
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
        if (residual) {
            fits::write(out+"_res.fits", ex.res);
        }
    }

    return true;
}

bool get_logdisp(int argc, char* argv[], const vec3d& cube) {
    bool bstrap = false;
    bool residual = false;
    std::string out = "";

    logdisp_extractor ex;

    {
        program_arguments pa(argc, argv);
        pa.read(arg_list(bstrap, residual, out));
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
        if (residual) {
            fits::write(out+"_res.fits", ex.res);
        }
    }

    return true;
}

bool run_batch(int argc, char* argv[], const std::string& file) {
    std::string op = argv[1];

    if (op == "flux") {
        bool bstrap = false;
        std::string out = "";

        flux_extractor ex;

        {
            program_arguments pa(argc-1, argv+1);
            pa.read(arg_list(bstrap, out));
            if (!ex.config(pa, true)) return false;
        }

        vec1s cfile;
        vec1d flux, bg, flux_err, bg_err;

        vec1s lines; {
            std::ifstream bfile(file);
            while (!bfile.eof()) {
                std::string line;
                std::getline(bfile, line);
                line = trim(line);

                if (line.empty()) continue;
                lines.push_back(line);
            }
        }

        auto pg = progress_start(lines.size());
        for (auto& line : lines) {
            vec1s args = trim(split(line, " "));
            cfile.push_back(args[0]);
            if (args.size() > 1) {
                args = args[uindgen(args.size()-1)+1];
                program_arguments tmp(args);
                ex.config(tmp);
            }

            vec3d cube; fits::read(cfile.data.back(), cube);
            flux.push_back(0.0); bg.push_back(0.0);
            ex.extract(cube, flux.data.back(), bg.data.back());

            if (bstrap) {
                flux_err.push_back(0.0); bg_err.push_back(0.0);
                ex.bootstrap(cube, flux_err.data.back(), bg_err.data.back());
            }

            progress(pg);
        }

        if (out.empty()) {
            out = "fluxcube.fits";
            warning("output file not set (out=...), using 'fluxcube.fits'");
        } else {
            file::mkdir(file::get_directory(out));
        }

        fits::write_table(out, ftable(cfile, flux, bg, flux_err, bg_err));

        return true;
    } else if (op == "logdisp") {
        bool bstrap = false;
        std::string out = "";

        logdisp_extractor ex;

        {
            program_arguments pa(argc-1, argv+1);
            pa.read(arg_list(bstrap, out));
            if (!ex.config(pa, true)) return false;
        }

        vec1s cfile;
        vec1d disp, bg, disp_err, bg_err;
        vec1d apcor;

        vec1s lines; {
            std::ifstream bfile(file);
            while (!bfile.eof()) {
                std::string line;
                std::getline(bfile, line);
                line = trim(line);

                if (line.empty()) continue;
                lines.push_back(line);
            }
        }

        auto pg = progress_start(lines.size());
        for (auto& line : lines) {
            vec1s args = trim(split(line, " "));
            cfile.push_back(args[0]);
            if (args.size() > 1) {
                args = args[uindgen(args.size()-1)+1];
                program_arguments tmp(args);
                ex.config(tmp);
            }

            vec3d cube; fits::read(cfile.data.back(), cube);
            disp.push_back(0.0); bg.push_back(0.0);
            ex.extract(cube, disp.data.back(), bg.data.back());
            if (!ex.raper.empty()) {
                apcor.push_back(ex.apcor);
            }

            if (bstrap) {
                disp_err.push_back(0.0); bg_err.push_back(0.0);
                ex.bootstrap(cube, disp_err.data.back(), bg_err.data.back());
            }

            progress(pg);
        }

        if (out.empty()) {
            out = "fluxcube.fits";
            warning("output file not set (out=...), using 'fluxcube.fits'");
        } else {
            file::mkdir(file::get_directory(out));
        }

        fits::write_table(out, ftable(cfile, disp, bg, disp_err, bg_err, apcor));

        return true;
    } else {
        error("unknown operation '"+op+"'");
        return false;
    }
}

void print_help() {
    using namespace format;

    print("fluxcube v1.0");
    header("Usage: fluxcube cube.fits [options]");
    header("Available operations:");
    bullet("flux", "perform mean or median stacking then PSF fitting to get the flux");
    bullet("logdisp", "perform MAD stacking then PSF fitting to get the lognormal dispersion");
    bullet("batch", "run an operation on a list of cubes provided in a batch file");

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
