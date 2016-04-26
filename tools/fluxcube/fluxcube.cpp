#include <phypp.hpp>
#include <phypp/astro/qstack.hpp>

void print_help();

bool get_flux(int argc, char* argv[], const std::string& file);
bool get_logdisp(int argc, char* argv[], const vec3d& cube);
bool run_batch(int argc, char* argv[], const std::string& file);

int phypp_main(int argc, char* argv[]) {
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
        return get_flux(argc-2, argv+2, file) ? 0 : 1;
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

template<std::size_t D, typename T>
rtype_t<T> fstddev(const vec<D,T>& v) {
    return stddev(v[where(is_finite(v))]);
}

template<std::size_t D, typename T>
rtype_t<T> frms(const vec<D,T>& v) {
    return rms(v[where(is_finite(v))]);
}

template<std::size_t D, typename T>
rtype_t<T> fmean(const vec<D,T>& v) {
    return mean(v[where(is_finite(v))]);
}

template<std::size_t D, typename T>
rtype_t<T> fmedian(const vec<D,T>& v) {
    return median(v[where(is_finite(v))]);
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
    double psf_err_factor;

    bool update_masks_() {
        ipsf = abs(psf) >= frac*max(abs(psf));
        vec1u id = where(ipsf);
        if (id.empty()) {
            psf_err_factor = dnan;
            error("no pixel to fit (frac=", frac,")");
            return false;
        }

        psf_err_factor = 1.0/sqrt(total(sqr(psf[id])) - sqr(total(psf[id]))/id.size());

        return true;
    }

    bool config(program_arguments& pa, bool init = false) {
        std::string tpsf;
        bool beam_smoothed = false;
        float smooth_radius = fnan;
        bool med = false;
        bool norm = false;

        pa.read(arg_list(
            name(tpsf, "psf"), name(mea, "mean"), name(med, "median"), frac, norm, nbstrap,
            name(tseed, "seed"), large_bg, residual, beam_smoothed, smooth_radius
        ));

        if (med) mea = false;

        if (!init && psf.empty() && tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        if (!tpsf.empty()) {
            fits::read(tpsf, psf);

            if (norm) {
                psf /= max(abs(psf));
            }

            if (beam_smoothed) {
                vec2d kernel;
                if (is_finite(smooth_radius)) {
                    kernel = gaussian_profile(psf.dims, smooth_radius);
                } else {
                    kernel = psf;
                }
                vec1u idz = where(kernel < 1e-5*max(abs(kernel)));
                kernel[idz] = 0.0;
                psf = convolve2d(psf, kernel);
                psf /= max(abs(psf));
            }

            psf_orig = psf;

            if (!update_masks_()) {
                return false;
            }
        }

        return true;
    }

    vec2d img, res;

    const vec1s names = {"flux", "bg", "res_rms", "res_err"};
    const uint_t nresult = names.size();

    void extract(const vec2d& timg, vec1d& result) {
        if (large_bg) {
            vec1u id = where(is_finite(timg));
            if (id.empty()) {
                result = replicate(dnan, nresult);
                error("no pixel to fit");
                return;
            }

            auto r = linfit(timg[id], 1.0, 1.0, psf[id]*ipsf[id]);
            result = replicate(dnan, nresult);
            result[0] = r.params[1];
            result[1] = r.params[0];
        } else {
            vec1u id = where(ipsf && is_finite(timg));
            if (id.empty()) {
                result = replicate(dnan, nresult);
                error("no pixel to fit");
                return;
            }

            auto r = linfit(timg[id], 1.0, 1.0, psf[id]);
            result = replicate(dnan, nresult);
            result[0] = r.params[1];
            result[1] = r.params[0];
        }

        if (residual) {
            res = timg - psf*result[0];
            result[2] = stddev(res);
            result[3] = result[2]*psf_err_factor;
        }
    }

    void extract(const vec3d& cube, vec1d& result) {
        if (psf.dims[0] != cube.dims[1] || psf.dims[1] != cube.dims[2]) {
            int_t hxsize = cube.dims[1]/2;
            int_t hysize = cube.dims[2]/2;
            vec1i mid = mult_ids(psf_orig, max_id(psf_orig));
            psf = subregion(psf_orig, {mid[0]-hxsize, mid[1]-hysize, mid[0]+hxsize, mid[1]+hysize});

            if (!update_masks_()) {
                result = replicate(dnan, nresult);
                return;
            }
        }

        if (mea) {
            img = qstack_mean(cube);
        } else {
            img = qstack_median(cube);
        }

        extract(img, result);
    }

    void bootstrap(const vec3d& cube, vec1d& result, vec1d& err) {
        auto seed = make_seed(tseed);

        vec2d rs; rs.reserve(nbstrap*nresult);
        uint_t nsrc = cube.dims[0];
        qstack_bootstrap(cube, nbstrap, nsrc/2, seed, [&](const vec3d& sub) {
            vec1d tr;
            extract(sub, tr);
            rs.push_back(tr);
        });

        result = replicate(dnan, nresult);
        err = replicate(dnan, nresult);

        uint_t tnr = residual ? nresult : nresult - 2;

        for (uint_t i : range(tnr)) {
            if (mea) {
                result[i] = fmean(rs(_,i));
            } else {
                result[i] = fmedian(rs(_,i));
            }
            rs(_,i) -= result[i];
            err[i] = frms(rs(_,i))/sqrt(2.0);
        }
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
            ipsf = abs(psf) >= frac*max(abs(psf));
            vec1u id = where(ipsf);
            if (id.empty()) {
                error("no pixel to fit (frac=", frac,")");
                return false;
            }

            opsf = abs(psf) < bfrac*max(abs(psf));
            id = where(opsf);
            if (id.empty()) {
                error("no pixel in background (bfrac=", bfrac,")");
                return false;
            }

            mpsf = (abs(psf) > 0 && abs(psf) < bfrac*max(abs(psf))) || ipsf;
            id = where(mpsf);
            if (id.empty()) {
                error("no pixel to fit (bfrac=", bfrac,", frac=", frac, ")");
                return false;
            }

            imax = max_id(psf);
        } else {
            int_t x0 = psf.dims[0]/2, y0 = psf.dims[1]/2;
            ipsf = generate_img(psf.dims, [&](int_t x, int_t y) {
                return sqr(x-x0) + sqr(y-y0) < sqr(raper[0]);
            });
            if (raper.size() == 3) {
                opsf = generate_img(psf.dims, [&](int_t x, int_t y) {
                    double r = sqr(x-x0) + sqr(y-y0);
                    return r > sqr(raper[1]) && r < sqr(raper[2]);
                });
            } else if (raper.size() == 2) {
                opsf = generate_img(psf.dims, [&](int_t x, int_t y) {
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

        ifpsf = abs(psf) >= ffrac*max(abs(psf));
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

        if (!is_finite(bfrac)) bfrac = ffrac;

        if (!init && psf.empty() && tpsf.empty()) {
            error("missing PSF file: psf=...");
            return false;
        }

        if (!tpsf.empty()) {
            fits::read(tpsf, psf);

            if (norm) {
                psf /= max(abs(psf));
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
            vec1i mid = mult_ids(psf_orig, max_id(psf_orig));
            psf = subregion(psf_orig, {mid[0]-hxsize, mid[1]-hysize, mid[0]+hxsize, mid[1]+hysize});

            if (!update_masks_()) {
                disp = dnan;
                bg = dnan;
                return;
            }
        }

        med.resize(cube.dims[1], cube.dims[2]);
        dsp.resize(cube.dims[1], cube.dims[2]);

        run_dim(0, cube, [&](uint_t i, vec1d d) {
            med[i] = inplace_median(d);
            for (uint_t j : range(d)) {
                d[j] = abs(d[j] - med[i]);
            }
            dsp[i] = inplace_median(d);
        });

        if (!raper.empty()) {
            vec1u idbg = where(opsf && is_finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel in aperture background");
                return;
            }

            bg = mean(sqr(dsp[idbg]));

            vec1u idd = where(ipsf && is_finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel in aperture");
                return;
            }

            disp = sqrt(apcor*total(sqr(dsp[idd]) - bg));
            bg = sqrt(bg);
        } else if (central) {
            vec1u idbg = where(opsf && is_finite(dsp));
            if (idbg.empty()) {
                disp = dnan;
                bg = dnan;
                error("no pixel to fit");
                return;
            }

            bg = median(dsp[idbg]);
            disp = sqrt(sqr(dsp[imax]) - sqr(bg))/psf[imax];
        } else if (large_bg) {
            vec1u idd = where(mpsf && is_finite(dsp));
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
            vec1u idd = where(ipsf && is_finite(dsp));
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
        vec1u idf = where(ifpsf && is_finite(med));
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

        disps = disps[where(is_finite(disps))];
        bgs = bgs[where(is_finite(bgs))];

        disp = stddev(disps)/sqrt(2);
        bg = stddev(bgs)/sqrt(2);
    }
};

bool get_flux(int argc, char* argv[], const std::string& file) {
    bool bstrap = false;
    bool residual = false;
    std::string out = "";

    flux_extractor ex;

    {
        program_arguments pa(argc, argv);
        pa.read(arg_list(bstrap, residual, out));
        if (!ex.config(pa)) return false;
    }


    vec1d res;
    if (fits::is_cube(file)) {
        vec3d cube; fits::read(file, cube);
        print("nsrc: ", cube.dims[0]);

        ex.extract(cube, res);
        if (bstrap) {
            vec1d bs, bs_err;
            ex.bootstrap(cube, bs, bs_err);
            auto do_print = [&](vec1u ids) {
                for (uint_t i : ids) {
                    print(ex.names[i]+": ", res[i], " +/- ", bs_err[i], ", sys: ", bs[i]);
                }
            };

            do_print({0,1});
            if (residual) {
                do_print({2,3});
            }
        } else {
            auto do_print = [&](vec1u ids) {
                for (uint_t i : ids) {
                    print(ex.names[i]+": ", res[i]);
                }
            };

            do_print({0,1});
            if (residual) {
                do_print({2,3});
            }
        }
    } else {
        vec2d img; fits::read(file, img);
        ex.extract(img, res);
        print("single image");

        auto do_print = [&](vec1u ids) {
            for (uint_t i : ids) {
                print(ex.names[i]+": ", res[i]);
            }
        };

        do_print({0,1});
        if (residual) {
            do_print({2,3});
        }
    }

    if (!out.empty()) {
        file::mkdir(file::get_directory(out));
        if (fits::is_cube(file)) {
            fits::write(out+".fits", ex.img);
        }
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
        file::mkdir(file::get_directory(out));
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
        vec2d rs, bs, bs_err;

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

        rs.reserve(lines.size()*ex.nresult);
        if (bstrap) {
            bs.reserve(lines.size()*ex.nresult);
            bs_err.reserve(lines.size()*ex.nresult);
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
            vec1d tr;
            ex.extract(cube, tr);
            rs.push_back(tr);

            if (bstrap) {
                vec1d tb, tbe;
                ex.bootstrap(cube, tb, tbe);
                bs.push_back(tb);
                bs_err.push_back(tbe);
            }

            progress(pg);
        }

        if (out.empty()) {
            out = "fluxcube.fits";
            warning("output file not set (out=...), using 'fluxcube.fits'");
        } else {
            file::mkdir(file::get_directory(out));
        }

        if (bstrap) {
            fits::write_table(out, ftable(cfile, ex.names, name(rs, "result"), bs, bs_err));
        } else {
            fits::write_table(out, ftable(cfile, ex.names, name(rs, "result")));
        }

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
    header("Usage: fluxcube <cube>.fits <operation> [options]");
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
