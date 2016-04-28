#include <phypp.hpp>
#include <phypp/astro/qstack.hpp>

void print_help() {
    using namespace format;

    print("qstack v2.0");
    paragraph("usage: qstack2 cat=\"...\" img=\"...\" out=\"...\" hsize=... [options=...]");

    paragraph(
        "The program will look for RA and Dec coordinates in 'cat' (FITS catalog), load WCS data "
        "from 'img' (FITS image), then stack all the sources of the catalog that are present in "
        "the image into a 2*hsize+1 square cutout, eventually saving the result into 'out' (FITS "
        "image)."
    );

    paragraph(
        "If 'cat' is empty, then 'img' is assumed to be a cube of cutouts that will be collapsed "
        "by the program to form a single stacked cutout."
    );

    header("List of available command line options:");
    bullet("cat", "[string] path to the source catalog (FITS file)");
    bullet("img", "[string] path to the flux map (FITS file)");
    bullet("out", "[string] path to the output file");
    bullet("hsize", "[unsigned integer] half size of the resulting cutout");
    bullet("wht", "[string, optional] path to the weight map (FITS file), giving the weight of "
        "each pixel in 'img'");
    bullet("err", "[string, optional] path to the error map (or RMS map) (FITS file), giving the "
        "estimated 1-sigma error on the measured flux in 'img' (will be converted to a weight by "
        "inverse square root)");
    bullet("pos", "[string, optional] prefix of the RA and Dec coordinates in 'cat'");
    bullet("ids", "[2 x integer, optional] inclusive range of IDs to stack in the catalog (or cube)");
    bullet("   ", "note: use '*' to symbolize the beginning or the end");
    bullet("subpixel", "[flag] perform sub-pixel interpolation (default: false)");
    bullet("mean", "[flag] perform mean stacking (default)");
    bullet("median", "[flag] perform median stacking");
    bullet("cube", "[flag] do not stack, just ouput the cube");
    bullet("keepnan", "[flag] do not reject sources with NaN pixels");
    bullet("bstrap", "[flag] perform bootstraping and save the resulting cube");
    bullet("nbstrap", "[unsigned integer, optional] number of boostraping realisations");
    bullet("sbstrap", "[unsigned integer, optional] size of a boostraping realisation");
    bullet("randomize", "[flag] when stacking from a catalog, randomize the positions inside the "
        "area that is covered by the selected sources");
    bullet("verbose", "[flag] print some information about the stacking process");
    print("");

    paragraph("Copyright (c) 2013 C. Schreiber (corentin.schreiber@cea.fr)");

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

bool read_ids(std::string sid, uint_t& id, uint_t pos, uint_t star, uint_t idmax) {
    if (sid == "*") {
        id = star;
    } else {
        if (!from_string(sid, id)) {
            error("could not read 'ids[", pos, "]': '", sid, "'");
            return false;
        }

        if (id >= idmax) {
            error("'ids[", pos, "]' is too large, only ", idmax, " sources available");
            return false;
        }
    }

    return true;
}

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string cat = "";
    vec1s img;
    vec1s wht;
    vec1s err;
    std::string pos = "";
    std::string out = "";
    uint_t hsize = npos;
    bool median = false;
    bool mean = false;
    bool verbose = false;
    bool bstrap = false;
    bool subpixel = false;
    uint_t randomize = 0;
    bool tcube = false;
    bool keepnan = false;
    uint_t nbstrap = 200;
    uint_t sbstrap = 0;
    uint_t tseed = 42;
    vec1s cids;

    read_args(argc, argv, arg_list(
        out, cat, img, wht, err, pos, hsize, median, mean, bstrap, nbstrap, sbstrap,
        randomize, name(tseed, "seed"), name(tcube, "cube"), subpixel, verbose, keepnan,
        name(cids, "ids")
    ));

    auto seed = make_seed(tseed);

    if (!median && !mean) {
        mean = true;
    } else if (mean && median) {
        median = false;
    }

    if (!cat.empty() && hsize == npos) {
        error("missing cutout size 'hsize'\n");
        print_help();
        return 1;
    }

    if (out.empty()) {
        error("missing output file name 'out'\n");
        print_help();
        return 1;
    }

    if (cat.empty()) {
        if (img.empty()) {
            error("missing catalog file 'cat'\n");
            print_help();
            return 1;
        }
    } else if (!file::exists(cat)) {
        error("cannot find '"+cat+"'\n");
        return 1;
    }

    if (img.empty()) {
        error("missing image file 'img'\n");
        print_help();
        return 1;
    } else {
        vec1u i = where(!file::exists(img));
        if (!i.empty()) {
            error("cannot find '"+img[i]+"'\n");
            return 1;
        }
    }

    if (!wht.empty()) {
        vec1u i = where(!file::exists(wht));
        if (!i.empty()) {
            error("cannot find '"+wht[i]+"'\n");
            return 1;
        }
    }

    if (!err.empty()) {
        vec1u i = where(!file::exists(err));
        if (!i.empty()) {
            error("cannot find '"+err[i]+"'\n");
            return 1;
        }
    }

    if (!cids.empty() && cids.size() != 2) {
        error("'ids' option needs exactly 2 elements (min and maximum ID)");
        return 1;
    }

    struct {
        vec1d ra, dec;
    } fcat;

    if (!cat.empty()) {
        std::string posh = pos.empty() ? "" : pos+".";
        fits::read_table(cat, toupper(posh+"ra"), fcat.ra, toupper(posh+"dec"), fcat.dec);

        if (!cids.empty()) {
            vec1u tcids(2);
            if (!read_ids(cids[0], tcids[0], 0, 0, fcat.ra.size())) return 1;
            if (!read_ids(cids[1], tcids[1], 1, fcat.ra.size()-1, fcat.ra.size())) return 1;

            vec1u tid = rgen(tcids[0], tcids[1]);
            fcat.ra = fcat.ra[tid];
            fcat.dec = fcat.dec[tid];
        }

        vec1u id = where(is_finite(fcat.ra) && is_finite(fcat.dec));
        fcat.ra = fcat.ra[id];
        fcat.dec = fcat.dec[id];

        if (randomize > 0) {
            if (randomize == 1) randomize = fcat.ra.size();

            vec1d rra = {min(fcat.ra), max(fcat.ra)};
            vec1d rdec = {min(fcat.dec), max(fcat.dec)};
            auto hull = build_convex_hull(fcat.ra, fcat.dec);

            randpos_uniform_box(seed, randomize, rra, rdec, fcat.ra, fcat.dec,
                [&](double tra, double tdec) {
                    return in_convex_hull(tra, tdec, hull);
                }
            );
        }
    }

    vec2f stack;
    vec3f cube;
    vec3f bs;

    auto do_subpixel = [&cube](const vec1d& dx, const vec1d& dy) {
        for (uint_t i : range(dx)) {
            cube(i,_,_) = translate(cube(i,_,_), dy[i], dx[i]);
        }
    };

    qstack_params params;
    params.keep_nan = keepnan;
    params.verbose = verbose;
    params.save_offsets = subpixel;

    if ((wht.empty() && err.empty()) || median) {
        if (cat.empty()) {
            fits::read(img[0], cube);

            for (uint_t i : range(1, img.size())) {
                vec3f tc;
                fits::read(img[i], tc);
                append<0>(cube, tc);
            }

            if (!cids.empty()) {
                vec1u tcids(2);
                if (!read_ids(cids[0], tcids[0], 0, 0, cube.dims[0])) return 1;
                if (!read_ids(cids[1], tcids[1], 1, cube.dims[0]-1, cube.dims[0])) return 1;

                vec1u tid = rgen(tcids[0], tcids[1]);
                cube = cube(tid, _, _);
            }

            if (verbose) print("stacking ", cube.dims[0], " sources");
        } else {
            if (img.size() != 1) {
                error("catalog stacking needs a single image");
                return 1;
            }

            vec1u ids;
            qstack_output qout = qstack(fcat.ra, fcat.dec, img[0], hsize, cube, ids, params);
            if (verbose) print("stacking ", ids.size(), "/", fcat.ra.size(), " sources");

            if (subpixel) {
                do_subpixel(qout.dx, qout.dy);
            }
        }

        if (!tcube) {
            if (mean) {
                stack = qstack_mean(cube);
            } else {
                stack = qstack_median(cube);
            }

            if (bstrap) {
                if (sbstrap == 0) sbstrap = cube.dims[0]/2;
                if (mean) {
                    bs = qstack_mean_bootstrap(cube, nbstrap, sbstrap, seed);
                } else {
                    bs = qstack_median_bootstrap(cube, nbstrap, sbstrap, seed);
                }
            }
        }
    } else {
        if (median) {
            error("cannot use weights when doing median stacking");
            return 1;
        }

        vec3f wcube;
        if (cat.empty()) {
            fits::read(img[0], cube);
            fits::read(wht.empty() ? err[0] : wht[0], wcube);

            for (uint_t i : range(1, img.size())) {
                vec3f tc;
                fits::read(img[i], tc);
                append<0>(cube, tc);
                fits::read(wht.empty() ? err[i] : wht[i], tc);
                append<0>(wcube, tc);
            }

            if (!cids.empty()) {
                vec1u tcids(2);
                if (!read_ids(cids[0], tcids[0], 0, 0, cube.dims[0])) return 1;
                if (!read_ids(cids[1], tcids[1], 1, cube.dims[0]-1, cube.dims[0])) return 1;

                vec1u tid = rgen(tcids[0], tcids[1]);
                cube = cube(tid, _, _);
                wcube = wcube(tid, _, _);
            }

            if (verbose) print("stacking ", cube.dims[0], " sources");
        } else {
            if (img.size() != 1) {
                error("catalog stacking needs a single image");
                return 1;
            }

            vec1u ids;
            qstack_output qout = qstack(fcat.ra, fcat.dec, img[0],
                wht.empty() ? err[0] : wht[0], hsize, cube, wcube, ids, params);
            if (verbose) print("stacking ", ids.size(), "/", fcat.ra.size(), " sources");

            if (subpixel) {
                do_subpixel(qout.dx, qout.dy);
            }
        }

        if (!tcube) {
            if (!wht.empty()) {
                stack = qstack_mean(cube, wcube);
            } else {
                stack = qstack_mean(cube, invsqr(wcube));
            }

            if (bstrap) {
                if (sbstrap == 0) sbstrap = cube.dims[0]/2;
                bs = qstack_mean_bootstrap(cube, nbstrap, sbstrap, seed);
            }
        }
    }

    out = trim(out);
    file::mkdir(file::get_directory(out));

    if (tcube) {
        fits::write(out, cube);
    } else {
        fits::write(out, stack);
    }

    if (bstrap) {
        std::string bstrap_out;
        if (end_with(out, ".fits")) {
            bstrap_out = out.substr(0, out.size()-5) + "_bs.fits";
        } else {
            bstrap_out = file::get_directory(out) + "bstrap.fits";
        }

        fits::write(bstrap_out, bs);
    }

    return 0;
}
