#include <phypp.hpp>

void header(const std::string& msg) {
    vec1s w = wrap("  "+msg, 80, "  ");
    for (auto& s : w) {
        print(s);
    }
}

void paragraph(const std::string& msg) {
    vec1s w = wrap("  "+msg, 80, "  ");
    for (auto& s : w) {
        print(s);
    }
    print("");
}

void bullet(const std::string& name, const std::string& desc) {
    std::string header = "    "+name+": ";
    vec1s w = wrap(header+desc, 80, std::string(header.size(), ' '));
    for (auto& s : w) {
        print(s);
    }
}

void print_help() {
    print("qstack v2.0");
    paragraph("usage: qstack2 cat=\"...\" img=\"...\" out=\"...\" hsize=... [options=...]");

    paragraph(
        "The program will look for RA and Dec coordinates in 'cat' (FITS catalog), load WCS data "
        "from 'img' (FITS image), then stack all the sources of the catalog that are present in "
        "the image into a 2*hsize+1 square cutout, eventually saving the result into 'out' (FITS "
        "image)."
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
    bullet("mean", "[flag] perform mean stacking (default)");
    bullet("median", "[flag] perform median stacking");
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string cat = "";
    std::string img = "";
    std::string wht = "";
    std::string err = "";
    std::string pos = "";
    std::string out = "";
    uint_t hsize = npos;
    bool median = false;
    bool mean = false;
    bool verbose = false;
    read_args(argc, argv, arg_list(out, cat, img, wht, err, pos, hsize, median, mean, verbose));

    if (!median && !mean) {
        mean = true;
    } else if (mean && median) {
        median = false;
    }

    if (hsize == npos) {
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
        error("missing catalog file 'cat'\n");
        print_help();
        return 1;
    } else if (!file::exists(cat)) {
        error("cannot find '"+cat+"'\n");
        return 1;
    }

    if (img.empty()) {
        error("missing image file 'img'\n");
        print_help();
        return 1;
    } else if (!file::exists(img)) {
        error("cannot find '"+img+"'\n");
        return 1;
    }

    if (!wht.empty() && !file::exists(wht)) {
        error("cannot find '"+wht+"'\n");
        return 1;
    }

    if (!err.empty() && !file::exists(err)) {
        error("cannot find '"+err+"'\n");
        return 1;
    }

    struct {
        vec1d ra, dec;
    } fcat;

    std::string posh = pos.empty() ? "" : pos+".";

    fits::read_table(cat, toupper(posh+"ra"), fcat.ra, toupper(posh+"dec"), fcat.dec);

    vec3f cube(0, (2*hsize+1), (2*hsize+1));
    vec1u fids;
    qstack(fcat.ra, fcat.dec, img, hsize, cube, fids);

    vec2f stack;
    if ((wht.empty() && err.empty()) || median) {
        if (verbose) print("stacking ", cube.dims[0], "/", fcat.ra.size(), " sources");
        if (mean) {
            stack = qstack_mean(cube);
        } else {
            stack = qstack_median(cube);
        }
    } else {
        if (median) {
            warning("cannot use weights when doing median stacking");
        }

        if (!wht.empty()) {
            vec3f wcube;
            vec1u wids;
            qstack(fcat.ra, fcat.dec, wht, hsize, wcube, wids);

            vec1u ids1, ids2;
            match(fids, wids, ids1, ids2);

            if (verbose) print("stacking ", ids1.size(), "/", fcat.ra.size(), " sources");

            stack = qstack_mean(cube(ids1,_,_), wcube(ids2,_,_));
        } else {
            vec3f ecube;
            vec1u eids;
            qstack(fcat.ra, fcat.dec, err, hsize, ecube, eids);

            vec1u ids1, ids2;
            match(fids, eids, ids1, ids2);

            if (verbose) print("stacking ", ids1.size(), "/", fcat.ra.size(), " sources");

            stack = qstack_mean(cube(ids1,_,_), pow(ecube(ids2,_,_), -2));
        }
    }

    file::mkdir(file::get_directory(out));
    fits::write(out, stack);

    return 0;
}
