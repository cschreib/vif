#include <vif.hpp>
#include <vif/astro/qxmatch.hpp>

void print_help() {
    using namespace terminal_format;

    print("qxmatch2 v2.0");
    paragraph("usage: qxmatch2 cats=[file1,file2] output=ofile [options=...]");

    paragraph(
        "'file1' and 'file2' must be FITS files, expected to contain an extension with "
        "at least two vector columns named 'RA' and 'DEC', each containing a single row. "
        "Both are assumed to contain galactic coordinates in degrees. "
        "They can be created in IDL using: "
    );
    print("        mwrfits, {ra:[...], dec:[...]}, \"filename.fits\", /create\n");
    paragraph(
        "The output file is written in a similar way: it contains one extension named "
        "'RESULT_BINARY', which contains a single column named 'ID', itself being a two "
        "dimensional vector containing the ID of the n'th closest sources. It can be  "
        "read in IDL using: "
    );
    print("        res = mrdfits(\"output.fits\", 1)");
    print("        best = res.id[*,0]");
    print("        second_best = res.id[*,1]\n");
    paragraph(
        "The program also outputs the crossmatch distance in a column called 'D', with a "
        "format similar to that of 'ID', in arcseconds."
    );

    paragraph(
        "It is also possible to compute an 'auto crossmatch' of a single catalog by only "
        "providing one file in the 'cats' parameter. The program will then compute the "
        "n'th nearest neighbors for each source within this catalog."
    );

    header("List of available command line options:");
    bullet("verbose", "set this flag to print additional information in the standard output");
    bullet("nth", "[number] set this value to the number of closest neighbors you want to retrieve "
        "(default: 1).");
    bullet("pos", "[string(array)]: defines the suffix of the RA and Dec variables inside the two "
        "catalogs (default: \"\")");
    bullet("radec1", "[string(array)]: defines the name of the RA and Dec variables in the first "
        "catalog (default: [ra,dec])");
    bullet("radec2", "[string(array)]: defines the name of the RA and Dec variables in the first "
        "catalog (default: [ra,dec])");
    bullet("thread", "[number]: set this value to the number of concurrent threads you want to run "
        "(default: 1).");
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

int vif_main(int argc, char* argv[]) {
    vec1s cats;
    std::string output;
    vec1s pos;

    vec1s radec1 = {"ra","dec"};
    vec1s radec2 = {"ra","dec"};
    uint_t nth = 1;
    uint_t thread = 1;
    bool   verbose = false;
    bool   quiet = false;
    bool   brute = false;
    bool   no_mirror = false;

    read_args(argc, argv, arg_list(
        cats, output, nth, thread, verbose, quiet, pos, radec1, radec2, brute, no_mirror
    ));

    if (quiet) verbose = false;

    if (output.empty()) {
        if (!quiet) print_help();
        return 0;
    }

    struct coord_t {
        vec1d ra, dec;
    };

    qxmatch_res res;

    if (cats.size() == 2) {
        if (pos.size() == 1) pos = replicate(pos[0], 2);
        else if (pos.empty()) pos = {"", ""};
        vec1u idne = where(!empty(pos));
        pos[idne] += ".";

        coord_t cat1, cat2;
        fits::read_table(cats[0], pos[0]+radec1[0], cat1.ra, pos[0]+radec1[1], cat1.dec);
        fits::read_table(cats[1], pos[1]+radec2[0], cat2.ra, pos[1]+radec2[1], cat2.dec);

        if (cat1.ra.size() == 0 || cat1.dec.size() == 0) {
            error("qxmatch: first catalog is empty");
        }

        if (cat2.ra.size() == 0 || cat2.dec.size() == 0) {
            error("qxmatch: second catalog is empty");
        }

        if (verbose) {
            print("qxmatch: crossmatching ", n_elements(cat1.ra), " sources from '", cats[0],
                "' with ", n_elements(cat2.ra), " sources from '", cats[1], "'...");
        }

        qxmatch_params p; p.nth = nth; p.thread = thread; p.verbose = verbose; p.no_mirror = no_mirror;
        p.brute_force = brute;
        res = qxmatch(cat1, cat2, p);
    } else if (cats.size() == 1) {
        if (pos.empty()) pos = {""};
        if (!pos[0].empty()) pos += ".";

        coord_t cat;
        fits::read_table(cats[0], pos[0]+radec1[0], cat.ra, pos[0]+radec1[1], cat.dec);

        if (cat.ra.size() == 0 || cat.dec.size() == 0) {
            error("qxmatch: catalog is empty");
        }

        if (verbose) {
            print("qxmatch: self matching ", n_elements(cat.ra), " sources from '", cats[0], "'...");
        }

        qxmatch_params p; p.nth = nth; p.thread = thread; p.verbose = verbose; p.no_mirror = no_mirror;
        p.brute_force = brute;
        res = qxmatch(cat, p);
    } else {
        if (!quiet) print_help();
        return 0;
    }

    qxmatch_save(output, res);

    return 0;
}
