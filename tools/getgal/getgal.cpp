#include <phypp.hpp>
#include <phypp/astro/qstack.hpp>

void print_help();

int main(int argc, char* argv[]) {
    vec1s tsrc;
    std::string out = "";
    std::string nbase = "";
    std::string dir = "";
    double radius = dnan;
    vec1u thsize;
    vec1s show;
    vec1s bands;
    bool show_rgb = false;
    bool verbose = false;

    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string clist = argv[1];

    read_args(argc-1, argv+1, arg_list(
        name(tsrc, "src"), out, name(nbase, "name"), dir, verbose, show, show_rgb, radius,
        name(thsize, "hsize"), bands
    ));

    if (!dir.empty()) {
        dir = file::directorize(dir);
    } else {
        dir = file::get_directory(clist);
    }

    if (!out.empty()) {
        out = file::directorize(out);
        if (!file::mkdir(out)) {
            warning("could not create directory '"+out+"'");
        }
    }

    vec1d ra, dec;
    vec1s name;
    if (tsrc.size() == 1 && end_with(tsrc[0], ".fits")) {
        struct {
            vec1d ra, dec;
            vec1s name;
            vec1u id;
        } tmp;

        fits::read_table_loose(tsrc[0], tmp);
        if (tmp.ra.empty() || tmp.dec.empty()) {
            error("missing RA and Dec coordinates in this FITS file");
            return 1;
        }

        if (tmp.name.empty()) {
            if (tmp.id.empty()) {
                tmp.id = uindgen(tmp.ra.size());
            }

            name = strna(tmp.id);
        } else {
            name = tmp.name;
        }

        name += "_";
        if (!nbase.empty()) name = nbase + "_" + name;

        ra = tmp.ra;
        dec = tmp.dec;
    } else if (tsrc.size() == 1 && end_with(tsrc[0], ".reg")) {
        file::read_table(tsrc[0], 0, ra, dec);
        name = strna(uindgen(ra.size())) + "_";
    } else if (tsrc.size() == 2) {
        name.resize(1);
        if (!nbase.empty()) name[0] = nbase + "_";
        ra.resize(1);
        dec.resize(1);
        if (find(tsrc[0], ":") == npos) {
            if (!from_string(tsrc[0], ra[0])) {
                error("could not convert '"+tsrc[0]+"' to a coordinate");
                return 1;
            }
            if (!from_string(tsrc[1], dec[0])) {
                error("could not convert '"+tsrc[1]+"' to a coordinate");
                return 1;
            }
        } else {
            sex2deg(tsrc[0], tsrc[1], ra[0], dec[0]);
        }
    } else {
        error("incorrect format for 'src' parameter");
        print_help();
        return 1;
    }

    if (!file::exists(clist)) {
        error("could not find '"+clist+"'");
        return 1;
    }

    vec1s mname;
    vec1s mfile;

    std::ifstream mlist(clist);
    uint_t l = 0;
    while (!mlist.eof()) {
        std::string line;
        std::getline(mlist, line);
        ++l;

        if (line.find_first_not_of(" \t") == line.npos) continue;
        line = trim(line);
        if (line[0] == '#') continue;
        vec1s spl = trim(split(line, "="));
        if (spl.size() != 2) {
            error("wrong format for l."+strn(l)+": '"+line+"'");
            note("expected: map_code_name = map_file");
            return 1;
        }

        spl[1] = dir+spl[1];
        if (!file::exists(spl[1])) {
            error("could not find '"+spl[1]+"'");
            return 1;
        }

        mname.push_back(spl[0]);
        mfile.push_back(spl[1]);
    }

    if (verbose) print("map list loaded successfully");

    if (thsize.size() > 1 && bands.size() > 1 && thsize.size() != bands.size()) {
        error("mismatch between number of bands and cutout sizes (",
            bands.size(), " vs. ", thsize.size(), ")");
        return 1;
    }

    vec1u idnm = where(!is_any_of(bands, mname));
    if (!idnm.empty()) {
        for (auto i : idnm) {
            warning("no band named '", bands[i], "', skipping");
        }

        bands = where(is_any_of(bands, mname));
    }

    uint_t ib = 0;
    for (uint_t b : range(mname)) {
        if (!bands.empty() && !is_any_of(mname[b], bands)) continue;

        ++ib;

        if (verbose) print(mname[b]);

        int_t hsize;
        if (finite(radius)) {
            // Convert radius to number of pixels
            double aspix;
            if (fits::get_pixel_size(mfile[b], aspix)) {
                hsize = ceil(radius/aspix);
                if (hsize < 10) hsize = 10;
            } else {
                warning("could not read WCS of '", mfile[b], "'");
                note("cutout size has been set to default (50 pixels)");
                hsize = 50;
            }
        } else if (thsize.size() > 1) {
            hsize = thsize[ib];
        } else if (thsize.size() == 1) {
            hsize = thsize[0];
        } else {
            // Use default cutout size
            hsize = 50;
        }

        qstack_params p;
        p.keep_nan = true;
        p.save_offsets = true;
        p.save_section = true;
        vec3d cube;
        vec1u ids;

        qstack_output qout = qstack(ra, dec, mfile[b], hsize, cube, ids, p);

        for (uint_t i : range(ids)) {
            fits::header nhdr = fits::read_header(mfile[b], qout.sect[i]);
            if (!fits::setkey(nhdr, "CRPIX1", hsize+1+qout.dx[i]) ||
                !fits::setkey(nhdr, "CRPIX2", hsize+1+qout.dy[i]) ||
                !fits::setkey(nhdr, "CRVAL1", ra[ids[i]]) ||
                !fits::setkey(nhdr, "CRVAL2", dec[ids[i]])) {
                warning("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
                note("WCS for the cutout will be wrong");
                note("parsing '"+mname[b]+"'");
            }

            std::string file_name = out+name[ids[i]]+mname[b]+".fits";
            if (verbose) print("writing ", file_name);
            fits::write(file_name, cube(i,_,_), nhdr);
        }

        // Create empty cutouts for non covered sources
        vec1u nids = complement(ra, ids);
        vec2d empty(2*hsize + 1, 2*hsize + 1);
        empty[_] = dnan;
        for (uint_t i : range(nids)) {
            fits::header nhdr = fits::read_header(mfile[b], 0);
            if (!fits::setkey(nhdr, "CRPIX1", hsize+1) ||
                !fits::setkey(nhdr, "CRPIX2", hsize+1) ||
                !fits::setkey(nhdr, "CRVAL1", ra[nids[i]]) ||
                !fits::setkey(nhdr, "CRVAL2", dec[nids[i]])) {
                warning("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
                note("WCS for the cutout will be wrong");
                note("parsing '"+mname[b]+"'");
            }

            std::string file_name = out+name[nids[i]]+mname[b]+".fits";
            if (verbose) print("writing ", file_name);
            fits::write(file_name, empty, nhdr);
        }
    }

    if (show.size() == 1 && show[0] == "1") {
        // No name specified: show all
        if (bands.empty()) {
            show = mname;
        } else {
            show = bands;
        }
    } else if (show.size() >= 1) {
        if (bands.empty()) {
            show = show[where(is_any_of(show, mname))];
        } else {
            show = show[where(is_any_of(show, bands))];
        }
    }

    if (show.size() > 0) {
        if (ra.size() != 1) {
            if (show.size() > 1) {
                warning("cannot display multiple bands for multiple sources with DS9");
                note("either ask for only one band to be shown, or extract sources one at a time");
            } else {
                spawn("ds9 -tile "+collapse(out+name[0]+show+".fits "));
            }
        } else {
            if (show_rgb) {
                if (show.size() > 3 && show_rgb) {
                    warning("cannot display more than 3 images at the same time in RGB mode, only "
                        "displaying the first 3");
                    show.resize(3);
                }

                vec1s chanels = {"red", "green", "blue"};
                chanels = chanels[uindgen(show.size())];

                spawn("ds9 -rgb "+collapse("-"+chanels+" "+out+name[0]+show+".fits "));
            } else {
                spawn("ds9 -tile "+collapse(out+name[0]+show+".fits "));
            }
        }
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("getgal v1.0");
    paragraph("usage: getgal field.param src=[...] out=... [options=...]");
    header("List of available command line options:");
    bullet("name", "[string] base name to give to individual cutouts (defalt: none)");
    bullet("dir", "[string] directory in which to look for maps (default: current)");
    bullet("verbose", "[flag] print some information about the process");
    bullet("show", "[string array] name of bands to show in DS9 (default: none)");
    bullet("show_rgb", "[flag] show then bands as RGB instead of tiles");
    bullet("radius", "[float] cutout radius in arcsec (if not provided, use the default cutout "
        "size from the parameter file)");
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
