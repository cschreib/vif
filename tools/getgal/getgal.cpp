#include <vif.hpp>
#include <vif/astro/qstack.hpp>

void print_help();

struct image_t {
    // Read from parameter file
    std::string filename;
    std::string short_name;
    std::string band;
    std::string psffile;
    double seeing = dnan;
    double zero_point = 23.9;

    // Work variables
    bool used = false;
};

bool read_image_list(const std::string& filename, vec<1,image_t>& imgs) {
    std::string idir = file::get_directory(filename);

    std::string line;
    std::ifstream file(filename);
    uint_t l = 0;

    while (std::getline(file, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto eqpos = line.find('=');
        if (eqpos == line.npos) {
            image_t img;
            if (line[0] == '/') {
                img.filename = line;
            } else {
                img.filename = idir+line;
            }
            imgs.push_back(std::move(img));
        } else {
            std::string key = trim(line.substr(0, eqpos));
            std::string val = trim(line.substr(eqpos+1));
            if (key == "band") {
                imgs.back().band = val;
            } else if (key == "seeing") {
                if (!from_string(val, imgs.back().seeing)) {
                    error("could not read seeing value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "zero_point") {
                if (!from_string(val, imgs.back().zero_point)) {
                    error("could not read zero point value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "short_name") {
                imgs.back().short_name = val;
            } else if (key == "psf") {
                if (!val.empty()) {
                    if (val[0] == '/') {
                        imgs.back().psffile = val;
                    } else {
                        imgs.back().psffile = idir+val;
                    }
                }
            } else {
                warning("unknown parameter '", key, "', ignored");
                warning("at ", filename, ":", l);
            }
        }
    }

    return true;
}

int vif_main(int argc, char* argv[]) {
    vec1s tsrc;
    std::string out = "";
    std::string nbase = "";
    std::string dir = "";
    double radius = dnan;
    vec1u thsize;
    vec1s bands;
    bool no_zero_point = false;
    bool make_list = false;
    bool verbose = false;

    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string clist = argv[1];

    read_args(argc-1, argv+1, arg_list(
        name(tsrc, "src"), out, name(nbase, "name"), dir, verbose, radius,
        name(thsize, "hsize"), bands, no_zero_point, make_list
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
    if (tsrc.size() == 1 && ends_with(tsrc[0], ".fits")) {
        struct {
            vec1d ra, dec;
            vec1s name;
            vec1u id;
        } tmp;

        fits::read_table_loose(tsrc[0], ftable(tmp.ra, tmp.dec, tmp.name, tmp.id));
        if (tmp.ra.empty() || tmp.dec.empty()) {
            error("missing RA and Dec coordinates in this FITS file");
            return 1;
        }

        if (tmp.name.empty()) {
            if (tmp.id.empty()) {
                tmp.id = uindgen(tmp.ra.size());
            }

            name = to_string_vector(tmp.id);
        } else {
            name = tmp.name;
        }

        name += "_";
        if (!nbase.empty()) name = nbase + "_" + name;

        ra = tmp.ra;
        dec = tmp.dec;
    } else if (tsrc.size() == 1 && ends_with(tsrc[0], ".reg")) {
        ascii::read_table(tsrc[0], ra, dec);
        name = to_string_vector(uindgen(ra.size())) + "_";
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

    vec<1,image_t> imgs;

    if (ends_with(clist, ".fits")) {
        // Single FITS file
        image_t img;
        img.filename = clist;
        img.short_name = split(erase_end(clist, ".fits"), "/").back();
        imgs.push_back(std::move(img));
    } else {
        // Map list
        if (!read_image_list(clist, imgs)) {
            return 1;
        }

        if (verbose) print("map list loaded successfully");
    }

    if (thsize.size() > 1 && imgs.size() > 1 && thsize.size() != imgs.size()) {
        error("mismatch between number of bands and cutout sizes (",
            bands.size(), " vs. ", thsize.size(), ")");
        return 1;
    }

    // vec1u idnm = where(!is_any_of(bands, mname));
    // if (!idnm.empty()) {
    //     for (auto i : idnm) {
    //         warning("no band named '", bands[i], "', skipping");
    //     }

    //     bands = bands[where(is_any_of(bands, mname))];
    // }

    vec1b covered(bands.size());
    uint_t ib = 0;
    for (auto& img : imgs) {
        if (!bands.empty()) {
            uint_t ibb = where_first(bands == img.short_name);
            if (ibb == npos) continue;
            covered[ibb] = true;
        }

        img.used = true;

        ++ib;

        if (verbose) print(img.short_name);

        // Find HDU containing data
        uint_t hdu = 0; {
            std::string tfile;
            if (ends_with(img.filename, ".sectfits")) {
                vec1s sects = fits::read_sectfits(img.filename);
                tfile = sects[0];
            } else {
                tfile = img.filename;
            }

            fits::input_image iimg(tfile);
            for (uint_t i : range(iimg.hdu_count())) {
                iimg.reach_hdu(i);
                if (iimg.axis_count() != 0) {
                    hdu = i;
                    break;
                }
            }
        }

        int_t hsize;
        if (is_finite(radius)) {
            // Convert radius to number of pixels
            double aspix;
            if (astro::get_pixel_size(img.filename, aspix)) {
                hsize = ceil(radius/aspix);
                if (hsize < 10) hsize = 10;
            } else {
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
        p.verbose = verbose;
        vec3d cube;
        vec1u ids;

        qstack_output qout = qstack(ra, dec, img.filename, hsize, cube, ids, p);

        if (!no_zero_point && is_finite(img.zero_point)) {
            // Apply zero point to convert map to uJy
            cube *= e10(0.4*(23.9 - img.zero_point));
            img.zero_point = 23.9;
        }

        for (uint_t i : range(ids)) {
            fits::header nhdr = astro::filter_wcs(fits::read_header_sectfits_hdu(
                img.filename, qout.sect[i], hdu
            ));

            if (!fits::setkey(nhdr, "CRPIX1", hsize+1+qout.dx[i]) ||
                !fits::setkey(nhdr, "CRPIX2", hsize+1+qout.dy[i]) ||
                !fits::setkey(nhdr, "CRVAL1", ra[ids[i]]) ||
                !fits::setkey(nhdr, "CRVAL2", dec[ids[i]])) {
                warning("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
                note("WCS for the cutout will be wrong");
                note("parsing '"+img.short_name+"'");
            }

            std::string filename = out+name[ids[i]]+img.short_name+".fits";

            // Make sure that we are not going to overwrite one of the images
            if (filename == img.filename) {
                error("this operation would overwrite the image '", filename, "'");
                note("aborting");
                return 1;
            }

            if (verbose) print("writing ", filename);
            fits::write(filename, cube(i,_,_), nhdr);
        }

        // Create empty cutouts for non covered sources
        vec1u nids = complement(ra, ids);
        vec2d empty(2*hsize + 1, 2*hsize + 1);
        empty[_] = dnan;
        for (uint_t i : range(nids)) {
            fits::header nhdr = astro::filter_wcs(fits::read_header_sectfits_hdu(
                img.filename, 0, hdu
            ));

            if (!fits::setkey(nhdr, "CRPIX1", hsize+1) ||
                !fits::setkey(nhdr, "CRPIX2", hsize+1) ||
                !fits::setkey(nhdr, "CRVAL1", ra[nids[i]]) ||
                !fits::setkey(nhdr, "CRVAL2", dec[nids[i]])) {
                warning("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
                note("WCS for the cutout will be wrong");
                note("parsing '"+img.short_name+"'");
            }

            std::string filename = out+name[nids[i]]+img.short_name+".fits";

            // Make sure that we are not going to overwrite one of the images
            if (filename == img.filename) {
                error("this operation would overwrite the image '", filename, "'");
                note("aborting");
                return 1;
            }

            if (verbose) print("writing ", filename);
            fits::write(filename, empty, nhdr);
        }
    }

    if (count(!covered) > 0) {
        for (auto i : where(!covered)) {
            warning("no band named '", bands[i], "', skipping");
        }
    }

    if (make_list) {
        for (uint_t i : range(name)) {
            std::ofstream olist(out+name[i]+"images.param");
            for (auto& img : imgs) {
                if (!img.used) continue;

                olist << name[i]+img.short_name+".fits\n";
                if (!img.band.empty())         olist << "band=" << img.band << "\n";
                if (is_finite(img.seeing))     olist << "seeing=" << img.seeing << "\n";
                if (is_finite(img.zero_point)) olist << "zero_point=" << img.zero_point << "\n";
                olist << "\n";
            }
        }
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("getgal v1.0");
    paragraph("usage: getgal field.param src=[...] out=... [options=...]");
    header("List of available command line options:");
    bullet("name", "[string] base name to give to individual cutouts (defalt: none)");
    bullet("dir", "[string] directory in which to look for maps (default: current)");
    bullet("verbose", "[flag] print some information about the process");
    bullet("make_list", "[flag] output a file containing the list of images, the "
        "corresponding band and other information to be used for flux extraction");
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
