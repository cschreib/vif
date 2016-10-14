#include "pixfit-common.hpp"
#include <phypp/astro/qstack.hpp>

int phypp_main(int argc, char* argv[]) {
    std::vector<map_info> maps;
    if (!read_maps(argv[1], maps)) return 1;

    double radius = 30.0;
    double ra = dnan, dec = dnan;
    read_args(argc-1, argv+1, arg_list(radius, ra, dec));

    if (!is_finite(ra) || !is_finite(dec)) {
        error("please provide the center coordinate with ra=... dec=...");
        return 1;
    }

    for (auto map : maps)
    for (std::string type : {"sci", "err"}) {
        std::string file = (type == "sci" ? map.img : map.err);

        int_t hsize;
        double aspix;
        if (astro::get_pixel_size(file, aspix)) {
            hsize = ceil(radius/aspix);
        } else {
            return 1;
        }

        qstack_params p;
        p.keep_nan = true;
        p.save_offsets = true;
        vec3d cube;
        vec1u ids;

        qstack_output qout = qstack(vec1d{ra}, vec1d{dec}, file, hsize, cube, ids, p);

        if (ids.empty()) {
            // Create empty cutouts for non covered sources
            cube.resize(1, 2*hsize + 1, 2*hsize + 1);
            cube[_] = dnan;
            qout.dx = {0.0};
            qout.dy = {0.0};
        }

        // Build new header
        fits::header nhdr = astro::filter_wcs(fits::read_header(file));
        if (!fits::setkey(nhdr, "CRPIX1", hsize+1+qout.dx[0]) ||
            !fits::setkey(nhdr, "CRPIX2", hsize+1+qout.dy[0]) ||
            !fits::setkey(nhdr, "CRVAL1", ra) ||
            !fits::setkey(nhdr, "CRVAL2", dec)) {
            error("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
            note("WCS for the cutout will be wrong");
            note("parsing band '"+map.band+"'");
            return 1;
        }

        fits::write(map.band+"-"+type+".fits", cube(0,_,_), nhdr);
    }

    return 0;
}
