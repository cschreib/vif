#include "pixfit-common.hpp"

int phypp_main(int argc, char* argv[]) {
    std::string map_file;
    std::string cat_file;
    double snr_min = dnan;

    read_args(argc, argv, arg_list(name(map_file, "maps"), name(cat_file, "cat"), snr_min));

    bool bad = false;
    if (map_file.empty()) {
        error("please provide the map list file in maps=...");
        bad = true;
    }
    if (cat_file.empty()) {
        error("please provide the flux catalog in cat=...");
        bad = true;
    }

    if (bad) return 1;

    std::vector<map_info> maps;
    if (!read_maps(map_file, maps)) return 1;

    vec1d ra, dec;
    vec2f flux;
    vec1s bands;
    fits::read_table(cat_file, ftable(ra, dec, flux, bands));

    vec1b sel = replicate(true, ra.size());
    if (is_finite(snr_min)) {
        vec1d lir, lir_err;
        fits::read_table(cat_file, ftable(lir, lir_err));
        sel = lir/lir_err > snr_min;
    }

    for (auto& map : maps) {
        // Find band in catalog
        uint_t b = where_first(bands == map.band);
        if (b == npos) {
            warning("no band named '", map.band, "' in the provided catalog");
            continue;
        }

        // Read map
        vec2d img;
        fits::header hdr;
        fits::read(map.band+"-sci.fits", img, hdr);

        // Read PSF
        int_t hsize;
        vec2d psf = read_psf(map, hsize);

        // Get pixel coordinates
        vec1d x, y;
        fits::ad2xy(fits::extast(hdr), ra, dec, x, y);
        x -= 1; y -= 1;

        // Select sources which fall on the map and have a flux measurement
        vec1u ids = where(is_finite(flux(_,b)) && sel &&
            x >= -hsize && x < img.dims[1]+hsize &&
            y >= -hsize && y < img.dims[0]+hsize);

        x = x[ids]; y = y[ids];
        vec1f flx = flux(ids,b);
        vec1i ix = round(x), iy = round(y);
        vec1d dx = x - ix, dy = y - iy;

        // Remove the sources from the map
        for (uint_t i : range(x)) {
            vec1f tpsf = flatten(translate(psf, dy[i], dx[i]));
            vec1u idi, idp;
            subregion(img, {iy[i]-hsize, ix[i]-hsize, iy[i]+hsize, ix[i]+hsize}, idi, idp);

            img[idi] -= (flx[i]/map.fconv)*tpsf[idp];
        }

        // Save the residual
        fits::write(map.band+"-fullres.fits", img, hdr);
    }

    return 0;
}
