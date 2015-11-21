#ifndef ASTRO_QSTACK_HPP
#define ASTRO_QSTACK_HPP

#include "phypp/astro.hpp"

namespace qstack_impl {
    struct image_workspace {
        int status = 0;
        fitsfile* fptr = nullptr;
        long width, height;
        fits::wcs astro;
        vec1d x, y;

        image_workspace() = default;

        explicit image_workspace(const std::string& file, const vec1d& ra, const vec1d& dec) {
            fits_open_image(&fptr, file.c_str(), READONLY, &status);
            fits::phypp_check_cfitsio(status, "cannot open file '"+file+"'");

            // Read the header as a string and read the WCS data
            char* hstr = nullptr;
            int nkeys  = 0;
            fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
            astro = fits::wcs(hstr);
            free(hstr);

            // Convert ra/dec to x/y
            fits::ad2xy(astro, ra, dec, x, y);

            // Get the dimensions of the image
            int naxis = 0;
            fits_get_img_dim(fptr, &naxis, &status);
            vec<1,long> naxes(naxis);
            fits_get_img_size(fptr, naxis, naxes.data.data(), &status);
            bool is2D = naxis == 2;
            if (is2D) {
                width = naxes[0];
                height = naxes[1];
            } else {
                uint_t found = 0;
                for (uint_t i : range(naxis)) {
                    if (naxes[i] > 1) {
                        if (found == 0) width = naxes[i];
                        if (found == 1) height = naxes[i];
                        ++found;
                    }
                }

                is2D = found == 2;
            }

            phypp_check(is2D, "cannot stack on image cubes (image dimensions: "+
                strn(naxes)+")");
        }

        image_workspace(const image_workspace&) = delete;
        image_workspace& operator=(const image_workspace&) = delete;

        image_workspace(image_workspace&& i) : status(i.status), fptr(i.fptr),
            width(i.width), height(i.height), astro(std::move(i.astro)),
            x(std::move(x)), y(std::move(y)) {
            i.fptr = nullptr;
        }

        image_workspace& operator=(image_workspace&& i) = delete;

        ~image_workspace() {
            if (fptr) fits_close_file(fptr, &status);
            fptr = nullptr;
        }
    };
}

struct qstack_params {
    bool keep_nan = false;
    bool save_offsets = false;
    bool save_section = false;
    bool verbose = false;
};

struct qstack_output {
    vec1d dx, dy;
    vec1u sect;
};

template<typename Type>
qstack_output qstack(const vec1d& ra, const vec1d& dec, const std::string& filename,
    uint_t hsize, vec<3,Type>& cube, vec1u& ids, qstack_params params = qstack_params()) {

    phypp_check(file::exists(filename), "cannot stack on inexistant file '"+filename+"'");
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    vec1s sects;
    if (end_with(filename, ".sectfits")) {
        sects = fits::read_sectfits(filename);
    } else {
        sects.push_back(filename);
    }

    std::vector<qstack_impl::image_workspace> imgs;
    imgs.reserve(sects.size());
    for (auto& s : sects) {
        imgs.emplace_back(s, ra, dec);
    }

    // Allocate memory to hold all the cutouts
    if (cube.empty()) {
        cube.dims[1] = cube.dims[2] = 2*hsize+1;
    }

    cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    ids.reserve(ids.size() + ra.size());

    vec1b found(ra.size());

    qstack_output out;
    if (params.save_offsets) {
        out.dx.reserve(ra.size());
        out.dy.reserve(ra.size());
    }

    if (params.save_section) {
        out.sect.reserve(ra.size());
    }

    // Loop over all images
    auto pg = progress_start(ra.size());
    for (uint_t iimg : range(imgs.size())) {
        auto& img = imgs[iimg];

        // Loop over all sources
        for (uint_t i = 0; i < ra.size(); ++i) {
            if (found[i]) continue;

            long p0[2] = {long(round(img.x[i]-hsize)), long(round(img.y[i]-hsize))};
            long p1[2] = {long(round(img.x[i]+hsize)), long(round(img.y[i]+hsize))};

            // Discard any source that falls out of the boundaries of the image
            if (p1[0] < 1 || p0[0] >= img.width || p1[1] < 1 || p0[1] >= img.height) {
                continue;
            }

            if (params.verbose) progress(pg);

            vec<2,Type> cut(2*hsize+1, 2*hsize+1);

            if (p0[0] < 1 || p1[0] >= img.width || p0[1] < 1 || p1[1] >= img.height) {
                // The source is overlapping with the edges of the map
                // Extract what we can
                cut[_] = fnan;

                long p0b[2] = {max(1, p0[0]),         max(1, p0[1])};
                long p1b[2] = {min(img.width, p1[0]), min(img.height, p1[1])};

                vec<2,Type> subcut(p1b[0]-p0b[0]+1, p1b[1]-p0b[1]+1);

                Type null = fnan;
                int anynul = 0;
                long inc[2] = {1, 1};
                fits_read_subset(img.fptr, fits::traits<Type>::ttype, p0b, p1b, inc, &null,
                    subcut.data.data(), &anynul, &img.status);

                long x0 = p0b[0]-p0[0];
                long x1 = p1b[0]-p0[0];
                long y0 = p0b[1]-p0[1];
                long y1 = p1b[1]-p0[1];
                cut(y0-_-y1,x0-_-x1) = subcut;
            } else {
                // The source is fully covered, easy
                Type null = fnan;
                int anynul = 0;
                long inc[2] = {1, 1};
                fits_read_subset(img.fptr, fits::traits<Type>::ttype, p0, p1, inc, &null,
                    cut.data.data(), &anynul, &img.status);
            }

            found[i] = true;

            // Discard any source that contains a bad pixel (either infinite or NaN)
            if (!params.keep_nan && count(!is_finite(cut)) != 0) {
                continue;
            }

            ids.push_back(i);
            cube.push_back(cut);

            if (params.save_offsets) {
                out.dx.push_back(img.x[i] - round(img.x[i]));
                out.dy.push_back(img.y[i] - round(img.y[i]));
            }

            if (params.save_section) {
                out.sect.push_back(iimg);
            }
        }
    }

    return out;
}

template<typename Type>
qstack_output qstack(const vec1d& ra, const vec1d& dec, const std::string& ffile,
    const std::string& wfile, uint_t hsize, vec<3,Type>& cube, vec<3,Type>& wcube,
    vec1u& ids, qstack_params params = qstack_params()) {

    qstack_output out;
    if (params.save_offsets) {
        out.dx.reserve(ra.size());
        out.dy.reserve(ra.size());
    }

    if (params.save_section) {
        out.sect.reserve(ra.size());
    }

    phypp_check(file::exists(ffile), "cannot stack on inexistant file '"+ffile+"'");
    phypp_check(file::exists(wfile), "cannot stack on inexistant file '"+wfile+"'");
    phypp_check(ra.size() == dec.size(), "need ra.size() == dec.size()");

    if (end_with(ffile, ".sectfits") || end_with(wfile, ".sectfits")) {
        uint_t norig = ids.size();
        qstack(ra, dec, ffile, hsize, cube, ids, params);
        uint_t nsci = ids.size() - norig;

        vec1u tids; // trash weight IDs: they are already collected in 'ids'
        qstack_params tparams = params; // trash output, for the same reason
        qstack(ra, dec, wfile, hsize, wcube, tids, tparams);
        uint_t nwht = tids.size();

        // Make sure that all sources that were extracted in one file match those that
        // are extracted in the other (could only keep the union of the two, WIP)
        phypp_check(nsci == nwht, "some sources are covered on '", ffile, "' but not '", wfile, "'");

        return out;
    }

    // Open the FITS file
    fitsfile* fptr;
    fitsfile* wfptr;
    int status = 0;

    fits_open_image(&fptr, ffile.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+ffile+"'");
    fits_open_image(&wfptr, wfile.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+wfile+"'");

    // Read the header as a string and read the WCS data
    char* hstr = nullptr;
    int nkeys  = 0;
    fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
    fits::wcs astro(hstr);
    free(hstr);

    // Get the dimensions of the image
    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);
    vec<1,long> naxes(naxis);
    fits_get_img_size(fptr, naxis, naxes.data.data(), &status);
    bool is2D = naxis == 2;
    long width, height;
    if (is2D) {
        width = naxes[0];
        height = naxes[1];
    } else {
        uint_t found = 0;
        for (uint_t i : range(naxis)) {
            if (naxes[i] > 1) {
                if (found == 0) width = naxes[i];
                if (found == 1) height = naxes[i];
                ++found;
            }
        }

        is2D = found == 2;
    }

    phypp_check(is2D, "cannot stack on image cubes (image dimensions: "+strn(naxes)+")");

    vec<1,long> wnaxes(naxis);
    fits_get_img_size(wfptr, naxis, wnaxes.data.data(), &status);
    phypp_check(naxes[0] == wnaxes[0] && naxes[1] == wnaxes[1], "image and weight map do not match");

    // Convert ra/dec to x/y
    vec1d x, y;
    fits::ad2xy(astro, ra, dec, x, y);

    // Allocate memory to hold all the cutouts
    if (cube.empty()) {
        cube.dims[1] = cube.dims[2] = 2*hsize+1;
    }
    if (wcube.empty()) {
        wcube.dims[1] = wcube.dims[2] = 2*hsize+1;
    }

    cube.reserve(cube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    wcube.reserve(wcube.size() + (2*hsize+1)*(2*hsize+1)*ra.size());
    ids.reserve(ids.size() + ra.size());

    // Loop over all sources
    for (uint_t i = 0; i < ra.size(); ++i) {
        long p0[2] = {long(round(x[i]-hsize)), long(round(y[i]-hsize))};
        long p1[2] = {long(round(x[i]+hsize)), long(round(y[i]+hsize))};

        // Discard any source that falls out of the boundaries of the image
        if (p0[0] < 1 || p1[0] >= width || p0[1] < 1 || p1[1] >= height) {
            continue;
        }

        vec<2,Type> cut(2*hsize+1, 2*hsize+1);
        vec<2,Type> wcut(2*hsize+1, 2*hsize+1);

        Type null = fnan;
        int anynul = 0;
        long inc[2] = {1, 1};

        fits_read_subset(fptr,  fits::traits<Type>::ttype, p0, p1, inc, &null,
            cut.data.data(),  &anynul, &status);
        fits_read_subset(wfptr, fits::traits<Type>::ttype, p0, p1, inc, &null,
            wcut.data.data(), &anynul, &status);

        // Discard any source that contains a bad pixel (either infinite or NaN)
        if (!params.keep_nan && count(!is_finite(cut) || !is_finite(wcut)) != 0) {
            continue;
        }

        ids.push_back(i);
        cube.push_back(cut);
        wcube.push_back(wcut);

        if (params.save_offsets) {
            out.dx.push_back(x[i] - round(x[i]));
            out.dy.push_back(y[i] - round(y[i]));
        }

        if (params.save_section) {
            out.sect.push_back(0);
        }
    }

    fits_close_file(fptr, &status);
    fits_close_file(wfptr, &status);

    return out;
}

template<typename Type>
vec<2,rtype_t<Type>> qstack_mean(const vec<3,Type>& fcube) {
    return partial_mean(0, fcube);
}

template<typename TypeF, typename TypeW>
vec<2,rtype_t<TypeF>> qstack_mean(const vec<3,TypeF>& fcube, const vec<3,TypeW>& wcube) {
    return partial_total(0, fcube*wcube)/partial_total(0, wcube);
}

template<typename Type>
vec<2,rtype_t<Type>> qstack_median(const vec<3,Type>& fcube) {
    return partial_median(0, fcube);
}

template<typename Type, typename TypeS, typename F>
void qstack_bootstrap(const vec<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed, F&& func) {

    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        auto tfcube = fcube(ids,_,_).concretise();
        func(tfcube);
    }
}

template<typename TypeF, typename TypeW, typename TypeS, typename F>
void qstack_bootstrap(const vec<3,TypeF>& fcube, const vec<3,TypeW>& wcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed, F&& func) {

    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        auto tfcube = fcube(ids,_,_).concretise();
        auto twcube = wcube(ids,_,_).concretise();
        func(tfcube, twcube);
    }
}

template<typename Type>
vec<3,rtype_t<Type>> qstack_bootstrap_apply_id_(const vec1u& ids, const vec<3,Type>& cube) {
    return cube(ids,_,_).concretise();
}

template<typename Type, typename ... Args>
uint_t qstack_bootstrap_get_size_(const vec<3,Type>& cube, const Args& ... cubes) {
    return cube.dims[0];
}

template<typename TypeS, typename F, typename ... Args>
void qstack_bootstrap(uint_t nbstrap, uint_t nsel, TypeS& seed, F&& func, const Args& ... cubes) {
    const uint_t nsrc = qstack_bootstrap_get_size_(cubes...);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, nsrc-1, nsel);
        func(qstack_bootstrap_apply_id_(ids, cubes)...);
    }
}

template<typename Type, typename TypeS>
vec<3,rtype_t<Type>> qstack_mean_bootstrap(const vec<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed) {

    vec<3,rtype_t<Type>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_mean(fcube(ids,_,_)));
    }

    return bs;
}

template<typename TypeF, typename TypeW, typename TypeS>
vec<3,rtype_t<TypeF>> qstack_mean_bootstrap(const vec<3,TypeF>& fcube,
    const vec<3,TypeW>& wcube, uint_t nbstrap, uint_t nsel, TypeS& seed) {

    vec<3,rtype_t<TypeF>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_mean(fcube(ids,_,_), wcube(ids,_,_)));
    }

    return bs;
}

template<typename Type, typename TypeS>
vec<3,rtype_t<Type>> qstack_median_bootstrap(const vec<3,Type>& fcube, uint_t nbstrap,
    uint_t nsel, TypeS& seed) {

    vec<3,rtype_t<Type>> bs;
    bs.reserve(nbstrap*fcube.dims[1]*fcube.dims[2]);
    for (uint_t i = 0; i < nbstrap; ++i) {
        vec1u ids = randomi(seed, 0, fcube.dims[0]-1, nsel);
        bs.push_back(qstack_median(fcube(ids,_,_)));
    }

    return bs;
}

#endif
