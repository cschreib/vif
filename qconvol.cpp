#include <vif.hpp>

using namespace vif;

void print_help() {
    using namespace terminal_format;

    print("qconvol v1.0");
    paragraph("usage: qconvol img.fits radius=1 kernel=\"\" out=output.fits");

    paragraph(
        "The program will convolve the provided image with a Gaussian beam of FWHM equal "
        "to 2 x radius [pixels] (or [arcsec] if the 'arcsec' keyword is provided), and "
        "save the result in a new FITS file.\n\n"
        "Alternatively, one may provide a 'kernel' image in FITS format which will be "
        "used directly to perform the convolution."
    );
}

int vif_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string fimg = argv[1];
    std::string fout = "";
    std::string kernel_file = "";
    double radius = 1.0;
    bool arcsec = false;

    read_args(argc-1, argv+1, arg_list(
        name(fout, "out"), radius, arcsec, name(kernel_file, "kernel")
    ));

    if (fout.empty()) {
        error("please provide the name of the output file in out=...");
        return 1;
    }


    vec2d beam;
    if (!kernel_file.empty()) {
        fits::read(kernel_file, beam);

        // Trim the PSF, make sure that the peak is at the center, and that the dimensions
        // are odds
        vec1i idm = mult_ids(beam, max_id(beam));
        int_t hsize = 0;
        int_t imax = std::min(
            std::min(idm[0], int_t(beam.dims[0])-1-idm[0]),
            std::min(idm[1], int_t(beam.dims[1])-1-idm[1])
        );

        for (int_t i = 1; i <= imax; ++i) {
            if (beam(idm[0]-i,idm[1]) == 0.0 &&
                beam(idm[0]+i,idm[1]) == 0.0 &&
                beam(idm[0],idm[1]-i) == 0.0 &&
                beam(idm[0],idm[1]+i) == 0.0) {
                hsize = i;
                break;
            }
        }

        if (hsize == 0) hsize = imax;
        //beam = subregion(beam, {idm[0]-hsize, idm[1]-hsize, idm[0]+hsize, idm[1]+hsize});
        beam = astro::subregion(beam, {idm[0]-hsize, idm[1]-hsize, idm[0]+hsize, idm[1]+hsize});
    } else {
        if (arcsec) {
            double aspix;
            if (!astro::get_pixel_size(fimg, aspix)) {
                error("could not read the pixel size of this image, please provide the beam "
                    "size in pixels and remove the 'arcsec' keyword");
                return 1;
            }

            radius /= aspix;
        }

        radius /= 1.117;
        uint_t nk = ceil(30*radius);
        if (nk % 2 == 0) ++nk;

        //beam = gaussian_profile({{nk,nk}}, radius);
        beam = astro::gaussian_profile({{nk,nk}}, radius);
        beam /= total(beam);
    }

    vec2d img;
    fits::header hdr;
    fits::read(fimg, img, hdr);

    //vec2d out = convolve2d(img, beam);
    vec2d out = astro::convolve2d(img, beam);
    fits::write(fout, out, hdr);

    return 0;
}
