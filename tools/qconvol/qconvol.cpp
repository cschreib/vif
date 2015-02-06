#include <phypp.hpp>

void print_help() {
    using namespace format;

    print("qconvol v1.0");
    paragraph("usage: qconvol img.fits radius=1 out=output.fits");

    paragraph(
        "The program will convolve the provided image with a Gaussian beam of FWHM equal "
        "to 2 x radius [pixels] (or [arcsec] if the 'arcsec' keyword is provided), and "
        "save the result in a new FITS file."
    );
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string fimg = argv[1];
    std::string fout = "";
    double radius = 1.0;
    bool arcsec = false;

    read_args(argc-1, argv+1, arg_list(
        name(fout, "out"), radius, arcsec
    ));

    if (fout.empty()) {
        error("please provide the name of the output file in out=...");
        return 1;
    }

    if (arcsec) {
        double aspix;
        if (!fits::get_pixel_size(fimg, aspix)) {
            error("could not read the pixel size of this image, please provide the beam "
                "size in pixels and remove the 'arcsec' keyword");
            return 1;
        }

        radius /= aspix;
    }

    vec2d img;
    fits::header hdr;
    fits::read(fimg, img, hdr);

    radius /= 1.117;
    uint_t nk = ceil(30*radius);
    if (nk % 2 == 0) ++nk;

    vec2d beam = gaussian_profile({{nk,nk}}, radius);

    vec2d out = convolve2d(img, beam);
    fits::write(fout, out, hdr);

    return 0;
}
