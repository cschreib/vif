#include <vif.hpp>

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 4) {
        print_help();
        return 0;
    }

    fits::header hdr = fits::read_header(argv[1]);
    astro::wcs astro(hdr);

    uint_t nc = argc - 2;
    if (nc % 2 != 0) {
        error("expected an even number of coordinates, RA and Dec");
        return 1;
    }

    nc /= 2;

    vec1d ra(nc), dec(nc);
    for (uint_t c : range(nc)) {
        std::string sra = argv[2*c+2];
        std::string sdec = argv[2*c+3];

        if (find(sra, ":") == npos) {
            if (!from_string(sra, ra[c])) {
                error("could not convert '"+sra+"' to a coordinate");
                return 1;
            }
            if (!from_string(sdec, dec[c])) {
                error("could not convert '"+sdec+"' to a coordinate");
                return 1;
            }
        } else {
            sex2deg(sra, sdec, ra[c], dec[c]);
        }
    }

    vec1d px, py;
    astro::ad2xy(astro, ra, dec, px, py);

    for (uint_t c : range(nc)) {
        print(px[c], "    ", py[c]);
    }

    return 0;
}

void print_help() {
    print("radec2pix v1.0");
    print("usage: radec2pix img.fits RA1 Dec1 RA2 Dec2 ...");
}
