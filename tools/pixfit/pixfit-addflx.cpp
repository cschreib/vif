#include <phypp.hpp>

int main(int argc, char* argv[]) {
    fits::table tbl(argv[1]);

    vec1u ids;
    tbl.read_column("id", ids);
    uint_t ngal = ids.size();

    vec1s bands;
    vec2f flux, flux_err;
    tbl.read_columns(fits::missing,
        "bands", bands, "flux", flux, "flux_err", flux_err
    );

    uint_t b = npos;
    uint_t id = npos;
    float flx = fnan, err = fnan;
    uint_t state = 0;
    for (int i : range(2, argc)) {
        vec1s spl = split(argv[i], "=");
        if (spl.size() == 2) {
            if (spl[0] == "band") {
                b = where_first(bands == spl[1]);
                if (b == npos) {
                    b = bands.size();
                    bands.push_back(spl[1]);
                    append<1>(flux, replicate(fnan, ngal, 1));
                    append<1>(flux_err, replicate(fnan, ngal, 1));
                }

                state = 0;
            } else {
                error("unknown parameter '", spl[0], "'");
                return 1;
            }
        } else if (spl.size() == 1) {
            if (state == 0) {
                ++state;
                if (!from_string(argv[i], id)) {
                    error("could not read ID in '", argv[i], "'");
                    return 1;
                }
            } else if (state == 1) {
                ++state;
                if (!from_string(argv[i], flx)) {
                    error("could not read flux in '", argv[i], "'");
                    return 1;
                }
            } else if (state == 2) {
                state = 0;
                if (!from_string(argv[i], err)) {
                    error("could not read flux error in '", argv[i], "'");
                    return 1;
                }

                uint_t p = where_first(ids == id);
                if (p == npos) {
                    error("could not find source with ID=", id);
                    return 1;
                }

                flux(p,b) = flx;
                flux_err(p,b) = err;
            }
        } else {
            error("ill formed parameter '", argv[i], "'");
            return 1;
        }
    }

    tbl.update_column("bands", bands);
    tbl.update_column("flux", flux);
    tbl.update_column("flux_err", flux_err);

    return 0;
}
