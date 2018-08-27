#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    fits::table tbl(argv[1]);

    vec1u ids;
    tbl.read_column("id", ids);
    uint_t ngal = ids.size();

    vec1f tdust_min, tdust_max;
    tbl.read_columns(fits::missing, "tdust_min", tdust_min, "tdust_max", tdust_max);

    if (tdust_min.empty()) {
        tdust_min = replicate(fnan, ids.size());
    }
    if (tdust_max.empty()) {
        tdust_max = replicate(fnan, ids.size());
    }

    uint_t id = npos;
    float mi = fnan, ma = fnan;

    auto add_cst = [&]() {
        if (is_finite(mi) || is_finite(ma)) {
            if (id == npos) {
                error("please specify the galaxy ID before min=... and max=...");
                return false;
            }

            uint_t p = where_first(ids == id);
            if (p == npos) {
                error("could not find source with ID=", id);
                return false;
            }

            tdust_min[p] = mi;
            tdust_max[p] = ma;

            mi = fnan; ma = fnan;
        }

        return true;
    };

    for (int i : range(2, argc)) {
        vec1s spl = split(argv[i], "=");
        if (spl.size() == 2) {
            if (spl[0] == "min") {
                if (!from_string(spl[1], mi)) {
                    error("could not read minimum temperature in '", argv[i], "'");
                    return 1;
                }
            } else if (spl[0] == "max") {
                if (!from_string(spl[1], ma)) {
                    error("could not read maximum temperature in '", argv[i], "'");
                    return 1;
                }
            } else {
                error("unknown parameter '", spl[0], "'");
                return 1;
            }
        } else if (spl.size() == 1) {
            if (!add_cst()) {
                return 1;
            }

            if (!from_string(argv[i], id)) {
                error("could not read ID in '", argv[i], "'");
                return 1;
            }
        } else {
            error("ill formed parameter '", argv[i], "'");
            return 1;
        }
    }

    if (!add_cst()) {
        return 1;
    }

    tbl.update_columns("tdust_min", tdust_min, "tdust_max", tdust_max);

    return 0;
}
