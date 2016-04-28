#include "pixfit-common.hpp"

int phypp_main(int argc, char* argv[]) {
    vec2d regs;
    vec1s text;

    bool physical = false;
    if (!read_ds9_region_circles(argv[1], regs, text, physical)) {
        return 1;
    }

    if (physical) {
        error("must be WCS coordinates");
        return 1;
    }

    vec1u id;
    vec1f z, m;
    for (std::string txt : text) {
        vec1s spl = split(txt,",");
        if (spl.size() != 3) {
            error("need three elements: ID, z and M*");
            return 1;
        }

        uint_t tid;
        float tz, tm;
        if (!from_string(spl[0], tid) || !from_string(spl[1], tz) || !from_string(spl[2], tm)) {
            error("could not convert ID, z and M* into numbers");
            return 1;
        }

        id.push_back(tid);
        z.push_back(tz);
        m.push_back(tm);
    }

    fits::write_table(file::remove_extension(argv[1])+".fits",
        "id", id, "ra", vec1d{regs(_,0)}, "dec", vec1d{regs(_,1)}, "z", z, "m", m
    );

    return 0;
}
