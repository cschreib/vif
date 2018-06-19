#include <phypp.hpp>

#define check(t, s) { \
    std::string st = to_string(t); \
    if (st == s) print("  checked: "+st); \
    else         print("  failed: "+std::string(#t)+" = "+st+" != "+s); \
    assert(st == s); \
}

int phypp_main(int argc, char* argv[]) {
    vec1u id = uindgen(100);
    vec1s sid = to_string_vector(reverse(id));

    {
        fits::write_table("col_1d.fits", ftable(id, sid));
    }

    {
        vec1u tid;
        vec1s tsid;
        fits::read_table("col_1d.fits", "id", tid, "sid", tsid);

        check(count(tid != id), "0");
        check(count(tsid != sid), "0");
    }

    {
        fits::output_table otbl("row_1d.fits");
        otbl.set_format(fits::output_format::row_oriented);
        otbl.write_columns(ftable(id, sid));
    }

    {
        vec1u tid;
        vec1s tsid;
        fits::read_table("row_1d.fits", "id", tid, "sid", tsid);

        check(count(tid != id), "0");
        check(count(tsid != sid), "0");
    }

    vec3u id2d = uindgen(3,5,100);
    vec3s sid2d = to_string_vector(id2d);

    {
        fits::write_table("col_2d.fits", ftable(id2d, sid2d));
    }

    {
        vec3u tid2d;
        vec3s tsid2d;
        fits::read_table("col_2d.fits", "id2d", tid2d, "sid2d", tsid2d);

        check(count(tid2d != id2d), "0");
        check(count(tsid2d != sid2d), "0");
    }

    {
        fits::output_table otbl("row_2d.fits");
        otbl.set_format(fits::output_format::row_oriented);
        otbl.write_columns(ftable(id2d, sid2d));
    }

    {
        vec3u tid2d;
        vec3s tsid2d;
        fits::read_table("row_2d.fits", "id2d", tid2d, "sid2d", tsid2d);

        check(count(tid2d != id2d), "0");
        check(count(tsid2d != sid2d), "0");
    }

    return 0;
}
