#include "pixfit-common.hpp"
#include <phypp/astro/qxmatch.hpp>

int main(int argc, char* argv[]) {
    // Read search circles
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

    // Read catalog
    fits::input_table tbl(argv[2]);

    std::string cra = "pos.ra";
    std::string cdec = "pos.dec";
    std::string cid = "id.cls";
    std::string cz = "z.mp.phot";
    std::string cm = "m.mp.del";
    uint_t nsrc = 100;
    double mref = 9, zref = 1.0, mmin = 10.0;
    read_args(argc-3, argv+3, arg_list(
        name(cra, "ra"), name(cdec, "dec"), name(cid, "id"), name(cz, "z"), name(cm, "m"),
        nsrc, mmin, mref, zref
    ));

    vec1u id;
    vec1d ra, dec;
    vec1f z, m;
    tbl.read_columns(fits::narrow, cra, ra, cdec, dec, cid, id, cz, z, cm, m);

    // Filter catalog
    auto cosmo = get_cosmo("std");
    vec1d d = lumdist(z, cosmo);
    double dref = lumdist(zref, cosmo);

    vec1u sid = where(m > mref + 2*log10(d/dref) || m > mmin);
    ra = ra[sid];
    dec = dec[sid];
    id = id[sid];
    z = z[sid];
    m = m[sid];

    // Cross match
    vec1d sra = regs(_,0);
    vec1d sdec = regs(_,1);

    qxmatch_params p; p.nth = nsrc; p.brute_force = true;
    auto res = qxmatch(sra, sdec, ra, dec, p);

    // Build region file
    std::ofstream out(argv[3]);
    out << "# Region file format: DS9 version 4.1\n";
    out << "global color=green dashlist=8 3 width=2 font=\"helvetica 10 normal roman\""
        "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n";
    out << "fk5\n";

    for (uint_t i : range(sra))
    for (uint_t j : range(nsrc)) {
        if (res.d(j,i) > regs(i,2)) break;

        uint_t k = res.id(j,i);
        std::string rra, rdec;
        deg2sex(ra[k], dec[k], rra, rdec);
        out << "circle(" << rra << "," << rdec << ",0.5\") # width=3 text={"
            << id[k] << ", " << z[k] << ", " << m[k] << "}\n";
    }

    return 0;
}
