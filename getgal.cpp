#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    vec1s tsrc;
    std::string out = "";
    std::string nbase = "";
    std::string dir = "";
    vec1s show;
    bool verbose = false;

    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string clist = argv[1];

    read_args(argc-1, argv+1, arg_list(
        name(tsrc, "src"), out, name(nbase, "name"), dir, verbose, show
    ));

    if (!dir.empty()) {
        dir = file::directorize(dir);
    }

    if (!out.empty()) {
        out = file::directorize(out);
        if (!file::mkdir(out)) {
            warning("could not create directory '"+out+"'");
        }
    }

    vec1d ra, dec;
    vec1s name;
    if (tsrc.size() == 1 && end_with(tsrc[0], ".fits")) {
        struct {
            vec1d ra, dec;
            vec1s name;
            vec1u id;
        } tmp;

        fits::read_table_loose(tsrc[0], tmp);
        if (tmp.ra.empty() || tmp.dec.empty()) {
            error("missing RA and Dec coordinates in this FITS file");
            return 1;
        }

        if (tmp.name.empty()) {
            if (tmp.id.empty()) {
                tmp.id = uindgen(tmp.ra.size());
            }

            name = strna(tmp.id);
        } else {
            name = tmp.name;
        }

        name += "_";
        if (!nbase.empty()) name = nbase + "_" + name;

        ra = tmp.ra;
        dec = tmp.dec;
    } else if (tsrc.size() == 2) {
        name.resize(1);
        if (!nbase.empty()) name[0] = nbase + "_";
        ra.resize(1);
        dec.resize(1);
        if (!from_string(tsrc[0], ra[0])) {
            error("could not convert '"+tsrc[0]+"' to a coordinate");
            return 1;
        }
        if (!from_string(tsrc[1], dec[0])) {
            error("could not convert '"+tsrc[1]+"' to a coordinate");
            return 1;
        }
    } else {
        error("incorrect format for 'src' parameter");
        print_help();
        return 1;
    }

    if (!file::exists(clist)) {
        error("could not find '"+clist+"'");
        return 1;
    }

    vec1s mname;
    vec1s mfile;
    vec1u msize;

    std::ifstream mlist(clist);
    uint_t l = 0;
    while (!mlist.eof()) {
        std::string line;
        std::getline(mlist, line);
        ++l;

        if (line.find_first_not_of(" \t") == line.npos) continue;
        line = trim(line);
        if (line[0] == '#') continue;
        vec1s spl = trim(split(line, "="));
        if (spl.size() != 2) {
            error("wrong format for l."+strn(l)+": '"+line+"'");
            note("expected: map_code_name = map_file, cutout_size");
            return 1;
        }

        vec1s spl2 = trim(split(spl[1], ","));
        if (spl2.size() != 2) {
            error("wrong format for l."+strn(l)+": '"+line+"'");
            note("expected: map_code_name = map_file, cutout_size");
            return 1;
        }

        spl[1] = dir+spl2[0];
        append(spl, spl2[rx(1,spl2.size()-1)]);

        if (!file::exists(spl[1])) {
            error("could not find '"+spl[1]+"'");
            return 1;
        }

        mname.push_back(spl[0]);
        mfile.push_back(spl[1]);

        uint_t hsize;
        if (!from_string(spl[2], hsize)) {
            error("could not convert '"+spl[2]+"' to a cutout size (in pixels)");
            return 1;
        }

        msize.push_back(hsize);
    }

    if (verbose) print("map list loaded successfully");

    for (uint_t b : range(mname)) {
        if (verbose) print(mname[b]);

        fits::header hdr = fits::read_header(mfile[b]);
        double rx = 0, ry = 0;
        fits::getkey(hdr, "CRPIX1", rx);
        fits::getkey(hdr, "CRPIX2", ry);

        vec1d x, y;
        fits::ad2xy(fits::wcs(hdr), ra, dec, x, y);
        x = round(x);
        y = round(y);
        vec1d crx = rx - x + msize[b] + 1;
        vec1d cry = ry - y + msize[b] + 1;

        qstack_params p;
        p.keepnan = true;
        vec3d cube;
        vec1u ids;

        qstack(ra, dec, mfile[b], msize[b], cube, ids, p);

        for (uint_t i : range(ids)) {
            fits::header nhdr = hdr;
            fits::setkey(nhdr, "CRPIX1", crx[ids[i]]);
            fits::setkey(nhdr, "CRPIX2", cry[ids[i]]);
            fits::write(out+name[ids[i]]+mname[b]+".fits", cube(i,_,_), nhdr);
        }

        vec1b ncov(ra.size());
        ncov[_] = true;
        ncov[ids] = false;

        ids = where(ncov);
        vec2d empty(2*msize[b] + 1, 2*msize[b] + 1);
        empty[_] = dnan;
        for (uint_t i : range(ids)) {
            fits::header nhdr = hdr;
            fits::setkey(nhdr, "CRPIX1", crx[ids[i]]);
            fits::setkey(nhdr, "CRPIX2", cry[ids[i]]);
            fits::write(out+name[ids[i]]+mname[b]+".fits", empty, nhdr);
        }
    }

    if (show.size() > 0) {
        if (ra.size() != 1) {
            warning("can only display sources with DS9 for single objects, not catalogs");
        } else {
            if (show.size() > 3) {
                warning("cannot display more than 3 images at the same time, ignoring");
                show.resize(3);
            }

            vec1s chanels = {"red", "green", "blue"};
            chanels = chanels[uindgen(show.size())];

            spawn("ds9 -rgb "+
                collapse("-"+chanels+" "+out+name[0]+show+".fits ")
            );
        }
    }

    return 0;
}

void print_help() {
    print("WIP");
}
