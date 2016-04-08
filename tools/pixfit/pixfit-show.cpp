#include "pixfit-common.hpp"

int main(int argc, char* argv[]) {
    std::vector<map_info> maps;
    if (!read_maps(argv[3], maps)) return 1;

    std::string cat_file = argv[4];
    std::string out_file = argv[5];

    vec1d ra, dec;
    fits::read_table(cat_file, ftable(ra, dec));
    fits::update_table(out_file, ftable(ra, dec));

    // Write down residual code
    std::ofstream oscript(argv[2]);
    oscript << "#!/bin/bash\n\n";
    for (uint_t i : range(maps)) {
        oscript << "subsrc \""+maps[i].band+"-sci.fits\" \""+out_file+"\" "
                << "bands="+maps[i].band+" psf=\"file,"+maps[i].psf+"\" "
                << "fconv="+strn(maps[i].fconv)+" "
                << "out="+maps[i].band+"-fullres.fits";

        oscript << "\n";
    }

    return 0;
}
