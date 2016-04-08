#include "pixfit-common.hpp"

int main(int argc, char* argv[]) {
    std::vector<map_info> maps;
    if (!read_maps(argv[2], maps)) return 1;

    std::string cat_file = argv[3];
    std::string cmb_file = argv[4];

    std::string catalogs = "\""+cmb_file+"\" \""+cat_file+"\"";
    std::ofstream oscript(argv[1]);
    oscript << "#!/bin/bash\n\n";
    for (uint_t i : range(maps)) {
        oscript << "pixfit-extract cat=\""+cat_file+"\" img=\""+maps[i].band+"-sci.fits\" "
                << "err=\""+maps[i].band+"-err.fits\" psf=\""+maps[i].psf+"\" "
                << "fconv="+strn(maps[i].fconv)+" "
                << "beam_flux="+strn(maps[i].beam_flux)+" "
                << "out=fit-"+maps[i].band+".fits "
                << "out_res="+maps[i].band+"-res.fits "
                << "out_model="+maps[i].band+"-mod.fits make_groups "
                << "out_gmap="+maps[i].band+"-grp.fits "
                << "group_fit_threshold="+strn(maps[i].group_fit_threshold)+" "
                << "group_aper_threshold="+strn(maps[i].group_aper_threshold)+" ";

        if (is_finite(maps[i].beam_smear)) {
            oscript << " beam_smeared beam_size=" << strn(maps[i].beam_smear);
        }

        oscript << "\n";

        catalogs += " "+maps[i].band+":fit-"+maps[i].band+".fits";
    }

    oscript << "\n";
    oscript << "pixfit-combine "+catalogs+"\n";

    return 0;
}
