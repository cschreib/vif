#include <phypp.hpp>

int main(int argc, char* argv[]) {
    struct map_info {
        std::string band;
        std::string img;
        std::string err;
        std::string psf;
        float fconv = 1.0;
        float group_fit_threshold = 1.0;
        float group_aper_threshold = 1.0;
        float beam_smear = fnan;
    };

    std::vector<map_info> maps;

    std::ifstream imaps(argv[3]);
    std::string line;
    map_info current;
    while (std::getline(imaps, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vec1s spl = trim(split(line, "="));
        if (spl.size() == 1) {
            if (!current.band.empty()) {
                if (current.img.empty()) {
                    error("missing 'img' parameter for band '", current.band, "'");
                    return 1;
                } else if (current.err.empty()) {
                    error("missing 'img' parameter for band '", current.band, "'");
                    return 1;
                } else if (current.psf.empty()) {
                    error("missing 'img' parameter for band '", current.band, "'");
                    return 1;
                }

                maps.push_back(current);
            }

            current = map_info{};
            current.band = line;
        } else if (spl.size() == 2) {
            if (spl[0] == "sci") {
                current.img = spl[1];
            } else if (spl[0] == "err") {
                current.err = spl[1];
            } else if (spl[0] == "psf") {
                current.psf = spl[1];
            } else if (spl[0] == "fconv") {
                if (!from_string(spl[1], current.fconv)) {
                    error("could not read flux conversion factor from '", spl[1], "'");
                    return 1;
                }
            } else if (spl[0] == "smooth_radius") {
                if (!from_string(spl[1], current.beam_smear)) {
                    error("could not read beam smoothing radius from '", spl[1], "'");
                    return 1;
                }
            } else if (spl[0] == "gfit_threshold") {
                if (!from_string(spl[1], current.group_fit_threshold)) {
                    error("could not read group fit threshold from '", spl[1], "'");
                    return 1;
                }
            } else if (spl[0] == "gaper_threshold") {
                if (!from_string(spl[1], current.group_aper_threshold)) {
                    error("could not read group aperture threshold radius from '", spl[1], "'");
                    return 1;
                }
            }
        } else {
            error("ill formed line: ", line);
            return 1;
        }
    }

    if (!current.band.empty()) {
        if (current.img.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return 1;
        } else if (current.err.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return 1;
        } else if (current.psf.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return 1;
        }

        maps.push_back(current);
    }

    std::string cat_file = argv[4];
    std::string out_file = argv[5];

    vec1d ra, dec;
    fits::read_table(cat_file, ftable(ra, dec));
    fits::update_table(out_file, ftable(ra, dec));

    // Write down residual code
    std::ofstream oscript(argv[2]);
    oscript << "#!/bin/bash\n\n";
    for (uint_t i : range(maps)) {
        oscript << "subsrc \""+maps[i].img+"\" \""+out_file+"\" "
                << "bands="+maps[i].band+" psf=\"file,"+maps[i].psf+"\" "
                << "fconv="+strn(maps[i].fconv)+" "
                << "out="+maps[i].band+"-fullres.fits";

        oscript << "\n";
    }

    return 0;
}
