#include <phypp.hpp>

struct map_info {
    std::string band;
    std::string img;
    std::string err;
    std::string psf;
    float fconv = 1.0;
    float group_fit_threshold = 0.97;
    float group_aper_threshold = 0.8;
    float group_aper_size = fnan;
    float beam_flux = 1.0;
    float beam_smear = fnan;
};

bool read_maps(std::string filename, std::vector<map_info>& maps) {
    std::ifstream imaps(filename);

    if (!imaps.is_open()) {
        error("could not find map file '", filename, "'");
        return false;
    }

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
                    return false;
                } else if (current.err.empty()) {
                    error("missing 'img' parameter for band '", current.band, "'");
                    return false;
                } else if (current.psf.empty()) {
                    error("missing 'img' parameter for band '", current.band, "'");
                    return false;
                }

                maps.push_back(current);
            }

            current = map_info{};
            current.band = line;
        } else if (spl.size() == 2) {
            auto fix_path = [&](std::string path) {
                if (!path.empty() && path[0] != '/') {
                    path = file::get_directory(filename)+path;
                }

                return path;
            };

            if (spl[0] == "sci") {
                current.img = fix_path(spl[1]);
            } else if (spl[0] == "err") {
                current.err = fix_path(spl[1]);
            } else if (spl[0] == "psf") {
                current.psf = fix_path(spl[1]);
            } else if (spl[0] == "fconv") {
                if (!from_string(spl[1], current.fconv)) {
                    error("could not read flux conversion factor from '", spl[1], "'");
                    return false;
                }
            } else if (spl[0] == "smooth_radius") {
                if (!from_string(spl[1], current.beam_smear)) {
                    error("could not read beam smoothing radius from '", spl[1], "'");
                    return false;
                }
            } else if (spl[0] == "beam_flux") {
                if (!from_string(spl[1], current.beam_flux)) {
                    error("could not read beam flux from '", spl[1], "'");
                    return false;
                }
            } else if (spl[0] == "gfit_threshold") {
                if (!from_string(spl[1], current.group_fit_threshold)) {
                    error("could not read group fit threshold from '", spl[1], "'");
                    return false;
                }
            } else if (spl[0] == "gaper_threshold") {
                if (!from_string(spl[1], current.group_aper_threshold)) {
                    error("could not read group aperture threshold radius from '", spl[1], "'");
                    return false;
                }
            } else if (spl[0] == "gaper_size") {
                if (!from_string(spl[1], current.group_aper_size)) {
                    error("could not read group flux aperture size from '", spl[1], "'");
                    return false;
                }
            }
        } else {
            error("ill formed line: ", line);
            return false;
        }
    }

    if (!current.band.empty()) {
        if (current.img.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return false;
        } else if (current.err.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return false;
        } else if (current.psf.empty()) {
            error("missing 'img' parameter for band '", current.band, "'");
            return false;
        }

        maps.push_back(current);
    }

    return true;
}

bool read_ds9_region_circles(std::string file_name, vec2d& regs, vec1s& text, bool& physical, std::string color = "") {
    std::ifstream file(file_name);

    if (!file.is_open()) {
        error("could not find region file '", file_name, "'");
        return false;
    }

    std::string global_color = "green";
    physical = false;
    std::string line;
    uint_t l = 0;
    while (std::getline(file, line)) {
        ++l;
        if (line.empty() || trim(line).empty() || trim(line)[0] == '#') continue;

        if (begins_with(line, "global")) {
            std::string key = "color=";
            auto pos = line.find(key);
            if (pos != line.npos) {
                pos += key.size();
                global_color = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
            }
            continue;
        }

        auto spos = line.find_first_of('(');
        if (spos == line.npos) {
            if (trim(line) == "fk5") physical = false;
            if (trim(line) == "physical") physical = true;
            continue;
        }

        std::string type = trim(line.substr(0, spos));
        if (type != "circle") continue;

        auto epos = line.find_first_of(')', spos+1);
        std::string targs = line.substr(spos+1, epos-(spos+1));
        vec1s args = split(targs, ",");
        if (args.size() != 3) {
            error(file_name, ":", l, ": ",
                "ill formed 'circle' line, expecting 3 arguments, got ", args.size());
            return false;
        }

        double ra, dec, rad;
        args = trim(args);
        if (args[0].find_first_of(':') != args[0].npos) {
            if (!sex2deg(args[0], args[1], ra, dec)) {
                error(file_name, ":", l, ": ",
                    "could not convert sexagesimal coordinates to degrees");
                return false;
            }

            if (!ends_with(args[2], "\"")) {
                error(file_name, ":", l, ": expected radius in arcsec");
                return false;
            }
        } else {
            if (!from_string(args[0], ra) || !from_string(args[1], dec)) {
                error(file_name, ":", l, ": ",
                    "could not read coordinates to ", (physical ? "(x,y)" : "degrees"));
                return false;
            }
        }

        if (physical) {
            if (!from_string(args[2], rad)) {
                error(file_name, ":", l, ": could not read radius in pixels");
                return false;
            }
        } else {
            args[2] = erase_end(args[2], "\"");
            if (!from_string(args[2], rad)) {
                error(file_name, ":", l, ": could not read radius in arcsec");
                return false;
            }
        }

        if (!color.empty()) {
            std::string rcol = global_color;
            spos = line.find_first_of('#', epos);
            if (spos != line.npos) {
                std::string key = "color=";
                auto pos = line.find(key, spos+1);
                if (pos != line.npos) {
                    pos += key.size();
                    rcol = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
                }
            }

            if (rcol != color) continue;
        }

        std::string txt;
        spos = line.find_first_of('#', epos);
        if (spos != line.npos) {
            std::string key = "text=";
            auto pos = line.find(key, spos+1);
            if (pos != line.npos) {
                pos += key.size()+1;
                txt = line.substr(pos, line.find_first_of("}", pos)-pos);
            }
        }

        text.push_back(txt);
        append<0>(regs, vec2d{{ra, dec, rad}});
    }

    return true;
}

vec2d read_psf(const std::string& psf_file, double beam_smear, int_t& hsize) {
    vec2d psf;
    fits::read(psf_file, psf);

    hsize = 0;

    vec1i idm = mult_ids(psf, max_id(psf));
    int_t imax = std::min(
        std::min(idm[0], int_t(psf.dims[0])-1-idm[0]),
        std::min(idm[1], int_t(psf.dims[1])-1-idm[1])
    );

    double psf_max = psf(idm[0], idm[1]);

    for (int_t i = imax; i > 0; --i) {
        if (abs(psf(idm[0]-i,idm[1])/psf_max) > 1e-4 ||
            abs(psf(idm[0]+i,idm[1])/psf_max) > 1e-4 ||
            abs(psf(idm[0],idm[1]-i)/psf_max) > 1e-4 ||
            abs(psf(idm[0],idm[1]+i)/psf_max) > 1e-4) {
            hsize = i;
            break;
        }
    }

    if (hsize == 0) hsize = imax;
    psf = subregion(psf, {idm[0]-hsize, idm[1]-hsize, idm[0]+hsize, idm[1]+hsize});

    if (is_finite(beam_smear)) {
        vec2d kernel;
        if (beam_smear > 0) {
            kernel = gaussian_profile(psf.dims, beam_smear);
        } else {
            kernel = psf;
        }

        vec1u idz = where(kernel < 1e-5*max(abs(kernel)));
        kernel[idz] = 0.0;
        psf = convolve2d(psf, kernel);
        psf /= psf(hsize,hsize);
    }

    return psf;
}

vec2d read_psf(const map_info& map, int_t& hsize) {
    return read_psf(map.psf, map.beam_smear, hsize);
}

// Find the first value in 'x' where 'y' equals zero, linearly interpolating between
// the gridded data points
template<typename TypeX, typename TypeY>
double find_zero(const vec<1,TypeX>& x, const vec<1,TypeY>& y) {
    phypp_check(x.size() == y.size(),
        "incompatible x and y array dimensions (", x.size(), " vs ", y.size(), ")");

    for (uint_t i : range(1, x.size())) {
        if (y.safe[i]*y.safe[i-1] <= 0) {
            return interpolate(x.safe[i-1], x.safe[i], y.safe[i-1], y.safe[i], 0.0);
        }
    }

    return dnan;
};
