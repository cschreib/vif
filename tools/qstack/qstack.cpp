//  Copyright (c) 2013 C. Schreiber
//
//  This software is provided 'as-is', without any express or implied warranty. In no event will the
//  authors be held liable for any damages arising from the use of this software.
//
//  Permission is granted to anyone to use this software for any purpose, including commercial
//  applications, and to alter it and redistribute it freely, subject to the following restrictions:
//
//      1. The origin of this software must not be misrepresented; you must not claim that you wrote the
//         original software. If you use this software in a product, an acknowledgment in the product
//         documentation would be appreciated but is not required.
//
//      2. Altered source versions must be plainly marked as such, and must not be misrepresented as
//         being the original software.
//
//      3. This notice may not be removed or altered from any source distribution.
//

#include <CCfits/CCfits>
#include <wcslib/wcshdr.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

std::string strtrim(const std::string& str) {
    if (str.empty()) return str;

    std::size_t is = 0;
    while (is != str.size() && (str[is] == ' ' || str[is] == '\t')) ++is;
    std::size_t ie = str.size()-1;
    while (ie != 0          && (str[ie] == ' ' || str[ie] == '\t')) --ie;
    return str.substr(is, ie-is+1);
}

template<typename T>
T from_string(const std::string& str) {
    std::istringstream ss(str);
    T t;
    ss >> t;
    return t;
}

template<typename T>
bool from_string(T& t, const std::string& str) {
    std::istringstream ss(str);
    return ss >> t;
}

template<typename T>
std::string to_string(const T& t) {
    std::ostringstream ss;
    ss << t;
    return ss.str();
}

template<typename T>
std::string to_string(const T& t, std::size_t n) {
    if (t == 0) {
        return std::string(floor(log10(n-1))+1, '0');
    } else {
        std::ostringstream ss;
        ss << t;
        std::size_t nz = floor(log10(n-1)) - floor(log10(t));
        if (nz > 0 && nz < 6) {
            return std::string(nz, '0') + ss.str();
        } else {
            return ss.str();
        }
    }
}

using arguments_t = std::map<std::string, std::vector<std::string>>;
arguments_t read_arguments(int argc, char** argv) {
    arguments_t map;
    for (int i = 1; i < argc; ++i) {
        std::string s = strtrim(argv[i]);
        if (s.empty()) continue;
        if (s[0] != '-') {
            map[""].push_back(s);
            continue;
        }

        s.erase(0,1);
        std::size_t eqpos = s.find_first_of('=');
        if (eqpos != s.npos) {
            map[strtrim(s.substr(0, eqpos))].push_back(strtrim(s.substr(eqpos+1)));
        } else {
            map[s].push_back("");
        }
    }

    map[""];

    return map;
}

template<class T>
bool set_value(T& arg, const std::string& name, const arguments_t& args) {
    auto iter = args.find(name);
    if (iter != args.end()) {
        return from_string(arg, iter->second.back());
    }

    return false;
}

bool set_value(bool& arg, const std::string& name, const arguments_t& args) {
    auto iter = args.find(name);
    if (iter != args.end()) {
        if (iter->second.empty() || iter->second.back().empty()) {
            arg = true;
            return true;
        } else {
            return from_string(arg, iter->second.back());
        }
    }

    return false;
}

bool argument_is_provided(const std::string& name, const arguments_t& args) {
    return args.find(name) != args.end();
}

bool argument_is_provided(const arguments_t& args) {
    return args.find("") != args.end();
}

const std::vector<std::string>& argument_list(const arguments_t& args) {
    return args.find("")->second;
}

namespace fits {
    using namespace CCfits;
}

struct vec2d {
    double x, y;
};

using vec2l = std::vector<long>;

std::istream& operator >> (std::istream& s, vec2d& v) {
    return s >> v.x >> v.y;
}

std::ostream& operator << (std::ostream& s, const vec2d& v) {
    return s << v.x << ", " << v.y;
}

struct vec2d_id {
    std::string id;
    vec2d v;
};

std::istream& operator >> (std::istream& s, vec2d_id& v) {
    return s >> v.id >> v.v;
}

const std::vector<long> strides = {1,1};

void subimage(fits::PHDU& hdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& data) {
    data.resize((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    hdu.read(data, p0, p1, strides);
}

void subimage_add(fits::PHDU& hdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& data) {
    std::valarray<double> slice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    hdu.read(slice, p0, p1, strides);
    data += slice;
}

void subimage_add_inv(fits::PHDU& hdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& data) {
    std::valarray<double> slice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    hdu.read(slice, p0, p1, strides);
    data += 1.0/slice;
}

void subimage_add_invsqr(fits::PHDU& hdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& data) {
    std::valarray<double> slice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    hdu.read(slice, p0, p1, strides);
    slice *= slice;
    data += 1.0/slice;
}

void subimage_add_sqr(fits::PHDU& hdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& data) {
    std::valarray<double> slice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    hdu.read(slice, p0, p1, strides);
    slice *= slice;
    data += slice;
}

void subimage_add_weight(fits::PHDU& vhdu, fits::PHDU& whdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& vdata, std::valarray<double>& wdata) {
    std::valarray<double> vslice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    std::valarray<double> wslice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    vhdu.read(vslice, p0, p1, strides);
    whdu.read(wslice, p0, p1, strides);
    vslice *= wslice;
    vdata += vslice;
    wdata += wslice;
}

void subimage_add_weightinv(fits::PHDU& vhdu, fits::PHDU& whdu, const vec2l& p0, const vec2l& p1, std::valarray<double>& vdata, std::valarray<double>& wdata) {
    std::valarray<double> vslice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    std::valarray<double> wslice((p1[0]-p0[0]+1)*(p1[1]-p0[1]+1));
    vhdu.read(vslice, p0, p1, strides);
    whdu.read(wslice, p0, p1, strides);
    wslice = 1.0/wslice;
    vslice *= wslice;
    vdata += vslice;
    wdata += wslice;
}

void print_help() {
    std::cout << "qstack v1.1\n";
    std::cout << "  usage: qstack [-options=...] config.cfg\n\n";

    std::cout << "  The content of 'config.cfg' must be as follow:\n" ;
    std::cout << "    1) a list of stacking parameters (see below), each line starting with a '#'\n";
    std::cout << "       followed by the name of the parameter, a semicolon ':', and the value.\n";
    std::cout << "       It is possible to write comments by starting the line with '##'.\n";
    std::cout << "    2) the list of all the sources to stack, each line containing the right\n";
    std::cout << "       ascension followed by the declination of the source (in degrees).\n";
    std::cout << "       Optionally, this list can be complemented by an ID collumn (must be\n";
    std::cout << "       the first collumn), which can be any string that does not contain any\n";
    std::cout << "       space. This ID will be used to name the individual cutouts.\n\n";


    std::cout << "  List of available command line options:\n";
    std::cout << "   -verbose: set this flag to print some information in the standard output\n";
    std::cout << "   -cov=[number]: set this value to the minimum allowed coverage for a source\n";
    std::cout << "                  to be stacked. When set, no stacking is performed: a new\n";
    std::cout << "                  source list is sent to the standard output.\n";
    std::cout << "   -cutouts: set this flag to simply output the cutouts. No stacking is performed.\n";

    std::cout << "  List of available parameters:\n";
    std::cout << "   - scifile: path to the file containing the image data (string, mandatory)\n";
    std::cout << "   - errfile: path to the file containing the error data (string, default: none)\n";
    std::cout << "   - whtfile: path to the file containing the weight data (string, default: none)\n";
    std::cout << "   - covfile: path to the file containing the coverage data (string, default: none)\n";
    std::cout << "   - output:  path to the base output file name (string, default: stout)\n";
    std::cout << "   - hsize:   half size of the resulting stack (unsigned integer, mandatory)\n";
    std::cout << "   - median:  do median stacking (1) or mean stacking (0) (boolean, default: 0)\n";
    std::cout << "   - noerror: compute the error on the stack (0) or not (1) (boolean, default: 0)\n\n";

    std::cout << "  It is only necessary to provide either 'errfile' or 'whtfile' in order to\n";
    std::cout << "  compute the error on the stack, the weight data being interpreted as the\n";
    std::cout << "  inverse of the error.\n\n";

    std::cout << "  When doing mean stacking (default), each source is weighted by the inverse\n";
    std::cout << "  of the error on the map. They are then all summed together and the result\n";
    std::cout << "  is divided by the sum of the weights. The returned error is computed as the\n";
    std::cout << "  square root of the mean of the squared errors among all sources.\n\n";

    std::cout << "  When doing median stacking (median:1), all the sources are kept into a data\n";
    std::cout << "  cube. The median is then computed among all sources for each pixel of the\n";
    std::cout << "  final stack. The error is computed the same way as for mean stacking.\n\n";

    std::cout << "  The output of the program is:\n";
    std::cout << "   - <output>_sci.fits: the stacked image.\n";
    std::cout << "   - <output>_err.fits: the stacked error.\n\n";

    std::cout << "  Copyright (c) 2013 C. Schreiber (corentin.schreiber@cea.fr)\n\n";
    std::cout << "  This software is provided 'as-is', without any express or implied warranty.\n";
    std::cout << "  In no event will the authors be held liable for any damages arising from the\n";
    std::cout << "  use of this software.\n";

    std::cout << "  Permission is granted to anyone to use this software for any purpose,\n";
    std::cout << "  including commercial applications, and to alter it and redistribute it\n";
    std::cout << "  freely, subject to the following restrictions:\n\n";

    std::cout << "    1. The origin of this software must not be misrepresented; you must not\n";
    std::cout << "       claim that you wrote the original software. If you use this software in\n";
    std::cout << "       a product, an acknowledgment in the product documentation would be\n";
    std::cout << "       appreciated but is not required.\n\n";

    std::cout << "    2. Altered source versions must be plainly marked as such, and must not be\n";
    std::cout << "       misrepresented as being the original software.\n\n";

    std::cout << "    3. This notice may not be removed or altered from any source distribution.\n\n";

    std::cout << std::flush;
}

class wcs_astro {
    wcsprm* wcs = nullptr;
    int nwcs    = 0;

public:
    wcs_astro(fits::PHDU& hdu) {
        fitsfile* fp = hdu.fitsPointer();

        // Read the header as a string
        char* hstr = nullptr;
        int nkeys  = 0;
        int status = 0;
        fits_hdr2str(fp, 0, nullptr, 0, &hstr, &nkeys, &status);

        // Feed it to WCSLib to extract the astrometric parameters
        int nreject = 0;
        wcspih(hstr, nkeys, WCSHDR_all, 0, &nreject, &nwcs, &wcs);
        free(hstr);
    }

    void convert(const std::vector<vec2d>& world, std::vector<vec2d>& pos) const {
        std::size_t ngal = world.size();

        std::vector<double> phi(ngal), theta(ngal);
        std::vector<vec2d>  itmp(ngal);
        std::vector<int>    stat(ngal);
        pos.resize(ngal);

        wcss2p(wcs, ngal, 2, reinterpret_cast<const double*>(world.data()), phi.data(), theta.data(),
            reinterpret_cast<double*>(itmp.data()), reinterpret_cast<double*>(pos.data()), stat.data());
    }

    ~wcs_astro() {
        wcsvfree(&nwcs, &wcs);
        wcs = nullptr;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    // Mandatory parameters
    std::string scifile = "";
    std::size_t hsize   = 0;
    std::string cfgfile = "";

    // Optional parameters
    std::string errfile;
    std::string whtfile;
    std::string covfile;
    std::string output  = "stout";
    bool        median  = false;
    bool        noerror = false;
    bool        errwht  = false;
    bool        verbose = false;
    bool        cutouts = false;
    double      mincov  = -1.0;

    // Read configuration from the command line arguments
    auto args = read_arguments(argc, argv);

    // Get the configuration file
    auto arglist = argument_list(args);
    if (arglist.empty()) {
        print_help();
        return 1;
    } else {
        cfgfile = arglist[0];
    }

    // Fetch optional keywords and parameters
    set_value(verbose, "verbose", args);
    set_value(mincov,  "cov", args);
    set_value(cutouts, "cutouts", args);

    // Read the configuration file
    std::vector<vec2d> world;
    std::vector<std::string> oids;
    bool list_found = false;

    std::ifstream cfg(cfgfile);
    std::istream& in = cfg;
    std::string line;
    while (!list_found) {
        std::getline(in, line);

        if (in.eof()) {
            std::cout << "error: ill formed configuration file: missing the source list" << std::endl;
            return 1;
        }

        if (strtrim(line).empty()) continue;

        if (line[0] == '#') {
            line = line.substr(1, line.size()-1);
            if (!line.empty() && line[0] == '#') continue;
            std::size_t seppos = line.find_first_of(':');
            if (seppos == line.npos) continue;

            std::string key = strtrim(line.substr(0, seppos));
            std::string val = strtrim(line.substr(seppos+1, line.size()-seppos-1));

            if        (key == "scifile") {
                scifile = val;
            } else if (key == "errfile") {
                errfile = val;
            } else if (key == "whtfile") {
                whtfile = val;
            } else if (key == "covfile") {
                covfile = val;
            } else if (key == "output") {
                output = val;
            } else if (key == "hsize") {
                from_string(hsize, val);
            } else if (key == "median") {
                from_string(median, val);
            } else if (key == "noerror") {
                from_string(noerror, val);
            } else {
                std::cout << "warning: unknown parameter '" << key << "'" << std::endl;
            }
        } else {
            list_found = true;
            std::size_t cols = 0;
            std::size_t pos = 0;
            while ((pos = line.find_first_not_of(" \t", pos)) != line.npos) {
                ++cols;
                pos = line.find_first_of(" \t", pos);
                if (pos == line.npos) break;
            }

            if (cols < 2) {
                std::cout << "error: source list must be composed of at least 2 collumns: RA and Dec coordinates" << std::endl;
            }

            if (cols == 2) {
                world.push_back(from_string<vec2d>(line));

                while (std::getline(in, line) and !in.eof()) {
                    world.push_back(from_string<vec2d>(line));
                }
            } else {
                vec2d_id v;
                from_string(v, line);
                oids.push_back(v.id);
                world.push_back(v.v);

                while (std::getline(in, line) and !in.eof()) {
                    from_string(v, line);
                    oids.push_back(v.id);
                    world.push_back(v.v);
                }
            }
        }
    }

    std::size_t ngal = world.size();

    try {
        if (mincov >= 0.0) {

            // -- Coverage test --

            if (covfile.empty()) {
                std::cout << "error: no coverage map is provided ('covfile'), cannot build the covered source list" << std::endl;
                return 1;
            }

            // Open the file
            fits::FITS cov(covfile, fits::Read);

            // Get the first HDU (extension)
            fits::PHDU& covhdu = cov.pHDU();

            if (covhdu.axes() != 2) {
                std::cout << "error: no image data (dimension: " << covhdu.axes() << ")" << std::endl;
                return 1;
            }

            long sx = covhdu.axis(0);
            long sy = covhdu.axis(1);

            // Convert world coordinates to pixel coordinates
            std::vector<vec2d> pos;
            wcs_astro(covhdu).convert(world, pos);

            // Check that all the sources are properly covered
            std::valarray<double> subcov(0.0, (2*hsize+1)*(2*hsize+1));
            std::vector<std::size_t> ids(ngal);
            for (std::size_t i = 0; i < ngal; ++i) ids[i] = i;
            auto iter = ids.begin();

            while (iter != ids.end()) {
                std::size_t i = *iter;
                if (pos[i].x - hsize < 1 || pos[i].x + hsize > sx ||
                    pos[i].y - hsize < 1 || pos[i].y + hsize > sy) {
                    iter = ids.erase(iter);
                } else {
                    vec2l p0 = {long(round(pos[i].x - hsize)), long(round(pos[i].y - hsize))};
                    vec2l p1 = {long(round(pos[i].x + hsize)), long(round(pos[i].y + hsize))};
                    subimage(covhdu, p0, p1, subcov);
                    double tcov = subcov[hsize + hsize*(2*hsize+1)]*4 + subcov[0] + subcov[2*hsize+1] + subcov[(2*hsize+1)*2*hsize] + subcov[(2*hsize+1)*(2*hsize+1)];
                    if (tcov >= 8.0*mincov) {
                        ++iter;
                    } else {
                        iter = ids.erase(iter);
                    }
                }
            }

            // Print the list of covered sources
            std::cout << "## " << ids.size() << " covered sources out of " << ngal << std::endl;
            for (auto i : ids) {
                std::cout << world[i].x << '\t' << world[i].y << '\n';
            }

            std::cout << std::flush;

        } else {

            // -- Actual source extraction --

            if (errfile.empty() && whtfile.empty()) {
                if (!noerror) {
                    std::cout << "warning: no error or weight map provided ('errfile', 'whtfile'), cannot compute the error on the stack" << std::endl;
                }
                noerror = true;
            } else {
                if (errfile.empty()) {
                    errfile = whtfile;
                    errwht = true;
                }
            }

            // Open the file
            fits::FITS sci(scifile, fits::Read);

            // Get the first HDU (extension)
            fits::PHDU& scihdu = sci.pHDU();

            if (scihdu.axes() != 2) {
                std::cout << "error: no image data (dimension: " << scihdu.axes() << ")" << std::endl;
                return 1;
            }

            long sx = scihdu.axis(0);
            long sy = scihdu.axis(1);

            // Convert the world positions to pixels positions
            std::vector<vec2d> pos;
            wcs_astro(scihdu).convert(world, pos);

            // Check that all the sources are properly covered
            std::vector<std::size_t> ids(ngal);
            for (std::size_t i = 0; i < ngal; ++i) ids[i] = i;

            auto iter = ids.begin();
            while (iter != ids.end()) {
                std::size_t i = *iter;
                if (pos[i].x - hsize < 1 || pos[i].x + hsize > sx ||
                    pos[i].y - hsize < 1 || pos[i].y + hsize > sy) {
                    iter = ids.erase(iter);
                } else {
                    ++iter;
                }
            }

            std::size_t ongal = ngal;
            ngal = ids.size();
            if (verbose) std::cout << "stacked sources: " << ngal << " out of " << ongal << std::endl;

            std::vector<vec2l> p0(ngal);
            std::vector<vec2l> p1(ngal);
            for (std::size_t i = 0; i < ngal; ++i) {
                p0[i] = {long(round(pos[ids[i]].x - hsize)), long(round(pos[ids[i]].y - hsize))};
                p1[i] = {long(round(pos[ids[i]].x + hsize)), long(round(pos[ids[i]].y + hsize))};
            }

            // Stack the images
            std::valarray<double> vstack(0.0, (2*hsize+1)*(2*hsize+1));
            std::valarray<double> estack(0.0, (2*hsize+1)*(2*hsize+1));

            if (noerror) {
                if (cutouts) {
                    // -- Output cutouts on the image only --
                    long isize[2]; isize[0] = isize[1] = 2*hsize+1;
                    std::valarray<double> vsub(0.0, (2*hsize+1)*(2*hsize+1));
                    for (std::size_t i = 0; i < ngal; ++i) {
                        subimage(scihdu, p0[i], p1[i], vsub);
                        std::string fname;
                        if (oids.empty()) {
                            fname = "!"+output+to_string(ids[i], ongal)+"_sci.fits";
                        } else {
                            fname = "!"+output+oids[ids[i]]+"_sci.fits";
                        }
                        fits::FITS vout(fname, DOUBLE_IMG, 2, isize);
                        vout.pHDU().write(1, vsub.size(), vsub);
                        vout.flush();
                    }
                } else {
                    if (median) {
                        // -- Median stacking on the image only --
                        // First build a data cube containing all the cutouts
                        std::vector<std::valarray<double>> vcube(ngal, vstack);
                        for (std::size_t i = 0; i < ngal; ++i) {
                            subimage_add(scihdu, p0[i], p1[i], vcube[i]);
                        }

                        // Then compute the median of each pixel
                        std::vector<double> vcollumn(ngal);
                        for (std::size_t y = 0; y < 2*hsize+1; ++y)
                        for (std::size_t x = 0; x < 2*hsize+1; ++x) {
                            for (std::size_t i = 0; i < ngal; ++i) {
                                vcollumn[i] = vcube[i][x+y*(2*hsize+1)];
                            }

                            std::nth_element(vcollumn.begin(), vcollumn.begin()+ngal/2, vcollumn.end());

                            vstack[x+y*(2*hsize+1)] = *(vcollumn.begin()+ngal/2);
                        }
                    } else {
                        // -- Mean stacking on the image only --
                        // Add up all the values of each cutout to the stack and divide by the number of images
                        std::valarray<double> wstack(0.0, (2*hsize+1)*(2*hsize+1));
                        for (std::size_t i = 0; i < ngal; ++i) {
                            subimage_add(scihdu, p0[i], p1[i], vstack);
                        }

                        vstack /= ngal;
                    }
                }
            } else {
                fits::FITS err(errfile, fits::Read);
                fits::PHDU& errhdu = err.pHDU();
                long sex = errhdu.axis(0);
                long sey = errhdu.axis(1);

                if (sex != sx || sey != sy) {
                    std::cout << "error: sci and " << std::string(errwht ? "wht" : "err") <<  " FITS files do not match" << std::endl;
                    return 1;
                }

                if (cutouts) {
                    // -- Output cutouts on the image and the error map --
                    std::valarray<double> vsub(0.0, (2*hsize+1)*(2*hsize+1));
                    long isize[2]; isize[0] = isize[1] = 2*hsize+1;
                    for (std::size_t i = 0; i < ngal; ++i) {
                        subimage(scihdu, p0[i], p1[i], vsub);
                        std::string fname;
                        if (oids.empty()) {
                            fname = "!"+output+to_string(ids[i], ongal)+"_sci.fits";
                        } else {
                            fname = "!"+output+oids[ids[i]]+"_sci.fits";
                        }
                        fits::FITS vout(fname, DOUBLE_IMG, 2, isize);
                        vout.pHDU().write(1, vsub.size(), vsub);
                        vout.flush();
                    }

                    std::valarray<double> esub(0.0, (2*hsize+1)*(2*hsize+1));
                    for (std::size_t i = 0; i < ngal; ++i) {
                        subimage(errhdu, p0[i], p1[i], esub);
                        if (errwht) esub = 1.0/esub;
                        std::string fname;
                        if (oids.empty()) {
                            fname = "!"+output+to_string(ids[i], ongal)+"_err.fits";
                        } else {
                            fname = "!"+output+oids[ids[i]]+"_err.fits";
                        }
                        fits::FITS eout(fname, DOUBLE_IMG, 2, isize);
                        eout.pHDU().write(1, esub.size(), esub);
                        eout.flush();
                    }
                } else {
                    if (errwht) {
                        if (median) {
                            // -- Median stacking on the image, mean on the weight map --
                            // First build a data cube containing all the cutouts
                            std::vector<std::valarray<double>> vcube(ngal, vstack);
                            for (std::size_t i = 0; i < ngal; ++i) {
                                subimage_add_invsqr(errhdu, p0[i], p1[i], estack);
                                subimage_add(scihdu, p0[i], p1[i], vcube[i]);
                            }

                            // Then compute the median of each pixel
                            std::vector<double> vcollumn(ngal);
                            std::vector<double> ecollumn(ngal);
                            for (std::size_t y = 0; y < 2*hsize+1; ++y)
                            for (std::size_t x = 0; x < 2*hsize+1; ++x) {
                                for (std::size_t i = 0; i < ngal; ++i) {
                                    vcollumn[i] = vcube[i][x+y*(2*hsize+1)];
                                }

                                std::nth_element(vcollumn.begin(), vcollumn.begin()+ngal/2, vcollumn.end());

                                vstack[x+y*(2*hsize+1)] = *(vcollumn.begin()+ngal/2);
                            }

                            estack /= ngal;
                            estack = sqrt(estack);
                        } else {
                            // -- Mean stacking on the image and the weight map --
                            // Add up all the values of each cutout to the stack, weighted by the error map
                            std::valarray<double> wstack(0.0, (2*hsize+1)*(2*hsize+1));
                            for (std::size_t i = 0; i < ngal; ++i) {
                                subimage_add_invsqr(errhdu, p0[i], p1[i], estack);
                                subimage_add_weight(scihdu, errhdu, p0[i], p1[i], vstack, wstack);
                            }

                            estack /= ngal;
                            estack = sqrt(estack);
                            vstack /= wstack;
                        }
                    } else {
                        if (median) {
                            // -- Median stacking on the image, mean on the error map --
                            // First build a data cube containing all the cutouts
                            std::vector<std::valarray<double>> vcube(ngal, vstack);
                            for (std::size_t i = 0; i < ngal; ++i) {
                                subimage_add_sqr(errhdu, p0[i], p1[i], estack);
                                subimage_add(scihdu, p0[i], p1[i], vcube[i]);
                            }

                            // Then compute the median of each pixel
                            std::vector<double> vcollumn(ngal);
                            for (std::size_t y = 0; y < 2*hsize+1; ++y)
                            for (std::size_t x = 0; x < 2*hsize+1; ++x) {
                                for (std::size_t i = 0; i < ngal; ++i) {
                                    vcollumn[i] = vcube[i][x+y*(2*hsize+1)];
                                }

                                std::nth_element(vcollumn.begin(), vcollumn.begin()+ngal/2, vcollumn.end());

                                vstack[x+y*(2*hsize+1)] = *(vcollumn.begin()+ngal/2);
                            }

                            estack /= ngal;
                            estack = sqrt(estack);
                        } else {
                            // -- Mean stacking on the image and error map --
                            // Add up all the values of each cutout to the stack, weighted by the error map
                            std::valarray<double> wstack(0.0, (2*hsize+1)*(2*hsize+1));
                            for (std::size_t i = 0; i < ngal; ++i) {
                                subimage_add_sqr(errhdu, p0[i], p1[i], estack);
                                subimage_add_weightinv(scihdu, errhdu, p0[i], p1[i], vstack, wstack);
                            }

                            estack /= ngal;
                            estack = sqrt(estack);
                            vstack /= wstack;
                        }
                    }
                }
            }

            if (!cutouts) {
                long isize[2]; isize[0] = isize[1] = 2*hsize+1;
                fits::FITS vout("!"+output+"_sci.fits", DOUBLE_IMG, 2, isize);
                vout.pHDU().write(1, vstack.size(), vstack);
                vout.flush();

                if (!noerror) {
                    fits::FITS eout("!"+output+"_err.fits", DOUBLE_IMG, 2, isize);
                    eout.pHDU().write(1, estack.size(), estack);
                    eout.flush();
                }
            }
        }
    } catch (fits::FitsException& e) {
        std::cout << "qstack: working on " << cfgfile << "\n";
        std::cout << "error: " << e.message() << std::endl;
        return 1;
    }

    return 0;
}

