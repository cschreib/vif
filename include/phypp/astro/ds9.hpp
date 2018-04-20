#ifndef PHYPP_ASTRO_DS9_HPP
#define PHYPP_ASTRO_DS9_HPP

#include "phypp/astro/wcs.hpp"

namespace phypp {
namespace astro {
namespace ds9 {
    struct region {
        std::string type, color, text;
        bool physical = true;

        vec1d params;
        // circle: params = {center_x,center_y,radius}
        // box:    params = {center_x,center_y,width,height,angle}
    };

    inline void read_regions(const std::string& filename, vec<1,region>& regs) {
        phypp_check(file::exists(filename), "could not open region file '"+filename+"'");

        std::ifstream file(filename);

        std::string global_color = "green";
        bool physical = false;
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
                if (trim(line) == "physical") physical = true;
                continue;
            }

            // New region
            region r;
            r.type = trim(line.substr(0, spos));
            r.color = global_color;
            r.physical = physical;

            // Isolate region parameters
            auto epos = line.find_first_of(')', spos+1);
            std::string targs = line.substr(spos+1, epos-(spos+1));
            vec1s args = trim(split(targs, ","));

            if (args.size() < 2) {
                warning("ill formed line, expecting at least 2 arguments, got ", args.size());
                warning("at ", filename, ":", l);
                continue;
            }

            // Read region center
            double x0, y0;
            if (args[0].find_first_of(':') != args[0].npos) {
                if (physical) {
                    warning("cannot work with sexagesimal coordinates in 'physical' mode");
                    warning("at ", filename, ":", l);
                    continue;
                }

                if (!sex2deg(args[0], args[1], x0, y0)) {
                    warning("could not convert sexagesimal coordinates to degrees");
                    warning("at ", filename, ":", l);
                    continue;
                }
            } else {
                if (!from_string(args[0], x0) || !from_string(args[1], y0)) {
                    warning("could not read coordinates to "+std::string(physical ? "(x,y)" : "degrees"));
                    warning("at ", filename, ":", l);
                    continue;
                }
            }

            if (physical) {
                // Subtract 1 from X and Y to go from FITS convention to C++ convention
                x0 -= 1.0; y0 -= 1.0;
            }

            // Helper function to read dimensions in pixels or arcsec
            auto read_dim = [&](const std::string& name, std::string from, double& to) {
                if (physical) {
                    if (!from_string(from, to)) {
                        warning("could not read "+name+" in pixels");
                        warning("at ", filename, ":", l);
                        return false;
                    }
                } else {
                    if (!ends_with(from, "\"")) {
                        warning("expected "+name+" in arcsec (\")");
                        warning("at ", filename, ":", l);
                        return false;
                    }

                    from = erase_end(from, "\"");
                    if (!from_string(from, to)) {
                        warning("could not read "+name+" in arcsec");
                        warning("at ", filename, ":", l);
                        return false;
                    }
                }

                return true;
            };

            if (r.type == "circle") {
                if (args.size() != 3) {
                    warning("ill formed 'circle' line, expecting 3 arguments, got ", args.size());
                    warning("at ", filename, ":", l);
                    continue;
                }

                double rad;
                if (!read_dim("radius", args[2], rad)) {
                    continue;
                }

                r.params = {x0, y0, rad};
            } else if (r.type == "box") {
                if (args.size() != 5) {
                    warning("ill formed 'box' line, expecting 5 arguments, got ", args.size());
                    warning("at ", filename, ":", l);
                    continue;
                }

                double width, height;
                if (!read_dim("width", args[2], width) || !read_dim("height", args[3], height)) {
                    continue;
                }

                double angle;
                if (!from_string(args[4], angle)) {
                    warning("could not read angle");
                    warning("at ", filename, ":", l);
                    continue;
                }

                r.params = {x0, y0, width, height, angle};
            } else {
                warning("support for regions of type '", r.type, "' is not implemented");
                continue;
            }

            // Read display options
            spos = line.find_first_of('#', epos);
            if (spos != line.npos) {
                std::string key = "color=";
                auto pos = line.find(key, spos+1);
                if (pos != line.npos) {
                    pos += key.size();
                    r.color = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
                }

                key = "text=";
                pos = line.find(key, spos+1);
                if (pos != line.npos) {
                    pos += key.size();
                    r.text = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
                }
            }

            regs.push_back(r);
        }
    }

    region physical_to_wcs(const astro::wcs& w, region r) {
        if (!r.physical) return r;

        double aspix = 1.0;
        astro::get_pixel_size(w, aspix);

        double ra, dec;
        astro::xy2ad(w, r.params[0]+1, r.params[1]+1, ra, dec);
        r.params[0] = ra;
        r.params[1] = dec;

        if (r.type == "circle") {
            r.params[2] *= aspix;
        } else if (r.type == "box") {
            r.params[2] *= aspix;
            r.params[3] *= aspix;
        } else
            warning("physical-to-wcs conversion for regions of type '",
                r.type, "' is not implemented");
            continue;
        }

        r.physical = false;
        return r;
    }

    region wcs_to_physical(const astro::wcs& w, region r) {
        if (r.physical) return r;

        double aspix = 1.0;
        astro::get_pixel_size(w, aspix);

        double x, y;
        astro::ad2xy(w, r.params[0], r.params[1], x, y);
        r.params[0] = x-1;
        r.params[1] = y-1;

        if (r.type == "circle") {
            r.params[2] /= aspix;
        } else if (r.type == "box") {
            r.params[2] /= aspix;
            r.params[3] /= aspix;
        } else {
            warning("wcs-to-physical conversion for regions of type '",
                r.type, "' is not implemented");
            continue;
        }

        r.physical = true;
        return r;
    }

    inline void read_regions_physical(const std::string& filename, const astro::wcs& w,
        vec<1,region>& regs) {

        read_regions(filename, regs);

        phypp_check(w.is_valid(),
            "invalid WCS, cannot convert regions to 'physical' coordinates");

        for (auto& r : regs) {
            r = wcs_to_physical(w, r);
        }
    }

    inline void read_regions_physical(const std::string& filename, vec<1,region>& regs) {
        read_regions(filename, regs);
        for (auto& r : regs) {
            phypp_check(r.physical, "expected regions in 'physical' coordinates (in ", filename, ")");
        }
    }

    inline void mask_regions(const vec<1,region>& regs, vec2b& mask) {
        phypp_check(!mask.empty(), "mask file must be initialized before calling this function");

        for (auto& r : regs) {
            phypp_check(r.physical, "regions must be in physical coordinates to create masks");

            if (r.type == "circle") {
                mask = mask || (circular_mask(mask.dims, r.params[2], r.params[1], r.params[0]) > 0.5);
            } else if (r.type == "box") {
                vec1d xo = {+0.5*r.params[2], +0.5*r.params[2], -0.5*r.params[2], -0.5*r.params[2]};
                vec1d yo = {+0.5*r.params[3], -0.5*r.params[3], -0.5*r.params[3], +0.5*r.params[3]};
                double ca = cos(r.params[4]*dpi/180.0);
                double sa = sin(r.params[4]*dpi/180.0);
                vec1d xr = xo*ca - yo*sa + r.params[0];
                vec1d yr = yo*ca + xo*sa + r.params[1];

                int_t x0 = floor(min(xr)), x1 = ceil(max(xr));
                int_t y0 = floor(min(yr)), y1 = ceil(max(yr));
                if (x1 < 0 || y1 < 0 || x0 > int_t(mask.dims[1])-1 || y0 > int_t(mask.dims[0])-1) {
                    // Not covered
                    continue;
                }

                x0 = max(x0, 0);
                y0 = max(y0, 0);
                x1 = min(x1, mask.dims[1]-1);
                y1 = min(y1, mask.dims[0]-1);

                auto hull = build_convex_hull(yr, xr);

                for (int_t iy = y0; iy <= y1; ++iy)
                for (int_t ix = x0; ix <= x1; ++ix) {
                    mask.safe(iy,ix) = mask.safe(iy,ix) || in_convex_hull(iy, ix, hull);
                }
            } else {
                warning("masks for regions of type '", r.type, "' is not implemented");
                continue;
            }
        }
    }

    inline void mask_regions(const std::string& filename, vec2b& mask) {
        vec<1,region> regs;
        read_regions_physical(filename, regs);
        mask_regions(regs, mask);
    }

    inline void mask_regions(const std::string& filename, const astro::wcs& w, vec2b& mask) {
        vec<1,region> regs;
        read_regions_physical(filename, w, regs);
        mask_regions(regs, mask);
    }
}
}

namespace impl {
namespace astro {
namespace ds9 {
    inline void write_region(std::ofstream& out, const std::string& ss, const region& r) {
        if (r.type == "circle") {
            out << r.type << "(" << r.params[0] << ", " << r.params[1] << ", " <<
                r.params[2] << ss << ")\n";
        } else if (r.type == "box") {
            out << r.type << "(" << r.params[0] << ", " << r.params[1] << ", " <<
                r.params[2] << ss << ", " <<
                r.params[3] << ss << ", " << r.params[4] << ")\n";
        } else {
            warning("writing regions of type '", r.type, "' is not implemented");
            continue;
        }
    }
}
}
}

namespace astro {
namespace ds9 {
    inline void write_regions_wcs(const std::string& filename, const vec<1,region>& regs,
        const astro::wcs& w) {

        std::ofstream regfile(filename);
        regfile << "# Region file format: DS9 version 4.1\n";
        regfile << "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal "
            "roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 "
            "source=1\n";

        bool wcs_radec = false;
        if (w.type[0] == axis_type::spatial) {
            wcs_radec = true;
        }

        std::string size_suffix;
        if (wcs_radec) {
            size_suffix = "\"";
            regfile << "fk5\n";
        } else {
            regfile << "wcs;\n";
        }

        regfile.precision(12);

        for (auto& r : regs) {
            if (r.physical) {
                region n = physical_to_wcs(w, r);
                impl::astro::ds9::write_region(regfile, size_suffix, n);
            } else {
                impl::astro::ds9::write_region(regfile, size_suffix, r);
            }
        }
    }

    inline void write_regions_physical(const std::string& filename, const vec<1,region>& regs,
        const astro::wcs& w) {

        double aspix;
        astro::get_pixel_size(w, aspix);

        std::ofstream regfile(filename);
        regfile << "# Region file format: DS9 version 4.1\n";
        regfile << "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal "
            "roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 "
            "source=1\n";
        regfile << "physical\n";

        regfile.precision(12);

        for (auto& r : regs) {
            if (r.physical) {
                impl::astro::ds9::write_region(regfile, "", r);
            } else {
                region n = wcs_to_physical(w, r);
                impl::astro::ds9::write_region(regfile, "", n);
            }
        }
    }

    inline void write_regions_physical(const std::string& filename, const vec<1,region>& regs) {
        std::ofstream regfile(filename);
        regfile << "# Region file format: DS9 version 4.1\n";
        regfile << "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal "
            "roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 "
            "source=1\n";
        regfile << "physical\n";

        regfile.precision(12);

        for (auto& r : regs) {
            phypp_check(r.physical, "no WCS provided to convert region to physical");
            impl::astro::ds9::write_region(regfile, "", r);
        }
    }
}
}
}

#endif
