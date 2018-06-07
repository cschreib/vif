#include <phypp.hpp>

bool read_ds9_region_circles(std::string file_name, vec2d& regs, bool& physical, std::string color) {
    if (!file::exists(file_name)) {
        error("could not open region file '", file_name, "'");
        return false;
    }

    std::ifstream file(file_name);

    std::string global_color = "green";
    physical = false;
    std::string line;
    uint_t l = 0;
    while (std::getline(file, line)) {
        ++l;
        if (line.empty() || trim(line).empty() || trim(line)[0] == '#') continue;

        if (start_with(line, "global")) {
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

            if (!end_with(args[2], "\"")) {
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

        append<0>(regs, vec2d{{ra, dec, rad}});
    }

    return true;
}

bool read_ds9_region_circles_physical(std::string file_name,
    std::string img_name, vec2d& regs, std::string color = "") {

    bool physical = false;
    if (!read_ds9_region_circles(file_name, regs, physical, color)) return false;

    if (physical) return true;

    astro::wcs w(fits::read_header(img_name));
    if (!w.is_valid()) {
        error("could not read astrometry from '", img_name, "'");
        return false;
    }

    for (uint_t i : range(regs.dims[0])) {
        double x, y;
        astro::ad2xy(w, regs(i,0), regs(i,1), x, y);
        regs(i,0) = x-1;
        regs(i,1) = y-1;
    }

    double aspix = 1.0;
    if (!astro::get_pixel_size(w, aspix)) {
        return false;
    }

    regs(_,2) /= aspix;

    return true;
}

int phypp_main(int argc, char* argv[]) {
    std::string aper_reg = argv[1];
    std::string bg_reg   = argv[2];
    std::string out_file = argv[3];

    if (argc <= 4) return 0;

    vec1s imgfile;
    vec1d apcor;
    vec1d flux;
    vec1d flux_err;
    vec1s band;
    vec1u eazy_band;
    vec1f lambda;

    // Read filter data base
    auto fdb = read_filter_db(data_dir+"fits/filter-db/db.dat");
    auto fmap = read_filter_map(data_dir+"fits/filter-db/fast.map");

    std::vector<std::pair<std::string,std::string>> bandmap = {
        {"kpno-u", "mosaic_u"},
        {"vimos-u", "vimos_u"},
        {"vimos-r", "vimos_r"},
        {"cfht-u", "cfht_u"},
        {"cfht-g", "cfht_g"},
        {"cfht-r", "cfht_r"},
        {"cfht-i", "cfht_i"},
        {"cfht-z", "cfht_z"},
        {"mosaic-z", "cfht_z"},
        {"subaru-B", "sub_B"},
        {"subaru-V", "sub_V"},
        {"subaru-g", "sub_g"},
        {"subaru-r", "sub_r"},
        {"subaru-R", "sub_r"},
        {"subaru-i", "sub_i"},
        {"subaru-z", "sub_z"},
        {"subaru-zpp", "sub_z"},
        {"subaru-ib427", "IB427"},
        {"subaru-ib445", "IB445"},
        {"subaru-ib464", "IB464"},
        {"subaru-ib484", "IB484"},
        {"subaru-ib505", "IB505"},
        {"subaru-ib527", "IB527"},
        {"subaru-ib550", "IB550"},
        {"subaru-ib574", "IB574"},
        {"subaru-ib598", "IB598"},
        {"subaru-ib624", "IB624"},
        {"subaru-ib651", "IB651"},
        {"subaru-ib679", "IB679"},
        {"subaru-ib709", "IB709"},
        {"subaru-ib738", "IB738"},
        {"subaru-ib767", "IB767"},
        {"subaru-ib827", "IB827"},
        {"subaru-ib856", "IB856"},
        {"wfi-U", "wfi_U"},
        {"wfi-U32", "wfi_U38"},
        {"wfi-B", "wfi_B"},
        {"wfi-V", "wfi_V"},
        {"wfi-R", "wfi_R"},
        {"wfi-I", "wfi_I"},
        {"acs-f435w", "f435w"},
        {"acs-f606w", "f606w"},
        {"acs-f775w", "f775w"},
        {"acs-f814w", "f814w"},
        {"acs-f850lp", "f850lp"},
        {"wfc3-f105w", "f105w"},
        {"wfc3-f125w", "f125w"},
        {"wfc3-f140w", "f140w"},
        {"wfc3-f160w", "f160w"},
        {"vista-Y", "vista_Y"},
        {"vista-J", "vista_J"},
        {"vista-nb118", "vista_NB118"},
        {"vista-H", "vista_H"},
        {"vista-Ks", "vista_Ks"},
        {"hawki-Y", "hawki_Y"},
        {"hawki-k", "hawki_Ks"},
        {"hawki-Ks", "hawki_Ks"},
        {"wircam-Y", "wircam_Y"},
        {"wircam-J", "wircam_J"},
        {"wircam-H", "wircam_H"},
        {"wircam-Ks", "wircam_Ks"},
        {"isaac-J", "isaac_J"},
        {"isaac-H", "isaac_H"},
        {"isaac-Ks", "isaac_Ks"},
        {"cfht-Ks", "wircam_Ks"},
        {"moircs-Y", "moircs_Y"},
        {"moircs-J", "moircs_J"},
        {"moircs-H", "moircs_H"},
        {"moircs-Ks", "moircs_Ks"},
        {"ukidss-J", "ukidss_J"},
        {"ukidss-H", "ukidss_H"},
        {"ukidss-K", "ukidss_K"},
        {"zfourge-J", "fourstar_J"},
        {"zfourge-J1", "fourstar_J1"},
        {"zfourge-J2", "fourstar_J2"},
        {"zfourge-J3", "fourstar_J3"},
        {"zfourge-H", "fourstar_H"},
        {"zfourge-H1", "fourstar_H1"},
        {"zfourge-H2", "fourstar_H2"},
        {"zfourge-Ks", "fourstar_Ks"},
        {"spitzer-irac-3p6", "i1"},
        {"spitzer-irac-4p5", "i2"},
        {"spitzer-irac-5p8", "i3"},
        {"spitzer-irac-8p0", "i4"}
    };

    // Sort bandmap by decreasing match pattern size, to avoid too early matches
    // (e.g.: "subaru-ib427" matching "subaru-i")
    std::sort(bandmap.begin(), bandmap.end(),
        [](const std::pair<std::string,std::string>& p1, const std::pair<std::string,std::string>& p2) {
        return p1.first.size() > p2.first.size();
    });

    // Extract fluxes
    for (int ii : range(argc-4)) {
        // Filter out special images
        std::string imgf = argv[ii+4];
        if (end_with(imgf, "-psf.fits")) continue;

        // Read image
        fits::input_image fimg(imgf);
        if (fimg.axis_count() != 2) continue;

        imgfile.push_back(imgf);
        print(imgf);

        fits::header hdr = fits::read_header(imgf);
        std::string base = file::remove_extension(file::get_basename(imgf));
        std::string tmp_dir = file::get_directory(imgf)+"apers/";
        file::mkdir(tmp_dir);

        vec2d img;
        fimg.read(img);
        vec2d xx = generate_img(img.dims, [](double, double tx) { return tx; });
        vec2d yy = generate_img(img.dims, [](double ty, double) { return ty; });

        // Read apertures
        vec2d aper, bg, bad;

        // Use custom aperture and background if they exist
        std::string aper_file = aper_reg;
        std::string bg_file = bg_reg;
        if (file::exists(file::remove_extension(imgf)+"-aper.reg")) {
            aper_file = file::remove_extension(imgf)+"-aper.reg";
        }
        if (file::exists(file::remove_extension(imgf)+"-bg.reg")) {
            bg_file = file::remove_extension(imgf)+"-bg.reg";
        }

        if (!read_ds9_region_circles_physical(aper_file, imgfile.back(), aper)) {
            return 1;
        }
        if (!read_ds9_region_circles_physical(bg_file, imgfile.back(), bg)) {
            return 1;
        }

        // Read badpixel masks if any
        std::string bad_file = file::remove_extension(imgf)+"-bad.reg";
        if (file::exists(bad_file) && !read_ds9_region_circles_physical(bad_file, imgfile.back(), bad)) {
            bad.clear();
        }

        // Read fit regions if any
        std::string fit_region_file = file::remove_extension(imgf)+"-fit.reg";

        vec2d fit_center, fit_self, fit_region;
        if (file::exists(fit_region_file)) {
            if (!read_ds9_region_circles_physical(fit_region_file, imgfile.back(), fit_center, "green")) {
                fit_center.clear();
            }
            if (!read_ds9_region_circles_physical(fit_region_file, imgfile.back(), fit_self, "yellow")) {
                fit_center.clear();
                fit_self.clear();
            }
            if (!read_ds9_region_circles_physical(fit_region_file, imgfile.back(), fit_region, "red")) {
                fit_center.clear();
                fit_self.clear();
                fit_region.clear();
            }
        }

        // Build PSF
        auto make_psf = [&](double y0, double x0) {
            vec2d psf;
            std::string psf_file = file::remove_extension(imgf)+"-psf.fits";
            if (!file::exists(psf_file)) {
                if (imgf.find("spitzer-irac") != imgf.npos) {
                    if (imgf.find("3p6") != imgf.npos) {
                        psf_file = "/data2/share/egg/psfs/spitzer-irac1.fits";
                    } else if (imgf.find("4p5") != imgf.npos) {
                        psf_file = "/data2/share/egg/psfs/spitzer-irac2.fits";
                    } else if (imgf.find("5p8") != imgf.npos) {
                        psf_file = "/data2/share/egg/psfs/spitzer-irac3.fits";
                    } else /* if (imgf.find("8p0") != imgf.npos) */ {
                        psf_file = "/data2/share/egg/psfs/spitzer-irac4.fits";
                    }
                }
            }

            if (file::exists(psf_file)) {
                int_t ix0 = round(x0), iy0 = round(y0);
                double dx = x0 - ix0, dy = y0 - iy0;

                vec2d tpsf = translate(fits::read(psf_file), dy, dx);
                int_t hy = tpsf.dims[0]/2, hx = tpsf.dims[1]/2;
                psf.resize(img.dims);

                vec1u ip, it;
                subregion(psf, {iy0-hy,ix0-hx,iy0+hy,ix0+hx}, ip, it);
                psf[ip] = tpsf[it];
            } else if (imgf.find("acs-f") != imgf.npos || imgf.find("wfc3-f") != imgf.npos) {
                // HST: assume ultra sharp PSF, hence no aperture correction
                psf.resize(img.dims);
                psf(int_t(round(y0)),int_t(round(x0))) = 1.0;
            } else {
                // Assume ground-based with seeing of ~1" FWHM
                double aspix = 1.0;
                astro::wcs w(fits::read_header(imgf));
                astro::get_pixel_size(w, aspix);

                double fwhm = 1.0;
                if (imgf.find("wfi-") != imgf.npos) {
                    fwhm = 1.2;
                } else if (imgf.find("hawki-") != imgf.npos) {
                    fwhm = 0.5;
                } else if (imgf.find("wircam-") != imgf.npos) {
                    fwhm = 0.8;
                }

                double width = fwhm/2.335/aspix;
                psf = exp(-(sqr(xx - x0) + sqr(yy - y0))/(2.0*sqr(width)))/(2.0*dpi*sqr(width));
            }

            return psf;
        };

        // Build background mask
        vec2d bg_mask(img.dims);
        for (uint_t i : range(bg.dims[0])) {
            bg_mask += circular_mask(img.dims, bg(i,2), bg(i,1), bg(i,0));
        }

        bg_mask = clamp(bg_mask, 0, 1);
        fits::write(tmp_dir+base+"-bg.fits", bg_mask, hdr);

        // Build badpixel mask
        vec2d bad_mask(img.dims);
        for (uint_t i : range(bad.dims[0])) {
            bad_mask += circular_mask(img.dims, bad(i,2), bad(i,1), bad(i,0));
        }

        bad_mask = clamp(bad_mask + vec2d{!is_finite(img)}, 0, 1);
        img[where(bad_mask > 0.1)] = 0;
        fits::write(tmp_dir+base+"-badpix.fits", bad_mask, hdr);

        // Make flux mask
        vec2d aper_mask(img.dims);
        for (uint_t i : range(aper.dims[0])) {
            aper_mask += circular_mask(img.dims, aper(i,2), aper(i,1), aper(i,0));
        }

        aper_mask = clamp(aper_mask, 0, 1);
        aper_mask *= (1.0 - bad_mask);
        fits::write(tmp_dir+base+"-aper.fits", aper_mask, hdr);

        double apertot = total(aper_mask);

        // Compute aperture correction
        // First find the centroid
        double y0 = total(yy*aper_mask)/apertot;
        double x0 = total(xx*aper_mask)/apertot;
        if (!is_finite(y0)) y0 = img.dims[0]/2;
        if (!is_finite(x0)) x0 = img.dims[1]/2;

        // Then build the sources PSF at the right position and compute correction
        vec2d psf = make_psf(y0, x0);
        apcor.push_back(1/total(psf*aper_mask));
        fits::write(tmp_dir+base+"-psf.fits", psf, hdr);

        // Fit and subtract neighbors, if any
        if (fit_center.dims[0] != 0) {
            // Build fit mask
            vec2d fit_mask(img.dims);
            for (uint_t i : range(fit_region.dims[0])) {
                fit_mask += circular_mask(img.dims, fit_region(i,2), fit_region(i,1), fit_region(i,0));
            }

            fit_mask = clamp(fit_mask, 0, 1);
            vec1u idf = where(fit_mask > 0.5 && bad_mask < 0.1);

            if (!idf.empty()) {
                vec2d psfs(2+fit_center.dims[0], psf.size());
                psfs(0,_) = 1;

                if (fit_self.empty()) {
                    psfs(1,_) = flatten(psf);
                } else {
                    for (uint_t i : range(fit_self.dims[0])) {
                        psfs(1,_) += flatten(make_psf(fit_self(i,1), fit_self(i,0)));
                    }

                    psfs(1,_) /= fit_self.dims[0];
                }

                for (uint_t i : range(fit_center.dims[0])) {
                    psfs(2+i,_) = flatten(make_psf(fit_center(i,1), fit_center(i,0)));
                }

                // Fit and discard sources until all fluxes are positive
                linfit_result res;
                do {
                    if (!res.params.empty()) {
                        vec1u idp = where(res.params > 0 || uindgen(psfs.dims[0]) <= 1);
                        psfs = psfs(idp,_);
                    }

                    res = linfit_pack(img[idf], replicate(1.0, idf.size()), psfs(_,idf));
                } while (count(res.params < 0 && uindgen(psfs.dims[0]) > 1) > 0);

                // Subtract
                for (uint_t i : range(psfs.dims[0]-2)) {
                    img -= reform(res.params[2+i]*psfs(2+i,_), img.dims);
                }

                fits::write_table(tmp_dir+base+"-fit.fits", "psfs", psfs, "flux", res.params);

                // Subtract source and make residual
                vec2f timg = img;
                timg -= reform(res.params[1]*psfs(1,_), img.dims);
                fits::write(tmp_dir+base+"-res.fits", timg, hdr);
            } else {
                fits::write(tmp_dir+base+"-res.fits", img, hdr);
            }
        }

        // Subtract background from image
        double bgtot = total(bg_mask);
        if (bgtot != 0.0) {
            vec1u idb = where(bg_mask > 0 && bad_mask < 0.1);
            if (!idb.empty()) {
                img -= median(img[idb]);
            }
        } else {
            error("no pixel belongs to the background");
            note("processing ", imgf);
            return 1;
        }

        fits::write(tmp_dir+base+"-sub.fits", img, hdr);

        // Compute RMS
        vec1u idbg = where(bg_mask > 0.5 && bad_mask < 0.1);
        double rms = stddev(img[idbg]);

        // Take into account the uncertainty on the background estimation
        vec1d bgs(bg.dims[0]);
        for (uint_t i : range(bg.dims[0])) {
            vec2d tbg_mask = circular_mask(img.dims, bg(i,2), bg(i,1), bg(i,0));
            vec1u ids = where(tbg_mask > 0 && bad_mask < 0.1);
            if (!ids.empty()) {
                bgs[i] = median(img[ids]);
            } else {
                bgs[i] = dnan;
            }
        }

        double bg_err = stddev(bgs[where(is_finite(bgs))]);

        // Compute flux
        flux.push_back(total(img*aper_mask));
        flux_err.push_back(sqrt(apertot*sqr(rms) + sqr(apertot*bg_err)));

        // Apply aperture correction
        flux.back() *= apcor.back();
        flux_err.back() *= apcor.back();

        // Try to find a band ID
        band.push_back("");
        eazy_band.push_back(-1);
        lambda.push_back(fnan);

        for (auto p : bandmap) {
            if (imgf.find(p.first) != imgf.npos) {
                band.back() = p.second;
                break;
            }
        }

        if (!band.back().empty()) {
            filter_t fil;
            if (!get_filter(fdb, band.back(), fil)) {
                warning("unknown filter '", band.back(), "'");
            } else {
                lambda.back() = fil.rlam;
                get_filter_id(fmap, band.back(), eazy_band.back());
            }
        }
    }

    fits::write_table(out_file, ftable(imgfile, apcor, flux, flux_err, band, lambda, eazy_band));

    return 0;
}
