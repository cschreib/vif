#include <phypp.hpp>
#include <phypp/astro/qxmatch.hpp>

void print_help();

int main(int argc, char* argv[]) {
    std::string cat_file; // Name of the catalog from which to take the sources
    std::string img_file; // Name of the image from which to measure upper limits
    std::string err_file; // Name of the associated error map
    std::string psf_file; // Name of the associated PSF profile
    std::string out_file; // Name of the output catalog
    std::string res_file; // Name of the output residual map (optional)
    std::string grp_file; // Name of the output group map (optional)
    std::string mod_file; // Name of the output model map (optional)

    // Conversion factor from map unit to uJy x PSF unit
    float fconv = 1.0;
    // Print some info about the process
    bool verbose = false;
    // Freeze the background to a given value
    float fixed_bg = fnan;
    // Convolve the PSF
    bool beam_smeared = false;
    // If 'beam_smeared' is 'true', then either use the PSF itself as the beam
    // (default, beam_size=dnan), or use a custom Gaussian beam of sigma 'beam_size'
    double beam_size = dnan;
    // Save the covariance matrix
    bool save_covariance = false;
    // Do not compute the whole covariance matrix, just invert some sub-matrices
    bool cell_approx = false;
    // Use a flux prior
    bool flux_prior = false;
    // Reduce image fit area to regions covered by the input prior list
    bool trim_image = false;
    // Group uncertain sources into single measurements
    bool make_groups = false;
    double group_fit_threshold = 0.97;
    double group_aper_threshold = 0.8;
    double group_aper_size = dnan;
    double group_cov_threshold = 0.4;
    bool group_post_process = false;
    // Number of threads the program can use
    uint_t nthread = 1;
    // Display help
    bool help = false;

    read_args(argc, argv, arg_list(
        name(cat_file, "cat"), name(img_file, "img"), name(psf_file, "psf"),
        name(err_file, "err"), name(out_file, "out"), name(res_file, "out_res"),
        name(mod_file, "out_model"), name(grp_file, "out_gmap"), help,
        fconv, verbose, fixed_bg, save_covariance, cell_approx, flux_prior,
        make_groups, group_fit_threshold, group_aper_threshold, group_aper_size,
        group_cov_threshold, group_post_process, name(nthread, "threads"),
        beam_smeared, beam_size, trim_image
    ));

    if (argc == 1 || help) {
        print_help();
        return 1;
    }

    if (cell_approx && save_covariance) {
        cell_approx = false;
    }

    if (out_file.empty()) {
        error("missing output file name (out=...)");
        return 1;
    }

    if (psf_file.empty()) {
        if (end_with(img_file, "_sci.fits")) {
            psf_file = erase_end(img_file, "_sci.fits")+"_psf.fits";

            if (!file::exists(psf_file)) {
                error("tried guessing the PSF file to '"+psf_file+"', but this file "
                    "does not exist.");
                note("please provide the name of the PSF file manually (psf=...)");
            }
        } else {
            error("please provide the name of the PSF file to use (psf=...)");
            return 1;
        }
    }

    if (err_file.empty()) {
        if (end_with(img_file, "_sci.fits")) {
            err_file = erase_end(img_file, "_sci.fits")+"_err.fits";

            if (!file::exists(err_file)) {
                error("tried guessing the error map file to '"+err_file+"', but this file "
                    "does not exist.");
                note("please provide the name of the error map file manually (err=...)");
            }
        } else {
            error("please provide the name of the error map file to use (err=...)");
            return 1;
        }
    }

    vec1d ra, dec;
    fits::read_table(cat_file, ftable(ra, dec));

    vec1f fprior, fprior_err;
    if (flux_prior) {
        fits::read_table(cat_file, ftable(fprior, fprior_err));
    }

    vec2d psf;
    fits::read(psf_file, psf);

    // Trim the PSF, make sure that the peak is at the center, and that the dimensions
    // are odds
    int_t hsize = 0; {
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

        if (beam_smeared) {
             vec2d kernel;
            if (is_finite(beam_size)) {
                kernel = gaussian_profile(psf.dims, beam_size);
            } else {
                kernel = psf;
            }
            vec1u idz = where(kernel < 1e-5*max(abs(kernel)));
            kernel[idz] = 0.0;
            psf = convolve2d(psf, kernel);
            psf /= psf(hsize,hsize);
        }
    }


    if (verbose) {
        print("will extract ", ra.size(), " sources from this map");
    }

    // Read image
    if (verbose) {
        print("reading map in memory...");
    }

    vec2d img, err;
    fits::header hdr;
    fits::read(img_file, img, hdr);
    fits::read(err_file, err);

    // If asked, identify grouped sources
    struct {
        vec1d ra, dec;
        vec1f fprior, fprior_err;
        vec1u group_fit_id;
        vec1u group_aper_id;
    } old_cat;

    struct {
        vec1u id;
        vec1d ra, dec;
        vec1f fprior, fprior_err;
        vec1f flux, flux_err;
        vec1b fit;
        vec1u nsrc;
    } group_cat;

    vec1u id_new;
    vec1u id_old;
    vec1u id_grfit;
    vec1u group_fit_id;
    vec1b is_grouped;
    bool has_groups = false;

    if (make_groups) {
        // Save original catalog
        old_cat.ra = ra;
        old_cat.dec = dec;
        old_cat.fprior = fprior;
        old_cat.fprior_err = fprior_err;

        // Compute physical size of the requested aperture
        double aspix;
        if (!fits::get_pixel_size(img_file, aspix)) {
            error("could not read pixel size from '", img_file, "', missing WCS?");
            return 1;
        }

        if (!is_finite(group_aper_size)) {
            // No aperture specified, chose the aperture as the radius outside of which
            // the amplitude of the PSF is less than 50%
            double ap_th = 0.5*max(psf);
            for (uint_t i : range(hsize)) {
                if (psf(hsize,i) > ap_th || psf(hsize,2*hsize-i) > ap_th ||
                    psf(i,hsize) > ap_th || psf(2*hsize-i,hsize) > ap_th) {
                    group_aper_size = hsize-i;
                    break;
                }
            }
        } else {
            group_aper_size /= aspix;
        }

        if (verbose) {
            print("group aperture size: ", group_aper_size, " pixels");
        }

        // Calibrate the relation between distance and PSF overlap
        vec1d d = uindgen(2*hsize+2);
        vec1d c(d.size());
        for (uint_t i : range(d)) {
            c[i] = total(psf*translate(psf, d[i], 0.0));
        }

        c /= c[0];
        d *= aspix;

        // Now build the groups
        old_cat.group_fit_id.resize(ra.size());
        old_cat.group_aper_id.resize(ra.size());
        uint_t last_group_fit_id = 0;
        uint_t last_group_aper_id = 0;

        auto group_pair = [](vec1u& grp, uint_t& last_id, uint_t i, uint_t j) {
            if (grp[i] != 0) {
                // Source 'i' has a group
                if (grp[j] != 0) {
                    // Source 'j' too
                    // Merge the two groups and keep only the one with the smallest ID
                    if (grp[i] < grp[j]) {
                        vec1u id = where(grp == grp[j]);
                        grp[id] = grp[i];
                    } else {
                        vec1u id = where(grp == grp[i]);
                        grp[id] = grp[j];
                    }
                } else {
                    // Source 'j' has no group
                    // Add 'j' to the group of 'i'
                    grp[j] = grp[i];
                }
            } else if (grp[j] != 0) {
                // Source 'i' has no group, but 'j' has one
                // Add 'i' to the group of 'j'
                grp[i] = grp[j];
            } else {
                // Both sources have no group
                // Create a new group for 'i' and 'j'
                ++last_id;
                grp[i] = grp[j] = last_id;
            }
        };

        for (uint_t i : range(ra))
        for (uint_t j : range(i+1, ra.size())) {
            double td = angdist(ra[i], dec[i], ra[j], dec[j]);
            double tc = interpolate(c, d, td);
            if (tc >= group_aper_threshold) {
                // These two sources are relatively close and the deblending is uncertain.
                // We therefore group them into a single object. In the following, we still
                // fit these objects separately, but we do not keep the measured fluxes.
                // Instead we point all the sources that make a single object into a
                // fictional source, added to the catalog, and whose flux is measured
                // inside an irregular aperture.
                group_pair(old_cat.group_aper_id, last_group_aper_id, i, j);

                if (tc >= group_fit_threshold) {
                    // These two sources are way too close and have no hope of being
                    // deblended, they will probably crash the fitting procedure.
                    // We therefore group them into a single object, and we assume in the
                    // following that it can be considered a single PSF.
                    group_pair(old_cat.group_fit_id, last_group_fit_id, i, j);
                }
            }
        }

        has_groups = last_group_aper_id != 0;

        if (last_group_fit_id != 0) {
            // Physically group 'group_fit' sources into a single PSF
            ra.clear();
            dec.clear();
            fprior.clear();
            fprior_err.clear();

            id_new.resize(old_cat.ra.size());

            vec1u gfid = old_cat.group_fit_id;
            for (uint_t i : range(old_cat.ra)) {
                if (gfid[i] == 0) {
                    id_new[i] = ra.size();
                    id_old.push_back(i);

                    ra.push_back(old_cat.ra[i]);
                    dec.push_back(old_cat.dec[i]);
                    is_grouped.push_back(old_cat.group_aper_id[i] != 0);
                    group_fit_id.push_back(npos);
                    if (flux_prior) {
                        fprior.push_back(old_cat.fprior[i]);
                        fprior_err.push_back(old_cat.fprior_err[i]);
                    }
                } else if (gfid[i] != npos) {
                    vec1u id = where(gfid == old_cat.group_fit_id[i]);
                    gfid[id] = npos;

                    id_new[id] = ra.size();
                    id_old.push_back(npos);
                    id_grfit.push_back(ra.size());

                    uint_t gid = group_cat.ra.size()+1;
                    double mra = mean(old_cat.ra[id]);
                    double mdec = mean(old_cat.dec[id]);

                    group_cat.id.push_back(gid);
                    old_cat.group_fit_id[id] = gid;

                    group_cat.ra.push_back(mra);
                    group_cat.dec.push_back(mdec);
                    group_cat.nsrc.push_back(id.size());
                    group_cat.fit.push_back(true);
                    ra.push_back(mra);
                    dec.push_back(mdec);
                    is_grouped.push_back(true);
                    group_fit_id.push_back(gid);
                    if (flux_prior) {
                        double mfp = total(old_cat.fprior[id]);
                        double mfpe = sqrt(total(sqr(old_cat.fprior_err[id])));
                        group_cat.fprior.push_back(mfp);
                        group_cat.fprior_err.push_back(mfpe);
                        fprior.push_back(mfp);
                        fprior_err.push_back(mfpe);
                    }
                }
            }

            if (verbose) {
                print("some sources are too close to one another");
                print("only ", count(group_fit_id == npos), "/", old_cat.ra.size(), " sources will "
                    "actually be fitted, and ", count(group_fit_id != npos), " fit groups were "
                    "created");
            }
        } else {
            is_grouped = old_cat.group_aper_id != 0;
            group_fit_id = replicate(npos, ra.size());
            id_old = uindgen(ra.size());
            id_new = uindgen(ra.size());
        }

        if (has_groups) {
            // Save the 'group_aper' groups
            uint_t ngrp = 0;
            vec1u gaid = old_cat.group_aper_id;
            for (uint_t i : range(old_cat.ra)) {
                if (gaid[i] == 0 || gaid[i] == npos) continue;

                vec1u id = where(gaid == old_cat.group_aper_id[i]);
                gaid[id] = npos;

                uint_t gid = group_cat.ra.size()+1;
                double mra = mean(old_cat.ra[id]);
                double mdec = mean(old_cat.dec[id]);

                group_cat.id.push_back(gid);
                old_cat.group_aper_id[id] = gid;

                group_cat.ra.push_back(mra);
                group_cat.dec.push_back(mdec);
                group_cat.nsrc.push_back(id.size());
                group_cat.fit.push_back(false);
                if (flux_prior) {
                    group_cat.fprior.push_back(0.0);
                    group_cat.fprior_err.push_back(0.0);
                }

                ++ngrp;
            }

            if (verbose) {
                print("created ", ngrp, " source groups");
            }
        }
    }

    // Subtract background if requested
    bool free_bg = !is_finite(fixed_bg);
    if (!free_bg) {
        img -= fixed_bg;
    }

    // Convert to flux unit
    img *= fconv;
    err *= fconv;

    double err_nocov; {
        vec1u idg = where(is_finite(err));
        err_nocov = 1e5*max(err[idg]);
    }

    vec2d snr = img/err; {
        vec1u idb = where(!is_finite(snr) || err == 0.0);
        snr[idb] = 0.0;
        err[idb] = err_nocov;
    }

    // Compute X and Y position of all the sources from embedded astrometry
    if (verbose) {
        print("compute astrometry...");
    }

    uint_t nsrc = ra.size();
    vec1d x, y;
    fits::ad2xy(fits::extast(hdr), ra, dec, x, y);
    x -= 1; y -= 1;

    // Only keep the sources which fall on the map
    vec1u idin = where(x >= -hsize && x < img.dims[1]+hsize &&
                       y >= -hsize && y < img.dims[0]+hsize);

    uint_t nobs = idin.size();

    if (verbose) {
        uint_t nout = nsrc - nobs;
        if (nout > 0) {
            note(nout, " sources (", 100.0*round(nout/nsrc), "%) are not observable on "
                "this map and will therefore have no flux measurement");
        }
    }

    x = x[idin]; y = y[idin];

    if (flux_prior) {
        fprior = fprior[idin];
        fprior_err = fprior_err[idin];
    }

    vec1i ix = round(x), iy = round(y);
    vec1d dx = x - ix, dy = y - iy;

    // Adjust image to flag out the areas non covered by the input prior catalog
    // Use this option if your catalog only covers a small part of the image, else
    // it is more stable to use the the full image and a complete prior list
    if (trim_image) {
        vec2b incov;
        if (x.size() > 2) {
            // Pick the pixels of the map that are within the coverage of our prior list
            auto hull = build_convex_hull(x, y);
            vec2d tix = replicate(dindgen(snr.dims[1]), snr.dims[0]);
            vec2d tiy = transpose(replicate(dindgen(snr.dims[0]), snr.dims[1]));
            incov = in_convex_hull(tix, tiy, hull);
        }

        if (count(incov) < 5) {
            // Use all pixels
            incov[_] = true;
        }

        vec1u idb = where(!incov);
        snr[idb] = 0;
        err[idb] = err_nocov;

        if (verbose) {
            print(count(incov), " pixels used in the image (",
                round(100*fraction_of(incov)), "% of the available area)");
        }
    }

    // Compute the matrix
    if (verbose) {
        print("compute matrix...");
    }

    uint_t nelem = nobs + (free_bg ? 1 : 0);

    vec2d alpha(nelem, nelem);
    vec1d beta(nelem);

    // Alpha is symmetric, so only compute one side
    // This matrix measures the overlap between the different fit components
    // Beta measures the product of each component with the actual data

    // Source terms
    vec1f local_error(nsrc);

    auto pg = progress_start(nobs*(nobs-1)/2);
    for (uint_t i : range(nobs)) {
        // TODO: for groups, build a combined PSF instead of just using a PSF at the center

        // Get the weighted PSF of source 'i'
        vec2d tpsf2 = translate(psf, dy[i], dx[i]);
        vec1u idi, idp;
        subregion(snr, {iy[i]-hsize, ix[i]-hsize, iy[i]+hsize, ix[i]+hsize}, idi, idp);

        vec1u idpc = complement(tpsf2, idp);
        vec1f terr = err[idi];
        tpsf2[idp] /= terr;
        tpsf2[idpc] = 0;
        vec1d tpsf = tpsf2[idp];

        // Compute the local RMS of the map
        local_error[idin[i]] = 1.0/sqrt(total(sqr(tpsf)));

        // Beta term: beta[i] = x[i]*image/err^2
        beta[i] = total(snr[idi]*tpsf);

        // Alpha terms
        // The source with itself: alpha(i,i) = (x[i]/err)^2
        alpha(i,i) = total(sqr(tpsf));

        if (flux_prior) {
            // If requested, add a prior on the flux of each source
            beta[i] += fprior[i]/sqr(fprior_err[i]);
            alpha(i,i) += 1.0/sqr(fprior_err[i]);
        }

        if (free_bg) {
            // Source x Background: alpha(i,bg) = x[i]/err^2
            alpha(i,nobs) = alpha(nobs,i) = total(tpsf/terr);
        }

        // Source x Source: alpha(j,i) = x[i]*x[j]/err^2
        for (uint_t j : range(i+1, nobs)) {
            int_t idx = ix[i]-ix[j], idy = iy[i]-iy[j];

            if (abs(idx) <= 2*hsize && abs(idy) <= 2*hsize) {
                vec1u pidi, pidj;
                subregion(psf, {idy, idx, idy+2*hsize, idx+2*hsize}, pidj, pidi);
                if (!pidi.empty()) {
                    // Get the weighted PSF of source 'j'
                    vec2d tpsf2_j; {
                        tpsf2_j = translate(psf, dy[j], dx[j]);
                        vec1u idi_j, idp_j;
                        subregion(snr, {iy[j]-hsize, ix[j]-hsize, iy[j]+hsize, ix[j]+hsize},
                            idi_j, idp_j);

                        vec1u idpc_j = complement(tpsf2_j, idp_j);
                        tpsf2_j[idp_j] /= err[idi_j];
                        tpsf2_j[idpc_j] = 0;
                    }

                    alpha(j,i) = alpha(i,j) = total(tpsf2_j[pidj]*tpsf2[pidi]);
                }
            }

            if (verbose) progress(pg, 1123);
        }
    }

    // Pure Background terms
    if (free_bg) {
        // Beta term: beta[bg] = image/err^2
        beta[nobs] = total(snr/err);

        // Background x Background: alpha(bg,bg) = 1/err^2
        alpha(nobs, nobs) = total(1.0/sqr(err));
    }

    // Solve the system
    vec1d best_fit, best_fit_err;
    vec1f flux(nsrc);
    vec1f flux_err(nsrc);
    double background, background_err;

    auto save_fit_basics = [&]() {
        if (make_groups) {
            if (!id_new.empty()) {
                // Save grouped fluxes
                group_cat.flux = flux[id_grfit];
                group_cat.flux_err = flux_err[id_grfit];

                // Ungroup the grouped sources
                ra = old_cat.ra;
                dec = old_cat.dec;
                flux = flux[id_new];
                flux_err = flux_err[id_new];
                local_error = local_error[id_new];

                // Give them no flux, to prevent confusion
                vec1u idg = where(old_cat.group_fit_id != 0);
                flux[idg] = 0.0;
                flux_err[idg] = 0.0;
            }
        }

        // Write the result to a file
        fits::write_table(out_file, ftable(
            ra, dec, background, background_err, flux, flux_err, local_error
        ));
    };

    vec2d covar;
    if (has_groups) {
        covar = alpha;
        matrix::symmetrize(covar);
    }

    if (cell_approx) {
        if (verbose) {
            print("compute approximated covariance errors...");
        }

        // Estimate the error for each source
        vec1u idn; idn.reserve(nobs);
        best_fit_err.resize(nelem);
        for (uint i : range(nobs)) {
            idn.clear();

            uint_t ii = npos;
            for (uint_t j : range(nobs)) {
                if (alpha(i,j) > 1e-3*sqrt(alpha(i,i)*alpha(j,j))) {
                    if (j == i) ii = idn.size();
                    idn.push_back(j);
                }
            }

            uint_t neib = idn.size();
            uint_t ntelem = neib + (free_bg ? 1 : 0);
            vec2d talpha(ntelem, ntelem);
            for (uint_t ti : range(neib)) {
                if (free_bg) {
                    talpha(ti,neib) = talpha(neib,ti) = alpha(nobs,idn[ti]);
                }

                for (uint_t tj : range(ti, neib)) {
                    talpha(tj,ti) = alpha(idn[tj],idn[ti]);
                    if (tj != ti) {
                        talpha(ti,tj) = talpha(tj,ti);
                    }
                }
            }

            if (free_bg) {
                talpha(neib,neib) += alpha(nobs,nobs);
            }

            if (!matrix::inplace_invert_symmetric(talpha)) {
                error("could not invert sub-matrix for source ", idin[i]);
                note("there are probably some prior source positions which are too "
                    "close and cannot be deblended");
                return 1;
            } else {
                double terr = sqrt(talpha(ii,ii));
                best_fit_err[i] = terr;
                flux_err[idin[i]] = terr;
            }
        }

        // Now solve the system to get the fluxes
        if (verbose) {
            print("solve system...");
        }

        if (!matrix::inplace_solve_symmetric(alpha, beta)) {
            error("could not solve the linear problem");
            note("there are probably some prior source positions which are too close and "
                "cannot be deblended");
            return 1;
        }

        // The best fit values are now stored in beta
        best_fit = beta;

        // Extract the background value if needed
        if (free_bg) {
            background = best_fit[nobs]/fconv;
        } else {
            background = fixed_bg;
        }

        background_err = 0.0;

        // Extract the fluxes
        flux[idin] = best_fit[_-(nobs-1)];

        // Save to disk
        save_fit_basics();
    } else {
        if (verbose) {
            print("invert matrix...");
        }

        if (!matrix::inplace_invert_symmetric(alpha)) {
            error("could not invert covariance matrix, it is singular");
            note("there are probably some prior source positions which are too close and "
                "cannot be deblended");
            return 1;
        }

        // Compute covariance matrix to get the errors
        matrix::symmetrize(alpha);

        // Multiply the inverted alpha with beta to get the best fit values
        best_fit = matrix::product(alpha, beta);
        best_fit_err = sqrt(matrix::diagonal(alpha));

        // Extract the background value if needed
        if (free_bg) {
            background = best_fit[nobs]/fconv;
            background_err = best_fit_err[nobs]/fconv;
        } else {
            background = fixed_bg;
            background_err = 0.0;
        }

        // Extract the fluxes and associated errors
        flux[idin] = best_fit[_-(nobs-1)];
        flux_err[idin] = best_fit_err[_-(nobs-1)];

        // Save to disk
        save_fit_basics();

        if (save_covariance) {
            vec2d covariance(nsrc, nsrc);
            covariance(idin,idin) = alpha(_-(nobs-1), _-(nobs-1));

            if (make_groups && !id_new.empty()) {
                // Ungroup grouped sources
                covariance = covariance(id_new,id_new);

                // TODO: what to do of grouped sources?
            }

            fits::update_table(out_file, ftable(covariance));
        }
    }

    if (has_groups) {
        // Post-process the fit result and identify inapropriate fits that are not grouped.
        // In particular: sources with significantly negative fluxes or with non-definite
        // errors, and their covarying neighbors.
        if (group_post_process) {
            vec1u ido = uindgen(nobs);
            vec1u idb = where(
                ((best_fit[ido] < 0 && (best_fit/best_fit_err)[ido] < -2) ||
                !is_finite(best_fit_err[ido]))
            );

            for (uint_t i : idb) {
                // Find a group to place this source and its neighbors in
                uint_t gid = npos;

                double bcov = group_cov_threshold;
                for (uint_t j : range(nobs)) {
                    double tcov = covar(i,j)/sqrt(covar(i,i)*covar(j,j));
                    if (tcov > bcov && is_grouped[j]) {
                        // We found one, but keep on going to make sure we pick the group
                        // that has the highest covariance
                        bcov = tcov;

                        if (group_fit_id[j] != npos) {
                            vec1u idg = where(old_cat.group_fit_id == group_fit_id[j]);
                            gid = old_cat.group_aper_id[idg[0]];
                        } else {
                            gid = old_cat.group_aper_id[id_old[j]];
                        }
                    }
                }

                bool new_group = gid == npos;
                if (new_group) {
                    // No group nearby, we have to create a new one...
                    gid = max(group_cat.id)+1;
                }

                // Notify sources of their new group
                vec1u nidg;
                for (uint_t j : range(nobs)) {
                    if (covar(i,j)/sqrt(covar(i,i)*covar(j,j)) > group_cov_threshold && !is_grouped[j]) {
                        if (group_fit_id[j] != npos) {
                            vec1u idg = where(old_cat.group_fit_id == group_fit_id[j]);
                            old_cat.group_aper_id[idg] = gid;
                            append(nidg, idg);
                        } else {
                            old_cat.group_aper_id[id_old[j]] = gid;
                            nidg.push_back(j);
                        }

                        is_grouped[j] = true;
                    }
                }

                if (new_group) {
                    group_cat.id.push_back(gid);
                    group_cat.ra.push_back(mean(old_cat.ra[nidg]));
                    group_cat.dec.push_back(mean(old_cat.dec[nidg]));

                    group_cat.nsrc.push_back(nidg.size());
                    group_cat.fit.push_back(false);
                    if (flux_prior) {
                        group_cat.fprior.push_back(0.0);
                        group_cat.fprior_err.push_back(0.0);
                    }
                } else {
                    uint_t iid = where(group_cat.id == gid)[0];
                    group_cat.nsrc[iid] += nidg.size();
                }
            }
        }

        fits::update_table(out_file,
            ftable(old_cat.group_fit_id, old_cat.group_aper_id)
        );

        // Produce a partial residual map where all non-grouped sources are removed
        auto tpg = progress_start(nobs);
        for (uint_t i : range(nobs)) {
            // Skip grouped sources
            if (is_grouped[i]) continue;

            // Subtract the rest
            vec2f tpsf = translate(psf, dy[i], dx[i]);
            vec1u idi, idp;
            subregion(img, {iy[i]-hsize, ix[i]-hsize, iy[i]+hsize, ix[i]+hsize}, idi, idp);

            img[idi] -= best_fit[i]*tpsf[idp];

            if (verbose) progress(tpg, 13);
        }

        // Remove the background, if not already done
        if (free_bg) {
            img -= background*fconv;
        }

        // Convert all sources to image coordinates
        vec1d tx, ty;
        fits::ad2xy(fits::extast(hdr), old_cat.ra, old_cat.dec, tx, ty);
        tx -= 1; ty -= 1;

        vec1i tix = round(tx), tiy = round(ty);
        vec1d tdx = tx - tix, tdy = ty - tiy;

        vec1f cov(old_cat.ra.size());

        vec2f aper = vec2f{circular_mask(psf.dims, group_aper_size) > 0.5};

        vec2i grp_map(img.dims);

        // For each group, measure an aperture flux within the area covered by
        // the grouped priors. The contribution of each prior to the total flux is then
        // studied outside of this program.

        // First build the masks of each group.
        // We want to make sure that the same pixels are not counted twice in different
        // groups, so we have to exclude the regions where two (or more) groups overlap.
        for (uint_t i : range(group_cat.ra)) {
            if (group_cat.fit[i]) {
                // This is a 'group_fit'
                // We don't need to care about it any more.
                continue;
            }

            // Locate the sources that are part of this group
            vec1u id = where(old_cat.group_aper_id == group_cat.id[i]);
            phypp_check(!id.empty(), "aper group ", group_cat.id[i], " is empty...");

            // Extract just what we need from the whole map
            uint_t xmi = max(0,             floor(min(tx[id]) - group_aper_size));
            uint_t xma = min(img.dims[1]-1, ceil(max(tx[id])  + group_aper_size));
            uint_t ymi = max(0,             floor(min(ty[id]) - group_aper_size));
            uint_t yma = min(img.dims[0]-1, ceil(max(ty[id])  + group_aper_size));

            vec2i tgrp = grp_map(ymi-_-yma,xmi-_-xma);

            // Convert coordinates to the local map
            tix[id] -= xmi;
            tiy[id] -= ymi;

            // Build the aperture mask
            vec2b mask(tgrp.dims);
            for (uint_t j : id) {
                // Create the aperture for this source
                vec2b taper = translate(aper, tdy[j], tdx[j]) > 0.5;
                vec1u idi, idp;
                subregion(mask, {tiy[j]-hsize, tix[j]-hsize, tiy[j]+hsize, tix[j]+hsize}, idi, idp);

                // Add it to the mask
                mask[idi] = mask[idi] || taper[idp];
            }

            // Flag pixels that already belong to another group
            vec1u ido = where(tgrp != 0 && mask);
            tgrp[ido] = -1;
            // Remove these pixels from the mask
            mask[ido] = false;
            if (count(mask) == 0) {
                warning("group ", group_cat.id[i], " has empty mask");
            }

            // Set pixels of this group
            tgrp[where(mask)] = group_cat.id[i];
            // Save back in the whole map
            grp_map(ymi-_-yma,xmi-_-xma) = tgrp;
        }

        // Now that the group map is built, compute the fraction of the flux of each prior
        // that is covered in the selected area and measure the aperture flux
        for (uint_t i : range(group_cat.ra)) {
            if (group_cat.fit[i]) {
                // This is a 'group_fit'
                // We don't need to care about it any more.
                continue;
            }

            // Locate the sources that are part of this group
            vec1u id = where(old_cat.group_aper_id == group_cat.id[i]);
            phypp_check(!id.empty(), "aper group ", group_cat.id[i], " is empty...");

            // Extract just what we need from the whole map
            uint_t xmi = max(0,             floor(min(tx[id]) - group_aper_size));
            uint_t xma = min(img.dims[1]-1, ceil(max(tx[id])  + group_aper_size));
            uint_t ymi = max(0,             floor(min(ty[id]) - group_aper_size));
            uint_t yma = min(img.dims[0]-1, ceil(max(ty[id])  + group_aper_size));

            vec2d timg = img(ymi-_-yma,xmi-_-xma);
            vec2d terr = err(ymi-_-yma,xmi-_-xma);
            vec2b mask = grp_map(ymi-_-yma,xmi-_-xma) == group_cat.id[i];

            if (count(mask) == 0) {
                warning("group ", group_cat.id[i], " has empty mask");
            }

            vec2d npsf = psf/total(psf);

            // Compute the fraction of flux contained in the mask for each source
            for (uint_t j : id) {
                // Create the model PSF for this source
                vec1f tpsf = flatten(translate(npsf, tdy[j], tdx[j]));
                vec1u idi, idp;
                subregion(mask, {tiy[j]-hsize, tix[j]-hsize, tiy[j]+hsize, tix[j]+hsize}, idi, idp);

                // Compute the coverage fraction
                vec1u idm = where(mask[idi]);
                cov[j] = total(tpsf[idp[idm]]);
            }

            // Compute flux
            vec1u idm = where(mask);
            group_cat.flux.push_back(total(timg[idm])/total(psf));

            // Error is just assuming Gaussian fluctations on each pixel
            group_cat.flux_err.push_back(sqrt(total(sqr(terr[idm])))/total(psf));
        }

        // Save the group map if requested
        if (!grp_file.empty()) {
            fits::write(grp_file, grp_map, hdr);
        }

        // Save to disk
        fits::update_table(out_file,
            "group_cov", cov, "group_id", group_cat.id,
            "group_ra", group_cat.ra, "group_dec", group_cat.dec,
            "group_flux", group_cat.flux, "group_flux_err", group_cat.flux_err,
            "group_fit", group_cat.fit, "group_nsrc", group_cat.nsrc
        );

        if (!res_file.empty()) {
            // Save the residual map
            // Convert to map units and add background
            img /= fconv;
            img += background;

            // Save to disk
            fits::write(res_file, img, hdr);
        }
    } else if (make_groups) {
        // No group was needed, just output empty maps & tables
        if (!grp_file.empty()) {
            vec2u grp_map(img.dims);
            fits::write(grp_file, grp_map, hdr);
        }

        fits::update_table(out_file,
            "group_cov", vec1f(old_cat.ra.size()), "group_id", vec1u(),
            "group_ra", vec1d(), "group_dec", vec1d(),
            "group_flux", vec1f(), "group_flux_err", vec1f(),
            "group_fit", vec1b(), "group_nsrc", vec1u()
        );
    }

    // If needed, produce the residual map
    if (!res_file.empty() && !has_groups) {
        if (verbose) {
            print("creating residual map");
        }

        // Remove the sources from the map.
        auto tpg = progress_start(nobs);
        for (uint_t i : range(nobs)) {
            // Subtract the source
            vec1f tpsf = flatten(translate(psf, dy[i], dx[i]));
            vec1u idi, idp;
            subregion(img, {iy[i]-hsize, ix[i]-hsize, iy[i]+hsize, ix[i]+hsize}, idi, idp);

            img[idi] -= best_fit[i]*tpsf[idp];

            if (verbose) progress(tpg, 13);
        }

        // Convert back to map units, and add back the background if it was subtracted
        img /= fconv;
        if (!free_bg) {
            img += fixed_bg;
        }

        // Save to disk
        fits::write(res_file, img, hdr);
    }

    if (!mod_file.empty()) {
        if (verbose) {
            print("creating model map");
        }

        vec2d mod = img*0;
        auto tpg = progress_start(nobs);
        for (uint_t i : range(nobs)) {
            vec1f tpsf = flatten(translate(psf, dy[i], dx[i]));
            vec1u idi, idp;
            subregion(mod, {iy[i]-hsize, ix[i]-hsize, iy[i]+hsize, ix[i]+hsize}, idi, idp);

            mod[idi] += best_fit[i]*tpsf[idp];

            if (verbose) progress(tpg, 13);
        }

        mod /= fconv;

        fits::write(mod_file, mod, hdr);
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("pixfit v1.0");
    paragraph("usage: pixfit cat=... img=... psf=... out=... [...]");

    header("List of mandatory command line options:");
    bullet("cat", "[string] path to the source catalog (FITS file)");
    bullet("img", "[string] path to the flux map (FITS file)");
    bullet("out", "[string] path to the output catalog (FITS file)");
    print("");

    header("List of optional command line options:");
    bullet("psf", "[string] path to the PSF (FITS file, default: guessed)");
    bullet("err", "[string] path to the error map (FITS file, default: guessed)");
    bullet("out_res", "[string] path to the output residual map (FITS file, default: none)");
    bullet("out_model", "[string] path to the output model map (FITS file, default: none)");
    bullet("fconv", "[float] map units to flux conversion factor (default: 1)");
    bullet("flux_prior", "[flag] use flux prior from the input catalog (default: no)");
    bullet("fixed_bg", "[float] fix the background to a given value (default: not fixed)");
    bullet("save_covariance", "[flag] also save the full covariance matrix in the output "
        "catalog (default: no)");
    bullet("cell_approx", "[flag] use a fitting approximation for non blended sources "
        "to make the fit significantly faster (default: no)");
    print("");
}
