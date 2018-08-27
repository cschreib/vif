#include <vif.hpp>
#include <vif/astro/template_fit.hpp>
#include "pixfit-common.hpp"

void print_help() {
    using namespace terminal_format;

    print("gfit v1.0");
    paragraph(
        "Perform SED fitting on the provided flux catalog with an infrared template library. "
        "If more than two 3-sigma measurements above 30um rest-frame are available, this "
        "procedure fits all the templates of the library by renormalizing them so each best-fit "
        "the data. Then, the best-fit among all the templates is chosen. Else, the behavior of the"
        "code depends on the chosen library."
    );

    header("The output file contains several quantities for each source:");
    bullet("nband", "total number of bands used in the fit");
    bullet("qflag", "-2: no fit, -1: no photometry above 30um observed, 0: no photometry above 30um "
        "rest-frame, 1: photometry above 30um rest-frame");
    bullet("chi2", "chi^2 of the best-fit");
    bullet("bfit_pah", "index of the best-fit template in the library");
    bullet("bfit_dust", "index of the best-fit template in the library");
    bullet("bfit_pah_err", "estimated error on this index");
    bullet("bfit_dust_err", "estimated error on this index");
    bullet("m_pah", "renormalization factor applied to the best-fit template");
    bullet("m_dust", "estimated error on the renormalization factor");
    bullet("flux", "best-fit flux for each fitted band");
    bullet("lir", "total IR luminosity estimated from the best-fit");
    bullet("lir_err", "estimated error on the total IR luminosity");
    print("");

    header("List of available command line options:");
    bullet("cat", "[string] path to the input catalog (FITS file)");
    bullet("out", "[string] path to the output catalog (FITS file) (default: [cat]_cefit.fits)");
    bullet("z", "[string] name of the column in the input catalog that contains the redshift");
    bullet("bands", "[string] POSIX regular expression to select the bands used in the fit");
    bullet("notes", "[string] POSIX regular expression to select the notes of the bands used in the "
        "fit");
    bullet("bands_notes", "[string] POSIX regular expression to select both bands and notes used in "
        "the fit");
    bullet("fix_fpah", "[flag] force fPAH to its redshift calibration (default: no)");
    bullet("fix_tdust", "[flag] force Tdust to its redshift calibration (default: no)");
    bullet("tdust_range", "[float,float] only allow Tdust values within the provided range (default: all)");
    bullet("ulim", "[flag] use or do not use upper limits on the photometry (default: 0)");
    bullet("cosmo", "[string] cosmological parameter set to use (wmap, std, ...) (default: std)");
    bullet("suffix", "[string] suffix to append to the output file (default: gfit)");
    bullet("thread", "[unsigned integer] number of threads to use to perform the fits (the more, "
        "the fastest, to some extent. Optimal value depends on the number of CPU your computer has)");
    bullet("verbose", "[flag] show a progress bar to estimate computing time");
}

int vif_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    // Default parameters
    std::string cat = "";
    std::string zvar = "z";
    std::string out = "";
    std::string scosmo = "std";
    std::string suffix = "gfit";
    std::string filter_db = data_dir+"fits/filter-db/db.dat";
    std::string irlib = data_dir+"fits/templates/s16-irlib/";
    uint_t thread = 1;
    uint_t tseed = 42;
    uint_t nsim = 100;
    double min_lam = 4.0;
    double snr_max = finf;
    double snr_max_group = finf;
    double max_fpah = dnan;
    bool fix_fpah = false;
    bool fix_tdust = false;
    bool no_negative = false;
    uint_t ntgrid = 10;
    uint_t ntrep = 2;
    vec1f tdust_range;
    bool verbose = false;
    vec1u only_ids;

    // Read parameters from command line
    read_args(argc, argv, arg_list(name(cat, "incat"), name(zvar, "z"), data_dir,
        verbose, thread, nsim, name(tseed, "seed"), name(out, "outcat"),
        name(scosmo, "cosmo"), suffix, snr_max, snr_max_group, ntgrid, ntrep,
        no_negative, filter_db, irlib, fix_tdust, fix_fpah, tdust_range,
        min_lam, max_fpah, only_ids
    ));

    irlib = file::directorize(irlib);

    auto seed = make_seed(tseed);

    if (!is_finite(max_fpah) || max_fpah > 1.0) {
        max_fpah = 1.0;
    } else if (max_fpah < 0) {
        max_fpah = 0.0;
    }

    if (!is_finite(snr_max_group)) {
        snr_max_group = snr_max;
    }

    if (cat.empty()) {
        error("gfit: missing flux catalog (cat=...)");
        print_help();
        return -1;
    }

    if (tdust_range.size() == 1 || tdust_range.size() > 2) {
        error("'tdust_range' must contain exactly two elements: [tdust_min, tdust_max]");
        return 1;
    }

    auto cosmo = get_cosmo(scosmo);

    // Bake output file name if not provided
    if (out.empty()) {
        vec1s spl = split(cat, ".");
        if (spl.size() > 1) {
            out = collapse(spl[_-(spl.size()-2)], ".")+"_"+suffix+".fits";
        } else {
            out = cat+"_"+suffix+".fits";
        }
    }

    if (verbose) {
        print("saving output in '"+out+"'");
    }

    // Read the source catalog
    struct {
        vec1u id;
        vec1d ra, dec;
        vec1d z;
        vec2d flux, flux_err, flux_group_cov;
        vec2u flux_group;
        vec1d group_flux, group_flux_err;
        vec1s bands;
        vec1d tdust_min, tdust_max;
    } fcat;

    fits::read_table_loose(cat, ftable(
        fcat.id, fcat.ra, fcat.dec, fcat.z, fcat.bands, fcat.flux, fcat.flux_err,
        fcat.flux_group_cov, fcat.flux_group, fcat.group_flux, fcat.group_flux_err,
        fcat.tdust_min, fcat.tdust_max
    ));

    if (fcat.z.empty()) {
        if (zvar.empty()) {
            error("gfit: no redshift information in '"+cat+"' (Z)");
            return 1;
        }

        fits::read_table(cat, zvar, fcat.z);
        if (fcat.z.empty()) {
            error("gfit: no redshift information in '"+cat+"' ("+to_upper(zvar)+")");
            return 1;
        }
    }

    vif_check(fcat.z.size() == fcat.flux.dims[0], "number of redshifts (", fcat.z.size(), ") does "
        "not match number of flux (", fcat.flux.dims, ")");

    if (fcat.bands.empty()) {
        error("gfit: no band to fit (missing 'BANDS')");
        return 1;
    } else {
        vif_check(fcat.bands.size() == fcat.flux.dims[1], "number of bands (", fcat.bands.size(),
            ") does not match number of flux (", fcat.flux.dims, ")");
    }

    if (fcat.flux_err.empty()) {
        warning("could not find any flux uncertainty, assuming SNR=5 for all bands!");
        fcat.flux_err = fcat.flux/5;
    }

    // Load filter response curves
    auto fdb = read_filter_db(filter_db);
    filter_bank_t fbank;
    if (!get_filters(fdb, fcat.bands, fbank)) {
        return 1;
    }

    uint_t nband = fbank.size();

    vec1f lam(nband);
    for (uint_t f = 0; f < nband; ++f) {
        lam[f] = fbank[f].rlam;
    }

    // Load template library
    struct lib_t {
        vec2d lam, sed;
        vec1d lir;
        vec1d umin;
        vec1d umean;
        vec1d l8;
        vec1d tdust;
    };

    auto read_lib = [](std::string libfile, lib_t& lib) {
        fits::read_table(libfile, ftable(
            lib.lam, lib.sed, lib.lir, lib.umin, lib.umean, lib.l8, lib.tdust
        ));
    };

    std::vector<lib_t> libs(2);
    read_lib(irlib+"s16_dust.fits", libs[0]);
    read_lib(irlib+"s16_pah.fits", libs[1]);
    const uint_t nlib = libs.size();

    // Apply constraints on Tdust
    vec1u itdust;
    if (!tdust_range.empty()) {
        uint_t onsed = libs[0].tdust.size();
        auto b = bounds(libs[0].tdust, tdust_range[0], tdust_range[1]);
        if (b[0] == npos) {
            warning("the IR library only goes as low as Tdust=", min(libs[0].tdust));
            b[0] = 0;
        }
        if (b[1] == npos) {
            warning("the IR library only goes as high as Tdust=", max(libs[0].tdust));
            b[1] = onsed-1;
        }

        itdust = uindgen(b[1]-b[0]+1) + b[0];

        for (uint_t l : range(2)) {
            libs[l].lam = libs[l].lam(itdust,_);
            libs[l].sed = libs[l].sed(itdust,_);
            libs[l].lir = libs[l].lir[itdust];
            libs[l].umin = libs[l].umin[itdust];
            libs[l].umean = libs[l].umean[itdust];
            libs[l].l8 = libs[l].l8[itdust];
            libs[l].tdust = libs[l].tdust[itdust];
        }

        if (verbose) {
            print("reducing IR library to Tdust=", min(libs[0].tdust),
                " to Tdust=", max(libs[0].tdust));
            print("using ", libs[0].tdust.size(), " SEDs out of ", onsed);
        }
    } else {
        itdust = uindgen(libs[0].tdust.size());
    }

    // Handle Tdust constraints
    if (fcat.tdust_min.empty()) {
        fcat.tdust_min = replicate(fnan, fcat.z.dims);
    }
    if (fcat.tdust_max.empty()) {
        fcat.tdust_max = replicate(fnan, fcat.z.dims);
    }

    for (uint_t i : range(fcat.z)) {
        if (fcat.tdust_min[i] < 0 or !is_finite(fcat.tdust_min[i])) {
            fcat.tdust_min[i] = libs[0].tdust.front();
        }

        if (fcat.tdust_max[i] < 0 or !is_finite(fcat.tdust_max[i])) {
            fcat.tdust_max[i] = libs[0].tdust.back();
        }
    }

    const uint_t nsed = libs[0].tdust.size();

    const uint_t ngal = fcat.flux.dims[0];
    if (verbose) {
        print(ngal, " sources in the catalog");
    }

    // Shortcuts...
    vec1f z = fcat.z;
    vec2f flux = fcat.flux;
    vec2f flux_err = fcat.flux_err;
    vec2u flux_group = fcat.flux_group;
    vec2f flux_cov = fcat.flux_group_cov;
    vec1f group_flux = fcat.group_flux;
    vec1f group_err = fcat.group_flux_err;
    vec1f tdustl = fcat.tdust_min;
    vec1f tdustu = fcat.tdust_max;

    // Adjust input catalog
    for (uint_t i : range(ngal)) {
        // Remove sources which will not contribute to their flux groups
        vec1b isg = flux_group(i,_) != npos;
        vec1b bad = lam/(1.0 + z[i]) < min_lam || flux_cov(i,_) == 0.0;
        vec1u idb = where(isg && bad);
        flux_group(i,idb) = npos;

        // Apply maximal SNR
        if (is_finite(snr_max)) {
            vec1u idm = where(flux(i,_)/flux_err(i,_) > snr_max && flux_err(i,_) > 0);
            flux_err(i,idm) = flux(i,idm)/snr_max;

            vec1u idg = flux_group(i, where(flux_group(i,_) != npos));
            idm = where(group_flux[idg]/group_err[idg] > snr_max_group);
            group_err[idg[idm]] = group_flux[idg[idm]]/snr_max_group;
        }

        // Also give flux coverage of 1.0 to single measurements
        vec1u idz = where(flux_cov(i,_) == 0.0 || !is_finite(flux_cov(i,_)));
        flux_cov(i,idz) = 1.0;
    }

    // Pick sources that are detected in the IR
    vec1u gifit; vec2b gff; vec2b ggf; vec2b ggg; {
        vec2b rlam = replicate(lam, ngal)/transpose(replicate(1.0 + z, nband)) >= min_lam;
        vec2b mindiv = is_finite(flux) && is_finite(flux_err) && flux_err > 0;
        if (no_negative) mindiv = mindiv && flux > 0;
        vec2b mgroup = flux_group != npos;

        gff = mindiv && rlam;
        ggf = (mindiv || mgroup) && rlam;
        ggg = mgroup && rlam;

        gifit = where(partial_count(1, gff) > 0);
    }

    const uint_t nfit = gifit.size();
    const uint_t ngroup = group_flux.size();
    if (verbose) {
        print(nfit, " sources to fit");
    }

    // Apply selection
    flux = flux(gifit,_);
    flux_err = flux_err(gifit,_);
    flux_group = flux_group(gifit,_);
    tdustl = tdustl[gifit];
    tdustu = tdustu[gifit];
    gff = gff(gifit,_);
    ggf = ggf(gifit,_);
    ggg = ggg(gifit,_);
    z = z[gifit];

    // Compute distance of each galaxy
    vec1f d = lumdist(z, cosmo);

    // Convolve each redshifted SED of the library with the filters.
    // Build a table of model fluxes in uJy/Mdust.
    // TODO: apply measure /= error here
    vec3d convd(nfit,nsed,nband);
    vec3d convp(nfit,nsed,nband);
    for (uint_t i : range(nfit)) {
        // Compute broad band fluxes
        convd(i,_,_) = template_observed(libs[0], z[i], d[i], fbank);
        convp(i,_,_) = template_observed(libs[1], z[i], d[i], fbank);

        // Apply aperture coverage
        for (uint_t s : range(nsed)) {
            convd(i,s,_) *= flux_cov(gifit[i],_);
            convp(i,s,_) *= flux_cov(gifit[i],_);
        }
    }

    // Build connected groups to only fit connected sources simultaneously
    struct fit_group {
        uint_t id;
        vec1u sids;
        vec1f z;
        vec1f measures, errors;
        vec2u measure_id;
        vec3d convd, convp;
    };

    std::vector<fit_group> fgroups;

    vec1b done(nfit);
    for (uint_t i : range(nfit)) {
        if (done[i]) continue;

        vec1u idm = where(ggg(i,_));
        if (idm.empty()) {
            // This source is not grouped and can be fit separately.
            // Just create a new group with it alone inside.
            fit_group f;
            f.id = fgroups.size();
            f.sids.push_back(i);
            f.z.push_back(z[i]);
            done[i] = true;

            f.measure_id = replicate(npos, 1, nband);
            f.convd.resize(1, nsed, nband);
            f.convp.resize(1, nsed, nband);

            // Add single measurements
            idm = where(gff(i,_));
            f.measures = flux(i,idm);
            f.errors = flux_err(i,idm);
            f.measure_id(0,idm) = uindgen(idm.size());
            f.convd(0,_,_) = convd(i,_,_);
            f.convp(0,_,_) = convp(i,_,_);

            fgroups.push_back(f);
        } else {
            // This source is grouped and has to be fit together with its neighbors.

            // Create a new group
            fit_group f;
            f.id = fgroups.size();

            // Add this source
            f.sids.push_back(i);
            f.z.push_back(z[i]);
            done[i] = true;

            // Add all connected sources
            for (uint_t j : range(nfit)) {
                if (done[j]) continue;
                for (uint_t b : range(nband)) {
                    if (flux_group(i,b) == flux_group(j,b) && flux_group(i,b) != npos) {
                        done[j] = true;
                        f.sids.push_back(j);
                        f.z.push_back(z[j]);
                        break;
                    }
                }
            }

            const uint_t tnfit = f.sids.size();
            f.measure_id = replicate(npos, tnfit, nband);

            f.convd = convd(f.sids,_,_);
            f.convp = convp(f.sids,_,_);

            // Add single measurements
            for (uint_t j : range(f.sids)) {
                uint_t k = f.sids[j];
                idm = where(gff(k,_));
                uint_t i0 = f.measures.size();
                append(f.measures, flux(k,idm));
                append(f.errors, flux_err(k,idm));
                f.measure_id(j,idm) = uindgen(idm.size()) + i0;
            }

            // Then add groups
            vec1u group_id = replicate(npos, ngroup);
            for (uint_t j : range(f.sids)) {
                uint_t k = f.sids[j];
                for (uint_t b : range(nband)) {
                    if (!ggg(k,b)) continue;

                    uint_t idg = flux_group(k,b);

                    if (group_id[idg] == npos) {
                        uint_t i0 = f.measures.size();
                        f.measures.push_back(group_flux[idg]);
                        f.errors.push_back(group_err[idg]);
                        group_id[idg] = i0;
                    }

                    f.measure_id(j,b) = group_id[idg];
                }
            }

            fgroups.push_back(f);
        }
    }

    if (verbose) {
        print("found ", fgroups.size(), " independent groups:");
        uint_t nms = 0;
        for (auto& f : fgroups) {
            nms = std::max(nms, f.sids.size());
            print(" - ", f.sids.size(), " sources, ", f.measures.size(), " measurements");
        }
    }

    auto dofit = [&](const fit_group& f, vec1d tdust, vec1d& mdust, vec1d& mdust_err,
        vec1d& fpah, vec1d& fpah_err, vec1d& residuals, vec1b fixpah) {

        const uint_t tnfit = f.sids.size();
        const uint_t nparam = 2;
        const uint_t np = nparam*tnfit;
        const uint_t nm = f.measures.size();

        // Interpolate library to the requested Tdust
        vec2d tconvd(tnfit,nband);
        vec2d tconvp(tnfit,nband);

        for (uint_t i : range(tnfit)) {
            auto ised = bounds(libs[0].tdust, tdust[i]);
            double x = 0.0;
            if (ised[0] == npos) {
                ised[0] = ised[1];
            } else if (ised[1] == npos) {
                ised[1] = ised[0];
            } else {
                x = (tdust[i] - libs[0].tdust[ised[0]]) /
                    (libs[0].tdust[ised[1]] - libs[0].tdust[ised[0]]);
            }

            tconvd(i,_) = (1.0-x)*f.convd(i,ised[0],_) + x*f.convd(i,ised[1],_);
            if (fixpah[i]) {
                vec1f cp = (1.0-x)*f.convp(i,ised[0],_) + x*f.convp(i,ised[1],_);
                tconvd(i,_) = (1.0 - fpah[i])*tconvd(i,_) + fpah[i]*cp;
            } else {
                tconvp(i,_) = (1.0-x)*f.convp(i,ised[0],_) + x*f.convp(i,ised[1],_);
            }
        }

        // Build matrix
        matrix::mat2d alpha(np,np);
        vec1d beta(np);

        auto tmp = f.measures/f.errors;
        for (uint_t i : range(tnfit))
        for (uint_t j : range(i, tnfit)) {
            for (uint_t b : range(nband)) {
                // Only consider overlapping measurements
                if (f.measure_id(i,b) != f.measure_id(j,b) || f.measure_id(i,b) == npos) continue;
                uint_t m = f.measure_id(i,b);
                double weight = 1.0/f.errors[m];
                double weight2 = sqr(weight);

                // alpha(i,j) = sum over all points of x[i]*x[j]/e^2

                // Dust component is always here
                alpha(nparam*i+0,nparam*j+0) += tconvd(i,b)*tconvd(j,b)*weight2;

                // PAH may not be there
                if (!fixpah[i]) {
                    alpha(nparam*i+1,nparam*j+0) += tconvp(i,b)*tconvd(j,b)*weight2;
                    if (!fixpah[j]) {
                        alpha(nparam*i+0,nparam*j+1) += tconvd(i,b)*tconvp(j,b)*weight2;
                        alpha(nparam*i+1,nparam*j+1) += tconvp(i,b)*tconvp(j,b)*weight2;
                    }
                } else if (!fixpah[j]) {
                    alpha(nparam*i+0,nparam*j+1) += tconvd(i,b)*tconvp(j,b)*weight2;
                }

                if (i == j) {
                    // beta[i] = sum over all points of x[i]*y/e^2
                    beta[nparam*i+0] += tconvd(i,b)*tmp[m]*weight;
                    if (!fixpah[i]) {
                        beta[nparam*i+1] += tconvp(i,b)*tmp[m]*weight;
                    } else {
                        // f_PAH is frozen for this source, put a dummy value in the
                        // matrix only in the diagonal so that this parameter does not
                        // influence the fit
                        alpha(nparam*i+1, nparam*i+1) = 1.0;
                    }
                }
            }
        }

        inplace_symmetrize(alpha);

        // Invert
        if (!inplace_invert_symmetric(alpha)) {
            return dnan;
        }

        // Extract fit parameters and errors
        inplace_symmetrize(alpha);

        vec2d bfit = reform(alpha*beta, tnfit, nparam);
        vec2d berr = reform(sqrt(diagonal(alpha)), tnfit, nparam);

        vec1u idf = where(fixpah);
        mdust[idf] = bfit(idf,0);
        mdust_err[idf] = berr(idf,0);
        fpah_err[idf] = 0;

        vec1u idnf = where(!fixpah);
        mdust[idnf] = bfit(idnf,0) + bfit(idnf,1);
        fpah[idnf] = bfit(idnf,1)/mdust[idnf];
        mdust_err[idnf] = sqrt(sqr(berr(idnf,0)) + sqr(berr(idnf,1)));
        fpah_err[idnf] = sqrt(sqr(berr(idnf,1)/mdust[idnf]) + sqr((mdust_err[idnf]/mdust[idnf])*fpah[idnf]));

        vec1d model(nm);
        for (uint_t i : range(tnfit)) {
            for (uint_t b : range(nband)) {
                uint_t m = f.measure_id(i,b);
                if (m != npos) {
                    if (fixpah[i]) {
                        model[m] += mdust[i]*tconvd(i,b);
                    } else {
                        model[m] += bfit(i,0)*tconvd(i,b) + bfit(i,1)*tconvp(i,b);
                    }
                }
            }
        }

        // Return the chi2
        model /= f.errors;
        residuals = model - tmp;

        return total(sqr(residuals));
    };


    // Initialize results
    struct {
        vec1u group;
        vec1i lir_qflag, l8_qflag, mdust_qflag, tdust_qflag;
        vec1b fixed_tdust, fixed_fpah;
        vec2u bfit;
        vec1f chi2;
        vec2f chi2_tdust_y;
        vec2f chi2_tdust_x;
        vec1f lir, lir_err;
        vec1f l8, l8_err;
        vec1f mdust, mdust_err;
        vec1f fpah, fpah_err;
        vec1f tdust, tdust_low, tdust_up;

        vec2f flux;
        vec1s bands;
        vec1f lambda;
    } res;

    res.group = replicate(npos, ngal);

    res.chi2 = replicate(fnan, ngal);
    res.chi2_tdust_x = replicate(fnan, ngal, ntrep*ntgrid);
    res.chi2_tdust_y = replicate(fnan, ngal, ntrep*ntgrid);

    res.lir_qflag = replicate(-1, ngal);
    res.mdust_qflag = replicate(-1, ngal);
    res.l8_qflag = replicate(-1, ngal);
    res.tdust_qflag = replicate(-1, ngal);

    res.fixed_tdust = replicate(fix_tdust, ngal);
    res.fixed_fpah = replicate(fix_fpah, ngal);

    res.bfit = replicate(npos, nlib, ngal);

    res.lir = res.lir_err = res.l8 = res.l8_err = res.mdust = res.mdust_err =
        res.fpah = res.fpah_err = res.tdust = res.tdust_low = res.tdust_up = replicate(fnan, ngal);

    res.flux = replicate(fnan, ngal, fcat.bands.size());
    res.bands = fcat.bands;
    res.lambda = lam;

    auto fpah_z = vectorize_lambda([](double tz) {
        return 0.035 + 0.035*(1.0-0.8*clamp(tz, 1.0, 2.0));
    });

    auto tdust_z = vectorize_lambda([](double tz) {
        return tz < 2 ? 24.6*pow(1.0+tz, 0.37) : 29.7*pow(1.0+tz, 0.2);
    });

    for (auto& f : fgroups) {
        vec1u gids = gifit[f.sids];

        if (!only_ids.empty()) {
            // Only fit the groups that contain the sources we are interested in
            if (count(is_any_of(fcat.id[gids], only_ids)) == 0) continue;
        }

        const uint_t tnfit = f.sids.size();
        if (verbose) {
            print("fitting simultaneously ", tnfit, " sources");
        }

        vec1u tlib_id = replicate(0, tnfit);

        vec1u lib_id;
        vec1d tdust;
        vec1d mdust, mdust_err;
        vec1d fpah, fpah_err;
        vec1b fixed_fpah;
        vec1d residuals;

        double bchi2 = dinf;

        vec1d tmdust(tnfit);
        vec1d tmdust_err(tnfit);
        vec1d tfpah(tnfit);
        vec1d tfpah_err(tnfit);
        vec1d tresiduals(tnfit);

        if (fix_fpah) {
            tfpah = fpah_z(f.z);
        }

        auto dofit_wrap = [&](vec1d ttdust) {
            double chi2;
            vec1b ffpah = replicate(fix_fpah, tnfit);

            bool fitok = true;

            do {
                fitok = true;

                chi2 = dofit(f, ttdust, tmdust, tmdust_err, tfpah, tfpah_err,
                    tresiduals, ffpah);

                // Check the fit parameters for invalid values
                if (!fix_fpah) {
                    vec1u idbf = where(tfpah < 0 || tfpah > max_fpah);
                    if (!idbf.empty()) {
                        // At least one f_PAH has gone wild, redo the fit by keeping
                        // it frozen to the average value at that redshift
                        // tfpah[idbf] = fpah_z(f.z[idbf]);
                        tfpah[idbf] = 0.0;
                        ffpah[idbf] = true;
                        fitok = false;
                    }
                }
            } while (!fitok);

            // Save fit if the chi2 has improved
            if (chi2 < bchi2) {
                bchi2 = chi2;

                lib_id = tlib_id;
                tdust = ttdust;
                mdust = tmdust;
                mdust_err = tmdust_err;
                fpah = tfpah;
                fpah_err = tfpah_err;
                fixed_fpah = ffpah;
                residuals = tresiduals;
            }

            return chi2;
        };

        if (fix_tdust) {
            vec1d ttdust = tdust_z(f.z);

            dofit_wrap(ttdust);

            if (verbose) {
                print("best chi2: ", bchi2, " (reduced: ",
                    bchi2/(f.measures.size() - count(!fixed_fpah) - tnfit), ")");
            }

            res.tdust_low[gids] = tdust;
            res.tdust_up[gids] = tdust;
        } else {
            double t0 = min(libs[0].tdust);
            double t1 = max(libs[0].tdust);

            vec2d tdust_grid = replicate(rgen(t0, t1, ntgrid), tnfit);

            vec3d chi2_tdust_x = replicate(finf, tnfit, ntrep, ntgrid);
            vec3d chi2_tdust_y = replicate(finf, tnfit, ntrep, ntgrid);

            chi2_tdust_x(_,0,_) = tdust_grid;

            uint_t niter = pow(ntgrid, tnfit);

            auto doloop = [&](uint_t r) {
                auto pg = progress_start(niter);
                for (uint_t i : range(niter)) {
                    vec1d ttdust = diagonal(matrix::as_matrix(tdust_grid(_,tlib_id)));

                    // Check this does not violate Tdust constraints
                    bool bad = false;
                    for (uint_t j : range(gids)) {
                        if (ttdust[j] < tdustl[gids[j]] || ttdust[j] > tdustu[gids[j]]) {
                            bad = true;
                        }
                    }

                    if (!bad) {
                        double chi2 = dofit_wrap(ttdust);

                        // Also save the chi2 if this is the best for a given Tdust
                        for (uint_t j : range(tnfit)) {
                            if (chi2 < chi2_tdust_y(j,r,tlib_id[j])) {
                                chi2_tdust_y(j,r,tlib_id[j]) = chi2;
                            }
                        }
                    }

                    increment_index_list(tlib_id, ntgrid);
                    if (verbose) progress(pg, 117);
                }

                chi2_tdust_y(_,r,_) /= bchi2;

                if (verbose) {
                    print("best chi2: ", bchi2, " (reduced: ",
                        bchi2/(f.measures.size() - count(!fixed_fpah) - tnfit*2), ")");
                }
            };

            // Do a first fit on a coarse Tdust grid covering the whole parameter space
            // to locate the region containing the global minimum.
            doloop(0);

            // Then repeat the fit with a reduced Tdust grid covering only the regions
            // close to the global minimum, but with an increasingly finer grid
            for (uint_t r : range(ntrep-1)) {
                for (uint_t i : range(tnfit)) {
                    uint_t i0 = lib_id[i] > 0     ? lib_id[i]-1 : 0;
                    uint_t i1 = lib_id[i] < ntgrid-1 ? lib_id[i]+1 : ntgrid-1;
                    tdust_grid(i,_) = rgen(tdust_grid(i,i0), tdust_grid(i,i1), ntgrid);
                }

                chi2_tdust_x(_,r+1,_) = tdust_grid;
                doloop(r+1);
            }

            // Save chi2(Tdust) and use that to estimate uncertainty on Tdust
            chi2_tdust_y *= bchi2;
            chi2_tdust_y /= (f.measures.size() - count(!fixed_fpah) - tnfit*2);
            for (uint_t i : range(tnfit)) {
                vec1d tx = flatten(chi2_tdust_x(i,_,_));
                vec1d ty = flatten(chi2_tdust_y(i,_,_));
                vec1u sid = sort(tx);
                tx = tx[sid];
                ty = ty[sid];

                res.chi2_tdust_x(gids[i],_) = tx;
                res.chi2_tdust_y(gids[i],_) = ty;

                // Estimate probability from chi2
                sid = unique_ids_sorted(tx);
                tx = tx[sid]; ty = ty[sid];
                ty = exp(-(ty - bchi2));

                // Find Tdust so that probability goes lower than 1-sigma (exp(-1))
                double sigma = exp(-1.0);
                uint_t k = max_id(ty);
                double low = find_zero(reverse(tx[_-k]), reverse(ty[_-k]) - sigma);
                double up = find_zero(tx[k-_], ty[k-_] - sigma);

                if (!is_finite(up)) {
                    up = tx.back();
                }
                if (!is_finite(low)) {
                    low = tx.front();
                }

                res.tdust_low[gids[i]] = low;
                res.tdust_up[gids[i]] = up;
            }
        }

        vec1u bfit = clamp(round(interpolate(findgen(nsed), libs[0].tdust, tdust)), 0, nsed-1);

        res.group[gids] = f.id;
        res.bfit(0,gids) = itdust[bfit];
        res.bfit(1,gids) = itdust[bfit];
        res.tdust[gids] = tdust;
        res.mdust[gids] = mdust;
        res.mdust_err[gids] = mdust_err;
        res.fpah[gids] = fpah;
        res.fpah_err[gids] = fpah_err;
        res.fixed_fpah[gids] = fixed_fpah;

        // NOTE: ignores uncertainty on Tdust & fPAH!
        // TODO: use the covariance matrix to get the right uncertainty
        res.lir[gids] = mdust*((1.0 - fpah)*libs[0].lir[bfit] + fpah*libs[1].lir[bfit]);
        res.lir_err[gids] = mdust_err*((1.0 - fpah)*libs[0].lir[bfit] + fpah*libs[1].lir[bfit]);
        res.l8[gids] = mdust*((1.0 - fpah)*libs[0].l8[bfit] + fpah*libs[1].l8[bfit]);
        res.l8_err[gids] = mdust_err*((1.0 - fpah)*libs[0].l8[bfit] + fpah*libs[1].l8[bfit]);

        // Compute individual chi2 and fluxes
        for (uint_t i : range(tnfit)) {
            uint_t j = gifit[f.sids[i]];

            vec1u idm = where(f.measure_id(i,_) != npos);
            res.chi2[j] = mean(sqr(residuals[f.measure_id(i,idm)]));

            res.flux(j,_) = mdust[i]*(
                (1.0 - fpah[i])*f.convd(i,bfit[i],_) + fpah[i]*f.convp(i,bfit[i],_)
            );
            res.flux(j,_) /= flux_cov(j,_);
        }

        // Save the result
        fits::write_table(out, ftable(
            fcat.id, fcat.ra, fcat.dec,
            res.group, res.lir_qflag, res.l8_qflag, res.mdust_qflag, res.tdust_qflag,
            res.fixed_tdust, res.fixed_fpah, res.bfit, res.chi2, res.chi2_tdust_y,
            res.chi2_tdust_x, res.lir, res.lir_err, res.l8, res.l8_err, res.mdust, res.mdust_err,
            res.fpah, res.fpah_err, res.tdust, res.tdust_low, res.tdust_up, res.flux,
            res.bands, res.lambda
        ));
    }

    return 0;
}
