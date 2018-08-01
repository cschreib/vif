#ifndef PHYPP_ASTRO_TEMPLATE_FIT_HPP
#define PHYPP_ASTRO_TEMPLATE_FIT_HPP

#include "phypp/astro/astro.hpp"
#include "phypp/math/mpfit.hpp"

namespace phypp {
namespace astro {
    // Convolve each SED with the response curve of the filters
    template<typename TLib, typename TFi>
    vec2d template_observed(const TLib& lib, const vec<1,TFi>& filters) {
        const uint_t nsed = lib.sed.dims[0];
        const uint_t nfilter = filters.size();

        vec2d flux(nsed, nfilter);
        for (uint_t f = 0; f < nfilter; ++f) {
            flux.safe(_,f) = sed2flux(filters[f], lib.lam, lib.sed);
        }

        return flux;
    }

    template<typename TLib, typename TFi>
    vec2d template_observed(TLib lib, double z, double d, const vec<1,TFi>& filters) {
        // Move each SED to the observed frame
        lib.sed = lsun2uJy(z, d, lib.lam, lib.sed);
        lib.lam *= (1.0 + z);

        // Convolve each SED with the response curve of the filters
        return template_observed(lib, filters);
    }

    template<typename TLib, typename TFi, typename TZ, typename TD>
    vec2d template_observed(const TLib& lib, const vec<1,TZ>& z, const vec<1,TD>& d,
        const vec<1,TFi>& filters) {

        phypp_check(z.size() == d.size(),
            "incompatible redshift and distance variables (", z.dims, " vs ", d.dims, ")");

        struct {
            vec2d lam, sed;
        } tlib;

        tlib.lam = lib.lam*(1.0 + mean(z));
        tlib.sed.resize(lib.sed.dims);

        // Combine the SEDs of each source to build an effective "redshift convolved" SED in the
        // observed frame
        for (uint_t i = 0; i < z.size(); ++i) {
            auto tlam = lib.lam*(1.0 + z[i]);
            auto tsed = lsun2uJy(z[i], d[i], lib.lam, lib.sed);
            tlib.sed += interpolate(tsed, tlam, tlib.lam);
        }

        tlib.sed /= z.size();

        // Convolve each SED with the response curve of the filters
        return template_observed(tlib, filters);
    }

    // Return the weight of an upper limit in a chi2 merit function.
    // Note that the analytical relation is numerically imprecise, and an asymptotic expression should
    // be used below -3. This function contains this approximation, which introduces at most an error of
    // 0.3%.
    // Assuming: d = (limit - model)/error
    // Also note that this function can also give the weight of a *lower* limit simply by using
    // d = (model - limit)/error.
    //
    // The underlying assumption is that, instead of a Gaussian weight
    //            exp(-(measure - true)^2/(2*error^2))
    // the upper limit weight is given by an error function
    //            0.5 + 0.5*erf(-(measure - limit)/(sqrt(2.0)*error))
    double limweight(double d) {
        return d < -3.0 ? d*d + 2.0*log(-2.0*sqrt(dpi/2.0)*d) : -2.0*log(0.5*(1.0 + erf(d/sqrt(2.0))));
    }

    template<std::size_t Dim, typename Type>
    vec<Dim,Type> limweight(vec<Dim,Type> d) {
        for (auto& x : d) {
            x = limweight(x);
        }

        return d;
    }

    struct template_fit_res_t {
        uint_t bfit; // index of the best fit template in the library
        vec1d chi2;  // chi^2 of each template
        vec1d amp;   // renormalization amplitude of each template
        vec2d flux;  // best-fit flux for each template

        vec1u sed_sim;       // index of each error realization's best fit template
        vec1d amp_bfit_sim;  // renormalization amplitude of the best fit for each error realization
        vec1d amp_sim;       // renormalization amplitude for each error realization's best fit
    };

    struct template_fit_params {
        uint_t nsim = 200;   // number of random realizations to perform to estimate errors
        bool renorm = false; // if true, allow templates to be renormalized when fitted
        bool ulim = false;   // if true, use upper limits to constrain the fit (negative errors)
        bool lib_obs = false; // if true, the input library is assumed to be in observer frame
    };

    // Note: Errors on the fit are computed by adding a random offset to the measured photometry (upper
    // limits are not touched) according to the provided error. The fit is performed on each of these
    // random realizations and the error on the parameters are computed as the standard deviation of the
    // fit results over all the realizations.
    //
    // Note: If 'ulim' is set to true, measurements with negative errors are considered as upper limits.
    // Since there is no analytical solution in this case, a numerical solver is used (mpfit). The
    // computation is thus slower (typically 3 times slower). To help the numerical solver, a classic
    // linear fit is performed only taking into account the measured values, and the best-fit amplitude
    // is used as a starting point for the solver.
    //
    // How and when to use upper limits:
    // The analysed source has been observed with an instrument, and the flux extraction procedure
    // produced a flux F and an error E on this flux. If the uncertainty on the measurement is really
    // gaussian with no systematic (to be tested with MC simulation with fake sources, then by analysing
    // the distribution of the extracted fluxes around the true value), then one should use the
    // extracted flux regardless on the signal to noise. Basically, saying that a source has a flux of
    // 1 mJy +/- 2 mJy (i.e. SNR = 0.5) provides more constraints than saying that is has an upper limit
    // of 6 mJy (3 sigma). Upper limits can still be used in this case, but they give less information
    // to the fitting procedure, and the best-fit will be more uncertain.
    // On the other hand, if the error distribution is NOT gaussian (typically because of confusion
    // noise, neighbor contamination or any other source of noise that may bias the result), it is safer
    // not to use any flux measurement below a given threshold L (to be determined from the MC
    // simulation). In this case, L should be considered as an upper limit, and the uncertainty on the
    // measurement E fed as the "error" on the limit (actually this is the error on the fact that the
    // flux is indeed below this limit, assuming gausssian statistics).
    template<typename TLib, typename TFi, typename TypeSeed,
        typename TZ, typename TD>
    template_fit_res_t template_fit(const TLib& lib, TypeSeed& seed, const TZ& z,
       const TD& d, vec1d flux, vec1d err,
       const vec<1,TFi>& filters, template_fit_params params = template_fit_params()) {

        template_fit_res_t res;
        if (params.lib_obs) {
            res.flux = template_observed(lib, filters);
        } else {
            res.flux = template_observed(lib, z, d, filters);
        }

        const uint_t nsed = res.flux.dims[0];
        const uint_t nfilter = filters.size();

        using ttype = decltype(err[0]*flux[0]*res.flux[0]);

        if (params.ulim) {
            vec1u idu = where(err < 0);
            vec1u idm = where(err > 0);
            err[idu] *= -1.0;
            flux /= err;

            for (uint_t i = 0; i < nsed; ++i) {
                res.flux(i,_) /= err;
            }

            if (params.renorm) {
                // Fit each template, compute chi2 and pick the best one
                res.chi2.resize(nsed);
                res.amp.resize(nsed);

                for (uint_t i = 0; i < nsed; ++i) {
                    vec1d model = res.flux(i,_);

                    auto lres = linfit(flux[idm], 1.0, model[idm]);
                    auto fres = mpfit([&](const vec1d& p) {
                        auto deviate = flux - p[0]*model;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, lres.params);

                    res.chi2[i] = fres.chi2;
                    res.amp[i] = fres.params[0];
                }

                // Find the best chi2 among all the SEDs
                res.bfit = min_id(res.chi2);

                // Compute the error with MC simulation
                res.amp_bfit_sim.resize(params.nsim);
                res.amp_sim.resize(params.nsim);
                res.sed_sim.resize(params.nsim);

                const uint_t nflux = flux.size();
                for (uint_t i = 0; i < params.nsim; ++i) {
                    auto fsim = flux;
                    fsim[idm] += randomn(seed, idm.size());

                    vec<1,ttype> amp(nsed), chi2(nsed);
                    for (uint_t i = 0; i < nsed; ++i) {
                        vec1d model = res.flux(i,_);
                        auto lres = linfit(fsim[idm], 1.0, model[idm]);
                        auto fres = mpfit([&](const vec1d& p) {
                            auto deviate = fsim - p[0]*model;
                            deviate[idu] = sqrt(limweight(deviate[idu]));
                            return deviate;
                        }, lres.params);

                        chi2[i] = fres.chi2;
                        amp[i] = fres.params[0];
                    }

                    auto ised = min_id(chi2);
                    res.sed_sim[i] = ised;
                    res.amp_sim[i] = amp[ised];
                    res.amp_bfit_sim[i] = amp[res.bfit];
                }
            } else {
                // Just compute chi2 and pick the best one
                res.chi2.resize(nsed);

                for (uint_t i = 0; i < nsed; ++i) {
                    auto deviate = flux - res.flux(i,_);
                    res.chi2[i] = total(sqr(deviate[idm])) + total(limweight(deviate[idu]));
                }

                // Find the best chi2 among all the SEDs
                res.bfit = min_id(res.chi2);

                // Then normalize the templates to fit the observation
                res.amp.resize(nsed);
                for (uint_t i = 0; i < nsed; ++i) {
                    vec1d model = res.flux(i,_);

                    auto lres = linfit(flux[idm], 1.0, model[idm]);
                    auto fres = mpfit([&](const vec1d& p) {
                        auto deviate = flux - p[0]*model;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, lres.params);

                    res.amp[i] = fres.params[0];
                }

                // Compute the error with MC simulation
                res.amp_bfit_sim.resize(params.nsim);
                res.amp_sim.resize(params.nsim);
                res.sed_sim.resize(params.nsim);

                const uint_t nflux = flux.size();
                for (uint_t i = 0; i < params.nsim; ++i) {
                    auto fsim = flux;
                    fsim[idm] += randomn(seed, idm.size());

                    vec<1,ttype> chi2(nsed);
                    for (uint_t i = 0; i < nsed; ++i) {
                        auto deviate = fsim - res.flux(i,_);
                        chi2[i] = total(sqr(deviate[idm])) + total(limweight(deviate[idu]));
                    }

                    auto ised = min_id(chi2);
                    res.sed_sim[i] = ised;

                    vec1d model = res.flux(ised,_);
                    auto lres = linfit(fsim[idm], 1.0, model[idm]);
                    auto fres = mpfit([&](const vec1d& p) {
                        auto deviate = fsim - p[0]*model;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, lres.params);

                    res.amp_sim[i] = fres.params[0];

                    model = res.flux(res.bfit,_);
                    lres = linfit(fsim[idm], 1.0, model[idm]);
                    fres = mpfit([&](const vec1d& p) {
                        auto deviate = fsim - p[0]*model;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, lres.params);

                    res.amp_bfit_sim[i] = fres.params[0];
                }
            }

            for (uint_t i = 0; i < nsed; ++i) {
                res.flux(i,_) *= err;
            }
        } else {
            // Compute chi2 & renormalization factor
            auto weight = invsqr(err);

            vec<1,ttype> tmp1(nsed), tmp2(nsed);
            for (uint_t i = 0; i < nsed; ++i) {
                auto tmp = res.flux(i,_);
                tmp1[i] = total(weight*flux*tmp);
                tmp2[i] = total(weight*tmp*tmp);
            }

            res.amp = tmp1/tmp2;

            res.chi2.resize(nsed);
            if (params.renorm) {
                for (uint_t i = 0; i < nsed; ++i) {
                    res.chi2.safe[i] = total(weight*sqr(flux - res.amp[i]*res.flux(i,_)));
                }
            } else {
                for (uint_t i = 0; i < nsed; ++i) {
                    res.chi2.safe[i] = total(weight*sqr(flux - res.flux(i,_)));
                }
            }

            // Find the best chi2 among all the SEDs
            res.bfit = min_id(res.chi2);

            // Compute the error with MC simulation
            res.amp_bfit_sim.resize(params.nsim);
            res.amp_sim.resize(params.nsim);
            res.sed_sim.resize(params.nsim);

            const uint_t nflux = flux.size();
            for (uint_t i = 0; i < params.nsim; ++i) {
                auto fsim = flux + randomn(seed, nflux)*err;
                for (uint_t t = 0; t < nsed; ++t) {
                    tmp1[t] = total(weight*fsim*res.flux(t,_));
                }

                auto amp = tmp1/tmp2;

                vec<1,ttype> chi2(nsed);
                if (params.renorm) {
                    for (uint_t t = 0; t < nsed; ++t) {
                        chi2.safe[t] = total(weight*sqr(fsim - amp[t]*res.flux(t,_)));
                    }
                } else {
                    for (uint_t t = 0; t < nsed; ++t) {
                        chi2.safe[t] = total(weight*sqr(fsim - res.flux(t,_)));
                    }
                }

                auto ised = min_id(chi2);

                res.amp_bfit_sim[i] = amp[res.bfit];
                res.amp_sim[i] = amp[ised];
                res.sed_sim[i] = ised;
            }
        }

        for (uint_t f = 0; f < nfilter; ++f) {
            res.flux(_,f) *= res.amp;
        }

        return res;
    }


    struct multi_template_fit_res_t {
        vec1u bfit; // index of the best fit template in each library
        vec1d amp;  // amplitude of the best fit of each library
        vec2d flux; // best-fit flux of each library
        double chi2;

        vec2u sed_sim; // index of each error realization's best fit templates
        vec2d amp_sim; // renormalization amplitude for each error realization's best fit
    };

    struct multi_template_fit_params {
        uint_t nsim = 200; // number of random realizations to perform to estimate errors
        bool ulim = false; // if true, use upper limits to constrain the fit (negative errors)
        bool lib_obs = false; // if true, the input library is assumed to be in observer frame
    };

    template<typename TLibs, typename TFi, typename TypeSeed,
        typename TZ, typename TD, typename TConst>
    multi_template_fit_res_t multi_template_fit(const TLibs& libs, TypeSeed& seed, const TZ& z,
       const TD& d, vec1d flux, vec1d err,
       const vec<1,TFi>& filters, multi_template_fit_params params = multi_template_fit_params(),
       TConst&& constraints = [](const vec1u&) { return true; }) {

        multi_template_fit_res_t res;

        uint_t nlib = libs.size();

        uint_t nseds = 1;
        vec1u nsed(nlib);
        std::vector<vec2f> convflux(nlib);
        for (uint_t l : range(nlib)) {
            nsed[l] = libs[l].sed.dims[0];
            nseds *= nsed[l];

            if (params.lib_obs) {
                convflux[l] = template_observed(libs[l], filters);
            } else {
                convflux[l] = template_observed(libs[l], z, d, filters);
            }
        }

        const uint_t nfilter = filters.size();

        using ttype = decltype(err[0]*flux[0]*res.flux[0]);

        if (params.ulim) {
            vec1u idu = where(err < 0);
            vec1u idm = where(err > 0);
            err[idu] *= -1.0;

            vec2f rflux(params.nsim, nfilter);
            vec1d rchi2 = replicate(dinf, params.nsim);
            res.sed_sim.resize(params.nsim, nlib);
            res.amp_sim.resize(params.nsim, nlib);

            for (uint_t s : range(params.nsim)) {
                rflux(s,_) = flux + err*randomn(seed, nfilter);
            }

            res.chi2 = dinf;

            vec1u ilib = replicate(0, nlib);
            for (uint_t i : range(nseds)) {
                if (!constraints(ilib)) continue;

                vec2f tpls(nlib, nfilter);
                for (uint_t l : range(nlib)) {
                    tpls(l,_) = convflux[l](ilib[l],_);
                }

                auto tres = linfit_pack(flux[idm], err[idm], tpls(_,idm));
                auto fres = mpfit([&](const vec1d& p) {
                    auto deviate = flux;
                    for (uint_t l : range(nlib)) {
                        deviate -= p[l]*tpls(l,_);
                    }

                    deviate /= err;
                    deviate[idu] = sqrt(limweight(deviate[idu]));
                    return deviate;
                }, tres.params);

                if (fres.chi2 < res.chi2) {
                    res.chi2 = fres.chi2;
                    res.bfit = ilib;
                    res.amp = fres.params;
                }

                for (uint_t s : range(params.nsim)) {
                    fres = mpfit([&](const vec1d& p) {
                        auto deviate = rflux(s,_);
                        for (uint_t l : range(nlib)) {
                            deviate -= p[l]*tpls(l,_);
                        }

                        deviate /= err;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, tres.params);

                    if (fres.chi2 < rchi2[s]) {
                        rchi2[s] = fres.chi2;
                        res.sed_sim(s,_) = ilib;
                        res.amp_sim(s,_) = fres.params;
                    }
                }

                increment_index_list(ilib, nsed);
            }
        } else {
            vec2f rflux(params.nsim, nfilter);
            vec1d rchi2 = replicate(dinf, params.nsim);
            res.sed_sim.resize(params.nsim, nlib);
            res.amp_sim.resize(params.nsim, nlib);

            flux /= err;

            for (uint_t s : range(params.nsim)) {
                rflux(s,_) = flux + randomn(seed, nfilter);
            }

            res.chi2 = dinf;

            for (uint_t l : range(nlib)) {
                for (uint_t i : range(nsed[l])) {
                    convflux[l](i,_) /= err;
                }
            }

            vec1u ilib = replicate(0, nlib);
            for (uint_t i : range(nseds)) {
                if (constraints(ilib)) {
                    matrix::mat2d alpha(nlib, nlib);
                    for (uint_t k1 : range(nlib))
                    for (uint_t k2 : range(k1, nlib)) {
                        alpha(k2,k1) = total(convflux[k1](ilib[k1],_)*convflux[k2](ilib[k2],_));
                    }

                    if (!inplace_invert_symmetric(alpha)) {
                        continue;
                    }

                    inplace_symmetrize(alpha);

                    vec1d beta(nlib);
                    for (uint_t k1 : range(nlib)) {
                        beta[k1] = total(flux*convflux[k1](ilib[k1],_));
                    }

                    vec1d amp = alpha*beta;

                    vec1d model(nfilter);
                    for (uint_t k1 : range(nlib)) {
                        model += amp[k1]*convflux[k1](ilib[k1],_);
                    }

                    double chi2 = total(sqr(flux - model));
                    if (chi2 < res.chi2) {
                        res.chi2 = chi2;
                        res.bfit = ilib;
                        res.amp = amp;
                    }

                    for (uint_t s : range(params.nsim)) {
                        for (uint_t k1 : range(nlib)) {
                            beta[k1] = total(rflux(s,_)*convflux[k1](ilib[k1],_));
                        }

                        amp = alpha*beta;

                        model[_] = 0;
                        for (uint_t k1 : range(nlib)) {
                            model += amp[k1]*convflux[k1](ilib[k1],_);
                        }

                        chi2 = total(sqr(rflux(s,_) - model));
                        if (chi2 < rchi2[s]) {
                            rchi2[s] = chi2;
                            res.sed_sim(s,_) = ilib;
                            res.amp_sim(s,_) = amp;
                        }
                    }
                }

                increment_index_list(ilib, nsed);
            }
        }

        res.flux.resize(nlib, nfilter);
        if (is_finite(res.chi2)) {
            for (uint_t l : range(nlib)) {
                res.flux(l,_) = res.amp[l]*convflux[l](res.bfit[l],_)*err;
            }
        } else {
            res.bfit = replicate(npos, nlib);
            res.amp = replicate(dnan, nlib);
            res.flux[_] = fnan;
        }

        return res;
    }
}
}

#endif
