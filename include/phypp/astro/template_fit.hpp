#ifndef ASTRO_TEMPLATE_FIT_HPP
#define ASTRO_TEMPLATE_FIT_HPP

#include "phypp/astro.hpp"
#include "phypp/mpfit.hpp"

template<typename TypeLib, typename TypeFi>
vec2d template_observed(const TypeLib& lib, double z, double d,
    const vec_t<1,TypeFi>& filters) {

    // Move each SED to the observed frame
    vec2d rflam = lib.lam*(1.0 + z);
    vec2d rfsed = lsun2uJy(z, d, lib.lam, lib.sed);

    // Convolve each SED with the response curve of the filters
    const uint_t nsed = rfsed.dims[0];
    const uint_t nfilter = filters.size();
    vec2d flux = dblarr(nsed, nfilter);
    for (uint_t f = 0; f < nfilter; ++f) {
        flux(_,f) = sed2flux(filters[f], rflam, rfsed);
    }

    return flux;
}

template<typename TypeLib, typename TypeFi, typename TypeZ, typename TypeD>
vec2d template_observed(const TypeLib& lib, const vec_t<1,TypeZ>& z, const vec_t<1,TypeD>& d,
    const vec_t<1,TypeFi>& filters) {

    phypp_check(z.size() == d.size(),
        "incompatible redshift and distance variables ("+strn(z.dims)+" vs "+strn(d.dims)+")");

    vec2d rflam = lib.lam*(1.0 + mean(z));
    vec2d rfsed = dblarr(lib.sed.dims);

    // Combine the SEDs of each source to build an effective "redshift convolved" SED in the
    // observed frame
    for (uint_t i = 0; i < z.size(); ++i) {
        auto tlam = lib.lam*(1.0 + z[i]);
        auto tsed = lsun2uJy(z[i], d[i], lib.lam, lib.sed);
        rfsed += interpolate(tsed, tlam, rflam);
    }

    rfsed /= z.size();

    // Convolve each SED with the response curve of the filters
    const uint_t nsed = rfsed.dims[0];
    const uint_t nfilter = filters.size();
    vec2d flux = dblarr(nsed, nfilter);
    for (uint_t f = 0; f < nfilter; ++f) {
        flux(_,f) = sed2flux(filters[f], rflam, rfsed);
    }

    return flux;
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
vec_t<Dim,Type> limweight(vec_t<Dim,Type> d) {
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
template<typename TypeLib, typename TypeFi, typename TypeSeed,
    typename TypeZ, typename TypeD>
template_fit_res_t template_fit(const TypeLib& lib, TypeSeed& seed, const TypeZ& z,
   const TypeD& d, vec1d flux, vec1d err,
   const vec_t<1,TypeFi>& filters, template_fit_params params = template_fit_params()) {

    template_fit_res_t res;
    res.flux = template_observed(lib, z, d, filters);

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

                res.chi2[i] = fres.chi2;
                res.amp[i] = fres.params[0];
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

                vec_t<1,ttype> amp(nsed), chi2(nsed);
                for (uint_t i = 0; i < nsed; ++i) {
                    vec1d model = res.flux(i,_);
                    auto lres = linfit(fsim[idm], 1.0, model[idm]);
                    auto fres = mpfit([&](const vec1d& p) {
                        auto deviate = fsim - p[0]*model;
                        deviate[idu] = sqrt(limweight(deviate[idu]));
                        return deviate;
                    }, lres.params);

                    chi2[i] = fres.chi2;
                    amp[i] = fres.params[0];
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

                res.amp[i] = fres.params[0];
            }

            // Compute the error with MC simulation
            res.amp_bfit_sim.resize(params.nsim);
            res.amp_sim.resize(params.nsim);
            res.sed_sim.resize(params.nsim);

            const uint_t nflux = flux.size();
            for (uint_t i = 0; i < params.nsim; ++i) {
                auto fsim = flux;
                fsim[idm] += randomn(seed, idm.size());

                vec_t<1,ttype> chi2(nsed);
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

        vec_t<1,ttype> tmp1(nsed), tmp2(nsed);
        for (uint_t i = 0; i < nsed; ++i) {
            auto tmp = res.flux(i,_);
            tmp1[i] = total(weight*flux*tmp);
            tmp2[i] = total(weight*tmp*tmp);
        }

        res.amp = tmp1/tmp2;
        if (params.renorm) {
            tmp1 *= res.amp;
            res.chi2 = total(weight*flux*flux) - tmp1;
        } else {
            res.chi2 = total(weight*flux*flux) - 2.0*tmp1 + tmp2;
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
            vec_t<1,ttype> chi2;
            if (params.renorm) {
                tmp1 *= amp;
                chi2 = total(weight*fsim*fsim) - tmp1;
            } else {
                chi2 = total(weight*fsim*fsim) - 2.0*tmp1 + tmp2;
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

#endif
