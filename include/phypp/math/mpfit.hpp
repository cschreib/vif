#ifndef PHYPP_MATH_MPFIT_HPP
#define PHYPP_MATH_MPFIT_HPP

#include "phypp/math/math.hpp"

namespace phypp {
    // Note:
    // The following code is a direct translation of MPFIT, a routine written in IDL language by
    // C. Markwardt, itself inspired from the MINPACK Fortran library.
    // See: http://cow.physics.wisc.edu/~craigm/idl/idl.html
    //
    // The main differences between the IDL and C++ versions are:
    //  - In IDL, the MPFIT routine does support not recursive calls (i.e. fitting a model that itself
    //    calls MPFIT). There is a way around this problem, although it is tedious. The C++ version
    //    supports recursion naturally.
    //  - For simplicity, the C++ version does not support:
    //      - specifying explicit derivatives,
    //      - registering a callback for each iteration,
    //      - tied parameters.
    //
    // Slight differences may appear in the result of the two codes due to different floating point
    // policies between IDL and C++.

    struct mpfit_result {
        bool success = true;

        enum reason_t {
            unknown = 0,
            // failure
            overflow, max_iter,
            // success or failure
            gtol, ftol, xtol, ftol_xtol
        } reason = unknown;

        double chi2;
        uint_t dof;
        uint_t iter;
        vec1d params;
        vec1d errors;
        vec2d covar;
    };

    struct mpfit_options {
        mpfit_options() : nparam(0) {}

        explicit mpfit_options(uint_t np) : nparam(np) {
            upper_limit.resize(nparam);
            upper_limit[_] = dnan;
            lower_limit.resize(nparam);
            lower_limit[_] = dnan;
            frozen.resize(nparam);
            max_step.resize(nparam);
            deriv.resize(nparam);
            deriv_step.resize(nparam);
            deriv_rstep.resize(nparam);
        }

        enum deriv_t {
            deriv_auto = 0,      // forward if possible, else backward
            deriv_forward = 1,   // f(x+dx) - f(x) / dx
            deriv_backward = -1, // f(x) - f(x-dx) / dx
            deriv_symmetric = 2  // f(x+dx) - f(x-dx) / 2*dx
        };

        uint_t nparam;

        // Per parameter options
        vec1d upper_limit; // upper limit constraint on fit parameter (nan = none, default: nan)
        vec1d lower_limit; // lower limit constraint on fit parameter (nan = none, default: nan)
        vec1b frozen;      // keep fit parameter fixed (default: false)
        vec1d max_step;    // maximum step to allow in a single iteration (0 = none, default: 0)
        vec1i deriv;       // how derivates should be computed (see deriv_t above, default: auto)
        vec1d deriv_step;  // what absolute step to use to compute derivatives (0 = auto, default: 0)
        vec1d deriv_rstep; // what relative step to use to compute derivatives (0 = auto, default: 0)

        // Global options
        uint_t max_iter = 200; // maximum number of iteration to reach
        double gtol = 1e-10;   // maximum absolute change of chi2 to select solution
        double ftol = 1e-10;   // maximum relative change of chi2 to select solution
        double xtol = 1e-10;   // target relative error on the solution
        bool nocovar = false;  // do not compute errors and covariance matrix (faster)
    };

    // Numerically stable sqrt(total(sqr(v)))
    double mpfit_enorm(const vec1d& v) {
        const double dwarf = sqrt(std::numeric_limits<double>::min()*1.5)*10.0;
        const double giant = sqrt(std::numeric_limits<double>::max())*0.1;

        auto mm = minmax(v);
        double mx = max(abs(mm.first), abs(mm.second));
        if (mx == 0.0) return 0.0;

        if (mx > giant/v.size() || mx < dwarf*v.size()) {
            double res = 0.0;
            for (double x : v) {
                res += sqr(x/mx);
            }

            return mx*sqrt(res);
        } else {
            double res = 0.0;
            for (double x : v) {
                res += sqr(x);
            }

            return sqrt(res);
        }
    }

    // Compute the Jacobian matrix with finite difference derivatives
    template<typename F>
    vec2d mpfit_fdjac2(F&& deviate, const vec1d& xall, const vec1u& ifree, const vec1d& x,
        const vec1d& fvec, const mpfit_options& options) {

        const double eps = sqrt(std::numeric_limits<double>::epsilon());

        const uint_t m = fvec.size();
        const uint_t n = x.size();

        vec2d fjac(n, m);

        // Calculate the step
        vec1d h(n);
        for (uint_t p : range(n)) {
            uint_t ip = ifree[p];

            if (options.deriv_rstep[ip] != 0.0) {
                // If relative step is given, use that
                h[p] = abs(options.deriv_rstep[ip]*x[p]);
            } else if (options.deriv_step[ip] != 0.0) {
                // If step is given, use that
                h[p] = options.deriv_step[ip];
            } else {
                h[p] = eps*abs(x[p]);
            }

            // Prevent zero step
            if (h[p] == 0.0) h[p] = eps;

            // Reverse the sign of the step if using backward derivative or if using forward derivative
            // and the function would be evaluated beyond the provided upper limit.
            if (options.deriv[ip] == mpfit_options::deriv_backward ||
                (!is_nan(options.upper_limit[ip]) && x[p] > options.upper_limit[ip] - h[p])) {
                h[p] = -h[p];
            }
        }

        // Compute the matrix for each parameter
        for (uint_t p : range(n)) {
            uint_t ip = ifree[p];

            vec1d xp = xall;
            xp[ip] += h[p];
            vec1d fp = flatten(deviate(xp));

            if (options.deriv[ip] == mpfit_options::deriv_backward ||
                options.deriv[ip] == mpfit_options::deriv_forward ||
                options.deriv[ip] == mpfit_options::deriv_auto) {
                // One sided derivative
                fjac(p,_) = (fp - fvec)/h[p];
            } else {
                // Two sided derivative
                xp[ip] = xall[ip] - h[p];
                vec1d fm = flatten(deviate(xp));
                fjac(p,_) = (fp - fm)/(2.0*h[p]);
            }
        }

        return fjac;
    }

    // Compute QR factorization of a matrix such that a(ipiv,_) = q*r
    void mpfit_qrfac(vec2d& a, vec1u& ipiv, vec1d& rdiag, vec1d& acnorm) {
        const double eps = std::numeric_limits<double>::epsilon();

        const uint_t n = a.dims[0];
        const uint_t m = a.dims[1];

        acnorm.resize(n);
        for (uint_t i : range(n)) {
            acnorm[i] = mpfit_enorm(a(i,_));
        }

        rdiag = acnorm;
        vec1d wa = rdiag;
        ipiv = uindgen(n);

        uint_t minmn = std::min(n, m);
        for (uint_t i : range(minmn)) {
            {
                vec1u r = rgen(i,n-1);
                uint_t kmax;
                double rmax = max(rdiag[r], kmax);
                if (!is_nan(rmax) && kmax != 0) {
                    // Exchange rows
                    kmax += i;
                    std::swap(ipiv[i], ipiv[kmax]);
                    rdiag[kmax] = rdiag[i];
                    wa[kmax] = wa[i];
                }
            }

            vec1u r = rgen(i,m-1);
            uint_t li = ipiv[i];
            vec1d aii = a(li,r);
            double ajnorm = mpfit_enorm(aii);

            if (ajnorm != 0.0) {
                if (a(li,i) < 0.0) ajnorm = -ajnorm;

                aii /= ajnorm;
                aii[0] += 1;
                a(li,r) = aii;

                if (i+1 < n) {
                    for (uint_t k : range(i+1,n)) {
                        uint_t lk = ipiv[k];
                        vec1d aik = a(lk,r);
                        if (a(li,i) != 0.0) {
                            a(lk,r) = aik - aii*total(aik*aii)/a(li,i);
                        }

                        if (rdiag[k] != 0.0) {
                            double temp = a(lk,i)/rdiag[k];
                            temp = 1.0 - sqr(temp);
                            if (temp < 0.0) temp = 0.0;
                            rdiag[k] *= sqrt(temp);
                            temp = rdiag[k]/wa[k];

                            if (0.05*sqr(temp) <= eps) {
                                rdiag[k] = mpfit_enorm(a(lk,rgen(i+1,n-1)));
                                wa[k] = rdiag[k];
                            }
                        }
                    }
                }
            }

            rdiag[i] = -ajnorm;
        }
    }

    // Solve the linear system (Q*R)*x = B
    void mpfit_qrsolv(vec2d& r, const vec1u& ipiv, const vec1d& diag, const vec1d& qtf,
        vec1d& x, vec1d& sdiag) {

        const uint_t n = r.dims[0];
        const uint_t m = r.dims[1];

        vec1u delm = uindgen(n)*(m+1); // Diagonal elements of r

        // Copy r and (q transpose)*b to preserve input and initialize s.
        // In particular, save the diagonal elements of r in x.
        for (uint_t j : range(n)) {
            vec1u rj = rgen(j,n-1);
            r(j,rj) = r(rj,j);
        }

        x = r[delm];
        vec1d wa = qtf;

        // Eliminate the diagonal matrix d using a givens rotation
        for (uint_t j : range(n)) {
            uint_t l = ipiv[j];
            if (diag[l] == 0.0) break;
            sdiag[rgen(j,n-1)] = 0.0;
            sdiag[j] = diag[l];

            // The transformations to eliminate the row of d modify only a
            // single element of (q transpose)*b beyond the first n, which
            // is initially zero.
            double qtbpj = 0.0;

            for (uint_t k : range(j,n)) {
                if (sdiag[k] == 0.0) continue;

                double sine = 0.0, cosine = 0.0;
                if (abs(r(k,k)) < abs(sdiag[k])) {
                    double cotan = r(k,k)/sdiag[k];
                    sine = 0.5/sqrt(0.25 + 0.25*sqr(cotan));
                    cosine = sine*cotan;
                } else {
                    double tang = sdiag[k]/r(k,k);
                    cosine = 0.5/sqrt(0.25 + 0.25*sqr(tang));
                    sine = cosine*tang;
                }

                // Compute the modified diagonal element of r and the
                // modified element of ((q transpose)*b,0).
                r(k,k) = cosine*r(k,k) + sine*sdiag[k];
                double temp = cosine*wa[k] + sine*qtbpj;
                qtbpj = -sine*wa[k] + cosine*qtbpj;
                wa[k] = temp;

                // Accumulate the transformation in the row of s
                if (k < n-1) {
                    vec1u rk = rgen(k+1,n-1);
                    vec1d tmp = cosine*r(k,rk) + sine*sdiag[rk];
                    sdiag[rk] = -sine*r(k,rk) + cosine*sdiag[rk];
                    r(k,rk) = tmp;
                }
            }

            sdiag[j] = r(j,j);
            r(j,j) = x(j);
        }

        // Solve the triangular system for z.  If the system is singular
        // then obtain a least squares solution
        uint_t nsing = n;
        for (uint_t j : range(n)) {
            if (sdiag[j] == 0.0) {
                nsing = j;
                wa[rgen(j,n-1)] = 0.0;
                break;
            }
        }

        if (nsing >= 1) {
            wa[nsing-1] /= sdiag[nsing-1]; // Degenerate case
            for (uint_t j = nsing-1; j >= 1; --j) {
                vec1u rj = rgen(j,nsing-1);
                wa[j-1] -= total(r(j,rj)*wa(rj));
                wa[j-1] /= sdiag[j];
            }
        }

        // Permute the components of z back to components of x
        x[ipiv] = wa;
    }

    // Determine the levenberg-marquardt parameter
    void mpfit_lmpar(vec2d& r, const vec1u& ipiv, const vec1d& diag, const vec1d& qtf,
        double delta, vec1d& x, vec1d& sdiag, double& par) {

        const double dwarf = sqrt(std::numeric_limits<double>::min()*1.5)*10.0;
        const double eps = std::numeric_limits<double>::epsilon();

        const uint_t n = r.dims[0];
        const uint_t m = r.dims[1];

        vec1u delm = uindgen(n)*(m+1); // Diagonal elements of r

        // Compute and store in x the gauss-newton direction.  If the
        // jacobian is rank-deficient, obtain a least-squares solution
        uint_t nsing = n;
        vec1d wa1 = qtf;
        double rthresh = max(abs(r[delm]))*eps;
        for (uint_t j : range(delm)) {
            if (abs(r[delm[j]]) < rthresh) {
                nsing = j;
                wa1[rgen(j,n-1)] = 0.0;
                break;
            }
        }

        if (nsing > 0u) {
            for (uint_t j = nsing; j >= 1; --j) {
                wa1[j-1] /= r(j-1,j-1);
                if (j >= 2u) {
                    vec1u rj = rgen(0,j-2);
                    wa1[rj] -= r(j-1,rj)*wa1[j-1];
                }
            }
        }

        x[ipiv] = wa1;

        // Evaluate the function at the origin, and test for acceptance of
        // the gauss-newton direction
        vec1d wa2 = diag*x;
        double dxnorm = mpfit_enorm(wa2);
        double fp = dxnorm - delta;
        if (fp <= 0.1*delta) {
            par = 0.0;
            return;
        }

        // If the jacobian is not rank deficient, the newton step provides a
        // lower bound, parl, for the zero of the function.  Otherwise set
        // this bound to zero.
        double parl = 0.0;
        if (nsing >= n) {
            wa1 = diag[ipiv]*wa2[ipiv]/dxnorm;

            wa1[0] /= r(0,0); // Degenerate case
            for (uint_t j : range(1u,n)) {
                vec1u rj = rgen(0,j-1);
                wa1[j] -= total(r(j,rj)*wa1(rj));
                wa1[j] /= r(j,j);
            }

            double temp = mpfit_enorm(wa1);
            parl = ((fp/delta)/temp)/temp;
        }

        // Calculate an upper bound, paru, for the zero of the function
        for (uint_t j : range(n)) {
            vec1u rj = rgen(0,j);
            wa1[j] = total(r(j,rj)*qtf[rj])/diag[ipiv[j]];
        }

        double gnorm = mpfit_enorm(wa1);
        double paru = gnorm/delta;
        if (paru == 0.0) {
            paru = dwarf/std::min(delta, 0.1);
        }

        // If the input par lies outside of the interval (parl,paru), set
        // par to the closer endpoint
        par = clamp(par, paru, parl);
        if (par == 0.0) {
            par = gnorm/dxnorm;
        }

        for (uint_t iter = 0u; iter < 10u; ++iter) {
            // Evaluate the function at the current value of par
            if (par == 0.0) {
                par = std::max(dwarf, paru*0.001);
            }

            double temp = sqrt(par);
            wa1 = temp*diag;
            mpfit_qrsolv(r, ipiv, wa1, qtf, x, sdiag);
            wa2 = diag*x;
            dxnorm = mpfit_enorm(wa2);
            temp = fp;
            fp = dxnorm - delta;

            if (abs(fp) <= 0.1*delta || (parl == 0.0 && fp <= temp && temp < 0.0)) {
                break;
            }

            // Compute the newton correction
            wa1 = (diag*wa2)[ipiv]/dxnorm;

            for (uint_t j : range(n-1)) {
                wa1[j] /= sdiag[j];
                vec1u rj = rgen(j+1,n-1);
                wa1[rj] -= r(j,rj)*wa1[j];
            }
            wa1[n-1] /= sdiag[n-1]; // Degenerate case

            temp = mpfit_enorm(wa1);
            double parc = ((fp/delta)/temp)/temp;

            // Depending on the sign of the function, update parl or paru
            if (fp > 0) {
                parl = std::max(parl, par);
            } else if (fp < 0) {
                paru = std::min(paru, par);
            }

            // Compute an improved estimate for par
            par = std::max(parl, par+parc);
        }
    }

    // Compute convariance matrix from Jacobian
    vec2d mpfit_covar(vec2d r, const vec1u& ipiv) {
        phypp_check(r.dims[0] == r.dims[1], "matrix must be squared (got ", r.dims, ")");
        const uint_t n = r.dims[0];

        // Form the inverse of r in the full upper triangle of r
        uint_t l = 0;
        const double tol = 1e-14;
        double tolr = tol*abs(r(0,0));
        for (uint_t k : range(n)) {
            if (abs(r(k,k)) <= tolr) break;
            r(k,k) = 1.0/r(k,k);
            for (uint_t j : range(k)) {
                vec1u rj = rgen(0,j);
                double temp = r(k,k)*r(k,j);
                r(k,j) = 0.0;
                r(k,rj) -= temp*r(j,rj);
            }

            l = k+1;
        }

        // Form the full upper triangle of the inverse of (r transpose)*r
        // in the full upper triangle of r
        for (uint_t k : range(l)) {
            for (uint_t j : range(k)) {
                vec1u rj = rgen(0,j);
                r(j,rj) += r(k,j)*r(k,rj);
            }

            vec1u rk = rgen(0,k);
            r(k,rk) *= r(k,k);
        }

        // Form the full lower triangle of the covariance matrix
        // in the strict lower triangle of r and in wa
        vec1d wa(n);
        for (uint_t j : range(n)) {
            uint_t jj = ipiv[j];
            uint_t sing = j+1 > l;
            for (uint_t i : range(j)) {
                if (sing) {
                    r(j,i) = 0.0;
                }

                uint_t ii = ipiv[i];
                if (ii > jj) {
                    r(jj,ii) = r(j,i);
                } else if (ii < jj) {
                    r(ii,jj) = r(j,i);
                }
            }

            wa[jj] = r(j,j);
        }

        // Symmetrize the covariance matrix in r
        for (uint_t j : range(n)) {
            vec1u rj = rgen(0,j);
            r(j,rj) = r(rj,j);
            r(j,j) = wa[j];
        }

        return r;
    }

    template<typename F>
    mpfit_result mpfit(F&& deviate, vec1d xall, mpfit_options options = mpfit_options()) {
        const double eps = std::numeric_limits<double>::epsilon();

        mpfit_result res;

        if (options.nparam == 0u) {
            // No option provided, use defaults
            options = mpfit_options(xall.size());
        } else {
            phypp_check(options.nparam == xall.size(), "incompatible number of elements in options "
                "with provided parameters ("+strn(options.nparam)+" vs "+strn(xall.size())+")");
        }

        // Note: be carefull that limits are NaN by default, thus !(a == b) != (a != b)
        vec1u ifree = where(!options.frozen && !(options.upper_limit == options.lower_limit));

        vec1d llim = options.lower_limit[ifree];
        vec1d ulim = options.upper_limit[ifree];
        vec1u ilow = where(!is_nan(llim));
        vec1u iup = where(!is_nan(ulim));
        bool anylim = !ilow.empty() || !iup.empty();

        vec1d mastep = options.max_step[ifree];
        vec1u imima = where(mastep != 0.0);
        bool anymima = !imima.empty();

        const uint_t n = ifree.size();

        // Initialization
        vec1d x = xall[ifree];
        vec1d fvec = flatten(deviate(xall));
        double fnorm = mpfit_enorm(fvec);
        double fnorm1 = fnorm;
        vec1d qtf(n);

        const uint_t npt = fvec.size();

        res.dof = npt - n;
        uint_t iter = 1u;
        double factor = 100.0;
        double delta = dnan;
        double par = 0.0;
        double xnorm = dnan;
        vec1d diag;
        vec2d fjac;
        vec1u ipiv;

        // l.3243 mpfit.pro
        while (true) {
            fjac = mpfit_fdjac2(deviate, xall, ifree, x, fvec, options);

            // Set derivatives of frozen parameters to zero
            for (uint_t p : range(n)) {
                if (!is_nan(llim[p]) && total(fvec*fjac(p,_)) > 0.0) {
                    fjac(p,_) = 0.0;
                }
                if (!is_nan(ulim[p]) && total(fvec*fjac(p,_)) < 0.0) {
                    fjac(p,_) = 0.0;
                }
            }

            // Compute QR factorization of the Jacobian
            // l.3338 mpfit.pro
            vec1d wa1, wa2;
            mpfit_qrfac(fjac, ipiv, wa1, wa2);

            if (iter == 1u) {
                diag = wa2;
                for (uint_t p : range(n)) {
                    if (diag[p] == 0.0) diag[p] = 1.0;
                }

                vec1d wa3 = diag*x;
                xnorm = mpfit_enorm(wa3);
                delta = factor*xnorm;
                if (delta == 0.0) delta = factor;
            }

            // Form (Q transpose)*fvec and store the first n components in qtf
            // l.3371 mpfit.pro
            vec1d wa4 = fvec;
            for (uint_t p : range(n)) {
                vec1u r = rgen(p, npt-1);
                uint_t lp = ipiv[p];
                if (fjac(lp,p) != 0.0) {
                    wa4[r] -= fjac(lp,r)*total(fjac(lp,r)*wa4[r])/fjac(lp,p);
                }

                fjac(lp,p) = wa1(p);
                qtf[p] = wa4[p];
            }

            // Reform the Jacobian matrix (only need square R factor)
            // l.3388 mpfit.pro
            fjac = fjac(ipiv,rgen(0,n-1));

            // Check for overflow
            bool stop = false;
            for (double d : fjac) {
                if (!is_finite(d)) {
                    res.success = false;
                    res.reason = mpfit_result::overflow;
                    stop = true;
                }
            }

            if (stop) break;

            // Compute the norm of the scaled gradient
            // l.3401 mpfit.pro
            double gnorm = 0.0;
            if (fnorm != 0.0) {
                for (uint_t p : range(n)) {
                    vec1u r = rgen(0,p);
                    uint_t l = ipiv[p];
                    if (wa2[l] != 0.0) {
                        gnorm = std::max(gnorm, abs(total(fjac(p,r)*qtf[r])/wa2[l]));
                    }
                }
            }

            // Possible success condition
            if (gnorm < options.gtol) {
                res.success = true;
                res.reason = mpfit_result::gtol;
                break;
            }

            // Rescale if necessary
            diag = max(diag, wa2);

            double ratio = 0.0;
            bool success = false;
            while (!success) {
                // Determine the levenberg-marquardt parameter
                // l.3429 mpfit.pro
                mpfit_lmpar(fjac, ipiv, diag, qtf, delta, wa1, wa2, par);

                // Store the direction p and x+p. Calculate the norm of p
                wa1 *= -1.0;

                double alpha = 1.0;
                if (!anylim && !anymima) {
                    // No parameter limits, so just move to new position WA2
                    wa2 = x + wa1;
                } else {
                    // Respect the limits.  If a step were to go out of bounds, then
                    // we should take a step in the same direction but shorter distance.
                    // The step should take us right to the limit in that case.
                    if (anylim) {
                        wa1[ilow] = max(wa1[ilow], 0.0);
                        wa1[iup] = min(wa1[iup], 0.0);

                        vec1b dwa1 = abs(wa1) > eps;
                        vec1u whl = where(dwa1 && x + wa1 < llim);
                        if (!whl.empty()) {
                            double mi = min((llim[whl] - x[whl])/wa1[whl]);
                            alpha = std::min(alpha, mi);
                        }

                        vec1u whu = where(dwa1 && x + wa1 > ulim);
                        if (!whu.empty()) {
                            double mi = min((ulim[whu] - x[whu])/wa1[whu]);
                            alpha = std::min(alpha, mi);
                        }
                    }

                    if (anymima) {
                        vec1d nwa1 = wa1*alpha;
                        double mrat = max(abs(nwa1[imima])/abs(mastep[imima]));
                        if (mrat > 1.0) {
                            alpha /= mrat;
                        }
                    }

                    // Scale the resulting vector
                    wa1 *= alpha;
                    wa2 = x + wa1;

                    if (anylim) {
                        // Adjust the final output values.  If the step put us exactly
                        // on a boundary, make sure we peg it there.
                        //                ... nonzero *LIM ...       ... zero *LIM ...
                        vec1d ulim1 = ulim*(1.0 - sign(ulim)*eps) - (ulim == 0.0)*eps;
                        vec1d llim1 = llim*(1.0 + sign(llim)*eps) + (llim == 0.0)*eps;

                        vec1u wh = where(wa2 >= ulim1);
                        wa2[wh] = ulim[wh];
                        wh = where(wa2 <= llim1);
                        wa2[wh] = llim[wh];
                    }
                }

                vec1d wa3 = diag*wa1;
                double pnorm = mpfit_enorm(wa3);

                // On the first iteration, adjust the initial step bound
                if (iter == 1u) {
                    delta = std::min(delta, pnorm);
                }

                xall[ifree] = wa2;

                // Evaluate the function at x+p and calculate its norm
                wa4 = flatten(deviate(xall));
                fnorm1 = mpfit_enorm(wa4);

                // Compute the scaled actual reduction
                double actred = -1.0;
                if (0.1*fnorm1 < fnorm) {
                    actred = 1.0 - sqr(fnorm1/fnorm);
                }

                // Compute the scaled predicted reduction and the scaled directional
                // derivative
                for (uint_t j : range(n)) {
                    wa3[j] = 0.0;
                    vec1u rj = rgen(0,j);
                    wa3[rj] += fjac(j,rj)*wa1[ipiv[j]];
                }

                // Remember, alpha is the fraction of the full LM step actually
                // taken
                double prered, dirder; {
                    double temp1 = mpfit_enorm(alpha*wa3)/fnorm;
                    double temp2 = (sqrt(alpha*par)*pnorm)/fnorm;
                    prered = sqr(temp1) + sqr(temp2)/0.5;
                    dirder = -(sqr(temp1) + sqr(temp2));
                }

                // Compute the ratio of the actual to the predicted reduction.
                if (prered != 0.0) {
                    ratio = actred/prered;
                }

                // Update the step bound
                if (ratio <= 0.25) {
                    double temp = 0.5;
                    if (actred < 0.0) {
                        temp = 0.5*dirder/(dirder + 0.5*actred);
                    }

                    if (0.1*fnorm1 >= fnorm || temp < 0.1) {
                        temp = 0.1;
                    }

                    delta = temp*std::min(delta, pnorm/0.1);
                    par /= temp;
                } else if (par == 0.0 || ratio >= 0.75) {
                    delta = pnorm/0.5;
                    par = 0.5*par;
                }

                if (ratio >= 0.0001) {
                    // Successful iteration.  Update x, fvec, and their norms
                    x = wa2;
                    wa2 = diag*x;

                    fvec = wa4;
                    xnorm = mpfit_enorm(wa2);
                    fnorm = fnorm1;
                    ++iter;

                    success = true;
                }

                // Tests for convergence
                bool ftol = abs(actred) <= options.ftol && prered <= options.ftol && 0.5*ratio <= 1.0;
                bool xtol = delta <= options.xtol*xnorm;
                if (ftol || xtol) {
                    res.success = true;
                    if (ftol && xtol) {
                        res.reason = mpfit_result::ftol_xtol;
                    } else if (ftol) {
                        res.reason = mpfit_result::ftol;
                    } else {
                        res.reason = mpfit_result::xtol;
                    }
                    stop = true;
                    break;
                }

                // Tests for termination and stringent tolerances
                if (iter >= options.max_iter) {
                    res.success = false;
                    res.reason = mpfit_result::max_iter;
                    stop = true;
                    break;
                }

                if (abs(actred) <= eps && prered <= eps && 0.5*ratio <= 1.0) {
                    res.success = false;
                    res.reason = mpfit_result::ftol;
                    stop = true;
                    break;
                }

                if (delta <= eps*xnorm) {
                    res.success = false;
                    res.reason = mpfit_result::xtol;
                    stop = true;
                    break;
                }

                if (gnorm <= eps) {
                    res.success = false;
                    res.reason = mpfit_result::gtol;
                    stop = true;
                    break;
                }
            }

            if (stop) break;

            // Check for over/underflow
            if (!is_finite(ratio) || count(!is_finite(wa1)) != 0 ||
                count(!is_finite(wa2)) != 0 || count(!is_finite(x)) != 0) {
                res.success = false;
                res.reason = mpfit_result::overflow;
                break;
            }
        }

        xall[ifree] = x;

        res.params = xall;
        fvec = flatten(deviate(xall));
        fnorm = mpfit_enorm(fvec);
        res.chi2 = sqr(std::max(fnorm, fnorm1));
        res.iter = iter;

        if (!options.nocovar) {
            // (very carefully) set the covariance matrix
            vec2d cv = mpfit_covar(fjac, ipiv);

            // Fill in actual covariance matrix, accounting for fixed
            // parameters.
            res.covar.resize(xall.size(), xall.size());
            res.covar(ifree,ifree) = cv;

            // Compute errors in parameters
            res.errors = matrix::diagonal(res.covar);
            vec1u wh = where(res.errors >= 0.0);
            res.errors[wh] = sqrt(res.errors[wh]);
        }

        return res;
    }

    // Wrapper around mpfit() for standard deviate (y - ytest)/yerr, where y and yerr are given and
    // ytest is compted from a model function taking as a first argument the position x at which to
    // compute the model (given) and the function parameters p (to be found).
    template<typename F, typename TX, typename TY, typename TYE>
    mpfit_result mpfitfun(const TY& y, const TYE& ye, const TX& x, F&& model,
        const vec1d& params, const mpfit_options& options = mpfit_options()) {

        bool bad = !same_dims_or_scalar(x, y, ye);
        if (bad) {
            phypp_check(same_dims_or_scalar(x, y), "incompatible dimensions between X and Y arrays "
                "("+strn(dim(x))+" vs. "+strn(dim(y))+")");
            phypp_check(same_dims_or_scalar(x, ye), "incompatible dimensions between X and YE arrays "
                "("+strn(dim(x))+" vs. "+strn(dim(ye))+")");
        }

        return mpfit([&](const vec1d& p) {
            return (y - model(x,p))/ye;
        }, params, options);
    }
}

#endif
