//
// Created by Mathieu GUIGUE on 10/11/2022.
//

#ifndef HATRECON_THATFITMETHODS_HXX
#define HATRECON_THATFITMETHODS_HXX

namespace ND {
namespace THATFitMethods {
/*! \brief Inversion of a symmetric positive matrix (Cholesky method)
 *
 *  @param a input matrix
 *  @param b output inverted matrix
 */
int matinv3(double *a, double *b);
/*! \brief Fit a circle on data assuming that there are errors on y
 *
 * Fits a circle using data (z,y) with errors on y.
 * Uses an iterative linear approach to reach the minimum of the chi2
 * Parameters are y(zc),phi(zc),sign/R
 * @param np Number of data points
 * @param z Pointer to the z data (size np)
 * @param y Pointer to the y data (size np)
 * @param sig Pointer to the errors on y (size np)
 * @param pari Initial guessed parameters
 * @param parf Final parameters extracted from the fit
 * @param chi2 Value of the Chi2 after the fit
 * @param C Covariance matrix projected as a vector of size 6 (only returning the upper triangular part of the
 * covariance matrix)
 */
int FitCirc(int np, double *z, double *y, double *sig, double *pari, double *parf, double &chi2, double *C);
/*! \brief Fit a parabola on data assuming that there are errors on y
 *
 * Fits a parabola using analytical formula for chi2 minimization and data (z,y) with errors on y.
 * Parameters are parf[2]+z*parf[1]*z*z*parf[0]
 * @param np Number of data points
 * @param z Pointer to the z data (size np)
 * @param y Pointer to the y data (size np)
 * @param sig Pointer to the errors on y (size np)
 * @param parf Final parameters extracted from the fit
 * @param chi2 Value of the Chi2 after the fit
 * @param C Covariance matrix projected as a vector of size 6 (only returning the upper triangular part of the
 * covariance matrix)
 */
int FitParab(int nz, double *z, double *y, double *sigy, double *parf, double &chi2, double *C);
void Shape(int np, double *z, double *y, double &phi, double *zp, double *up, double &length, double &width, int &ipmin,
           int &ipmax, int &ipmed);
int FitHelix(int np, double *z, double *y, int *dirm, double *sig, double *zp, double *up, double &zfit, double &yfit,
             double &phifit, double &Rfit, double *param, double &chi2, double *cov, bool forceCircFit = false);

} // namespace THATFitMethods
} // namespace ND

#endif // HATRECON_THATFITMETHODS_HXX
