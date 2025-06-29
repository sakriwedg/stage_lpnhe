//
// Created by Mathieu GUIGUE on 10/11/2022.
//

#include "THATFitMethods.hxx"

#include <iostream>
#include <cmath>

int ND::THATFitMethods::matinv3(double *a, double *b)
{
   double r11, r12, r22, r13, r23, r33, s11, s12, s22, s13, s23, s33;

   if (a[0] <= 0)
      return 1;
   r11 = sqrt(a[0]);
   r12 = a[1] / r11;
   r13 = a[3] / r11;
   if ((r22 = a[2] - r12 * r12) <= 0)
      return 2;
   r22 = sqrt(r22);
   r23 = (a[4] - r12 * r13) / r22;
   if ((r33 = a[5] - r13 * r13 - r23 * r23) <= 0)
      return 3;
   r33 = sqrt(r33);

   s11 = 1 / r11;
   s12 = -s11 * r12 / r22;
   s13 = -(s11 * r13 + s12 * r23) / r33;
   s22 = 1 / r22;
   s23 = -(s22 * r23) / r33;
   s33 = 1 / r33;

   b[0] = s11 * s11 + s12 * s12 + s13 * s13;
   b[1] = s12 * s22 + s13 * s23;
   b[2] = s22 * s22 + s23 * s23;
   b[3] = s13 * s33;
   b[4] = s23 * s33;
   b[5] = s33 * s33;

   return 0;
}
/*------------------------------------------------------------------------------------------*/
int ND::THATFitMethods::FitCirc(int np, double *z, double *y, double *sig, double *pari, double *parf, double &chi2,
                                double *C)
{

   double T[3], W[6], dpar[3], der[3];
   int NITER = 2;

   for (int i = 0; i < 3; i++)
      parf[i] = pari[i];
   for (int iter = 0; iter < NITER; iter++) {
      for (int i = 0; i < 3; i++)
         T[i] = 0;
      for (int i = 0; i < 6; i++)
         W[i] = 0;
      chi2 = 0;
      for (int ip = 0; ip < np; ip++) {
         double dz, dy, dzy, dr;
         dz = z[ip] - parf[0];
         dy = y[ip] - parf[1];
         dzy = sqrt(dz * dz + dy * dy);
         dr = dzy - parf[2];
         double wgt = 1. / (sig[ip] * sig[ip]);
         chi2 += wgt * dr * dr;
         // summation of coefficients of the linear system ----------------
         der[0] = dz / dzy;
         der[1] = dy / dzy;
         der[2] = 1;
         for (int i = 0; i < 3; i++)
            T[i] += wgt * der[i] * dr;
         W[0] += wgt * der[0] * der[0];
         W[1] += wgt * der[0] * der[1];
         W[2] += wgt * der[1] * der[1];
         W[3] += wgt * der[0] * der[2];
         W[4] += wgt * der[1] * der[2];
         W[5] += wgt * der[2] * der[2];
      }
      int flag = matinv3(W, C);
      if (flag)
         return 0;
      dpar[0] = C[0] * T[0] + C[1] * T[1] + C[3] * T[2];
      dpar[1] = C[1] * T[0] + C[2] * T[1] + C[4] * T[2];
      dpar[2] = C[3] * T[0] + C[4] * T[1] + C[5] * T[2];
      for (int i = 0; i < 3; i++)
         parf[i] += dpar[i];
      double dchi2;
      dchi2 = W[0] * dpar[0] * dpar[0] + W[2] * dpar[1] * dpar[1] + W[5] * dpar[2] * dpar[2] +
              2 * (W[1] * dpar[0] * dpar[1] + W[3] * dpar[0] * dpar[2] + W[4] * dpar[1] * dpar[2]);
      chi2 -= dchi2;
   }
   return 1;
}
/*----------------------------------------------------------------------------------------------*/
int ND::THATFitMethods::FitParab(int nz, double *z, double *y, double *sigy, double *parf, double &chi2, double *C)
{
   double T[3], W[6];

   for (int i = 0; i < 3; i++)
      T[i] = 0;
   for (int i = 0; i < 6; i++)
      W[i] = 0;
   chi2 = 0;
   for (int iz = 0; iz < nz; iz++) {
      double wgt = 1 / (sigy[iz] * sigy[iz]);
      double zi = z[iz], zi2 = zi * zi, yi = y[iz];
      chi2 += wgt * yi * yi;
      T[0] += wgt * yi;
      T[1] += wgt * zi * yi;
      T[2] += wgt * zi2 * yi;
      W[0] += wgt;
      W[1] += wgt * zi;
      W[2] += wgt * zi2;
      W[3] += wgt * zi2;
      W[4] += wgt * zi * zi2;
      W[5] += wgt * zi2 * zi2;
   }
   int flag = matinv3(W, C);
   if (flag) {
      printf("!!!invert %d\n", flag);
      return 0;
   }
   parf[0] = C[0] * T[0] + C[1] * T[1] + C[3] * T[2];
   parf[1] = C[1] * T[0] + C[2] * T[1] + C[4] * T[2];
   parf[2] = C[3] * T[0] + C[4] * T[1] + C[5] * T[2];
   double dchi2;
   dchi2 = W[0] * parf[0] * parf[0] + W[2] * parf[1] * parf[1] + W[5] * parf[2] * parf[2] +
           +2 * (W[1] * parf[0] * parf[1] + W[3] * parf[0] * parf[2] + W[4] * parf[1] * parf[2]);
   chi2 -= dchi2;
   return 1;
}
/*------------------------------------------------------*/
void ND::THATFitMethods::Shape(int np, double *z, double *y, double &phi, double *zp, double *up, double &length,
                               double &width, int &ipmin, int &ipmax, int &ipmed)
{
   // defining an average direction and global shape
   double zpad, ypad, sum, sumz, sumy, sumzz, sumyy, sumzy, sdelta, zavg, yavg;
   double cphi, sphi;
   sum = sumz = sumy = sumzz = sumyy = sumzy = 0;
   for (int ip = 0; ip < np; ip++) {
      zpad = z[ip];
      ypad = y[ip];
      sum += 1;
      sumz += zpad;
      sumy += ypad;
      sumzz += zpad * zpad;
      sumzy += zpad * ypad;
      sumyy += ypad * ypad;
   }
   zavg = sumz / sum;
   yavg = sumy / sum;
   sumzz /= sum;
   sumyy /= sum;
   sumzy /= sum;
   sumzz -= zavg * zavg;
   sumyy -= yavg * yavg;
   sumzy -= zavg * yavg;
   phi = atan2(2 * sumzy, sumzz - sumyy) / 2;
   cphi = cos(phi);
   sphi = sin(phi);
   sdelta = sqrt(pow(sumzz - sumyy, 2) + 4 * pow(sumzy, 2));
   length = sqrt((sumzz + sumyy + sdelta) / 2);
   width = sqrt((sumzz + sumyy - sdelta) / 2);

   // finding the first, last and medium point along this direction
   double zmin, zmax, zmed, dzmed;
   ipmin = ipmax = 0;
   zmin = zmax = zp[0] = cphi * z[0] + sphi * y[0];
   for (int ip = 0; ip < np; ip++) {
      zp[ip] = cphi * z[ip] + sphi * y[ip];
      if (zp[ip] < zmin) {
         ipmin = ip;
         zmin = zp[ip];
      }
      if (zp[ip] > zmax) {
         ipmax = ip;
         zmax = zp[ip];
      }
      up[ip] = -sphi * z[ip] + cphi * y[ip];
   }
   zmed = (zmin + zmax) / 2;
   ipmed = 0;
   dzmed = fabs(zp[0] - zmed);
   for (int ip = 0; ip < np; ip++) {
      double dz = fabs(zp[ip] - zmed);
      if (dz < dzmed) {
         ipmed = ip;
         dzmed = dz;
      }
   }
}

// input ------------------
// np: number of points
//     z,y: coordinates
//     dirm: local direction for clustering: 0 for horiz, 1 for diag-up, 2 for vert, 3 for diag-down
//     sig: error on measurement
// return value -----------------
// 0 if rejected (mainly because of short length)
// 1 if parabolic fit
// 2 if circle fit
// output -----------------
// zp,up: coordinates along the longitudinal/transverse directions of the principal axes
// param,chi2: result of the fit
// param:
//- circle fit: Z,Y (center) , R
//- parabolic fit: p0,p1,p2 in the equation up = p0+p1*(zp-zm)+p2*(zp-zm)^2 with zm = middle point
//     zfit,yfit,phifit,R: parameters of the "reference point" on the trajectory
//  cov[9] : 3x3 covariance matrix on Z Y and Rfit. Option lower down here to use transverse coordinates instead
int ND::THATFitMethods::FitHelix(int np, double *z, double *y, int *dirm, double *sig, double *zp, double *up,
                                 double &zfit, double &yfit, double &phifit, double &Rfit, double *param, double &chi2,
                                 double *cov, bool forceCircFit)

#define LENGTH_MIN 10.
#define WIDTH_OVER_LENGTH_MAX .05
{
   if (np == 0) {
      std::cout <<"No points in the track" << std::endl;
      return 0;
   }

   double chi2_circ, param_circ[3], cov_circ[6], chi2_parab, param_parab[3], cov_parab[6];
   double phi, cphi, sphi, length, width;
   double zm, ym, z1, y1, z2, y2, dz1, dy1, dz2, dy2, det, dzy1, dzy2, dz, dy;
   int ipmin, ipmax, ipmed;
   double param_ini[3];
   // double phidir[4] = {0., M_PI / 4, M_PI / 2, 3 * M_PI / 4};
   double phidir[4] = {M_PI / 2, 3 * M_PI / 4, 0, M_PI / 4};
   // define the principal axes (phi is the direction of the longitudinal axis)
   Shape(np, z, y, phi, zp, up, length, width, ipmin, ipmax, ipmed);

   if (length < LENGTH_MIN)
      return 0;
   cphi = cos(phi);
   sphi = sin(phi);

   // if(width/length<WIDTH_OVER_LENGTH_MAX and not forceCircFit) {

   // parabolic fit in a rotated frame -----------------
   // projecting the error onto the transverse axis
   for (int ip = 0; ip < np; ip++)
      sig[ip] *= fabs(cos(phidir[dirm[ip]] - phi));
   // shift the coordinates on the longitudinal axis
   double zmed;
   zmed = zp[ipmed];
   for (int ip = 0; ip < np; ip++)
      zp[ip] -= zmed;
   // -------------------
   auto parab_fit_res = FitParab(np, zp, up, sig, param_parab, chi2_parab, cov_parab);
   if (parab_fit_res == 0) {
       std::cout << "Parabolic fit failed"<< std::endl;
      return 0;
   }
   double z0 = -param_parab[1] / param_parab[2] / 2;
   double u0 = param_parab[0] + param_parab[1] * z0 + param_parab[2] * z0 * z0;
   zfit = cphi * (z0 + zmed) - sphi * u0;
   yfit = sphi * (z0 + zmed) + cphi * u0;
   phifit = phi + param_parab[1] + 2 * param_parab[2] * z0;
   Rfit = 1. / (2 * param_parab[2]);
   chi2 = chi2_parab; // Setting the chi2 to parabola chi2
                      //    if ( Rfit < 0 ) {
                      //        Rfit *= -1;
                      //        phifit += M_PI;
                      //    }
   for (int i = 0; i < 3; i++)
      param[i] = param_parab[i];
   for (int i = 0; i < 6; i++)
      cov[i] = cov_parab[i];

    //std::cout <<"parabola output: p1 = " << param[0] << ", p2 = " << param[1] << ", p3 = " << param[2]<< std::endl;

   // Here we compute the proper radius and center of the circle without changing the Rfit value extracted from
   // the parabolic fit: this was done above by
   //// if ( Rfit < 0 ) {
   ////        Rfit *= -1;
   ////        phifit += M_PI;
   ////    }
   double center_circle_Z, center_circle_Y;
   if (Rfit < 0) {
      Rfit *= -1.;
      center_circle_Z = zfit - Rfit * sin(phifit + M_PI);
      center_circle_Y = yfit + Rfit * cos(phifit + M_PI);
   } else {
      center_circle_Z = zfit - Rfit * sin(phifit);
      center_circle_Y = yfit + Rfit * cos(phifit);
   }


   // circular fit ------------------------------/
   // define a circle with the first,median and last points
   z1 = z[ipmin];
   y1 = y[ipmin];
   zm = z[ipmed];
   ym = y[ipmed];
   z2 = z[ipmax];
   y2 = y[ipmax];
   dz1 = z1 - zm;
   dy1 = y1 - ym;
   dz2 = z2 - zm;
   dy2 = y2 - ym;
   det = dz1 * dy2 - dz2 * dy1;
   if (fabs(det) < 0.5 || dz1 == 0 || dy1 == 0 || dz2 == 0 || dy2 == 0) {
      z1 = z[ipmin];
      y1 = y[ipmin];
      zm = z[ipmed] + 1.0;
      ym = y[ipmed] + 1.0;
      z2 = z[ipmax];
      y2 = y[ipmax];
      dz1 = z1 - zm;
      dy1 = y1 - ym;
      dz2 = z2 - zm;
      dy2 = y2 - ym;
      det = dz1 * dy2 - dz2 * dy1;
   }
   if (det == 0) {
       std::cout <<"Circular fit failed because of det=0"<< std::endl;
      return 0;
   }
   dzy1 = dz1 * dz1 + dy1 * dy1;
   dzy2 = dz2 * dz2 + dy2 * dy2;
   param_ini[0] = zm + (dy2 * dzy1 - dy1 * dzy2) / det / 2;
   param_ini[1] = ym + (-dz2 * dzy1 + dz1 * dzy2) / det / 2;
   dz = param_ini[0] - zm;
   dy = param_ini[1] - ym;
   param_ini[2] = sqrt(dz * dz + dy * dy);

   param_ini[0] = center_circle_Z;
   param_ini[1] = center_circle_Y;
   param_ini[2] = Rfit;

   // projecting the error onto the radial direction
   for (int ip = 0; ip < np; ip++) {
      double phip = atan2(param_ini[1] - y[ip], param_ini[0] - z[ip]);
      sig[ip] *= fabs(sin(phidir[dirm[ip]] - phip));
   }
   if (!FitCirc(np, z, y, sig, param_ini, param_circ, chi2_circ, cov_circ)) {
      if (parab_fit_res == 1) {
         return 1;
      }
      return 0;
   }
   double Z, Y, R;
   Z = param_circ[0];
   Y = param_circ[1];
   R = param_circ[2];
   chi2 = chi2_circ; // Setting the chi2 to parabola chi2
   // reverse phi if the center is on the right side of the main axis
   double phim, cphis, sphis;
   cphis = cphi;
   sphis = sphi;
   phim = atan2(-dz, dy);
   phifit = phi;
   if (fabs(fmod(phim - phi + 3 * M_PI, 2 * M_PI) - M_PI) > M_PI / 2) {
      phifit += M_PI;
      cphis = -cphis;
      sphis = -sphis;
   }
   zfit = Z + R * sphis;
   yfit = Y - R * cphis;
   Rfit = R;
   cov[0] = cov_circ[0];
   cov[1] = cov_circ[1];
   cov[2] = cov_circ[2];
   cov[3] = cov_circ[3];
   cov[4] = cov_circ[4];
   cov[5] = cov_circ[5];
   // curvfit = 1 / R;
   // // transformation of the covariance matrix from (Z,Y,R) to (u,phi,1/R)
   // double duZ, duY, duR, dfZ, dfY, dfR, dcZ, dcY, dcR;
   // duZ = -sphis;
   // duY = cphis;
   // duR = -1;
   // dfZ = -cphis / R;
   // dfY = -sphis / R;
   // dfR = 0;
   // dcZ = dcY = 0;
   // dcR = -1 / R / R;
   // cov[0] = duZ * duZ * cov_circ[0] + duY * duY * cov_circ[2] + duR * duR * cov_circ[5] +
   //          2 * (duZ * duY * cov_circ[1] + duZ * duR * cov_circ[3] + duY * duR * cov_circ[4]);
   // cov[1] = duZ * dfZ * cov_circ[0] + (duZ * dfY + duY * dfZ) * cov_circ[1] + duY * dfY * cov_circ[2] +
   //          (duZ * dfR + duR * dfZ) * cov_circ[3] + (duY * dfR + duR * dfY) * cov_circ[4] + duR * dfR * cov_circ[5];
   // cov[2] = dfZ * dfZ * cov_circ[0] + dfY * dfY * cov_circ[2] + dfR * dfR * cov_circ[5] +
   //          2 * (dfZ * dfY * cov_circ[1] + dfZ * dfR * cov_circ[3] + dfY * dfR * cov_circ[4]);
   // cov[3] = duZ * dcZ * cov_circ[0] + (duZ * dcY + duY * dcZ) * cov_circ[1] + duY * dcY * cov_circ[2] +
   //          (duZ * dcR + duR * dcZ) * cov_circ[3] + (duY * dcR + duR * dcY) * cov_circ[4] + duR * dcR * cov_circ[5];
   // cov[4] = dfZ * dcZ * cov_circ[0] + (dfZ * dcY + dfY * dcZ) * cov_circ[1] + dfY * dcY * cov_circ[2] +
   //          (dfZ * dcR + dfR * dcZ) * cov_circ[3] + (dfY * dcR + dfR * dcY) * cov_circ[4] + dfR * dcR * cov_circ[5];
   // cov[5] = dcZ * dcZ * cov_circ[0] + dcY * dcY * cov_circ[2] + dcR * dcR * cov_circ[5] +
   //          2 * (dcZ * dcY * cov_circ[1] + dcZ * dcR * cov_circ[3] + dcY * dcR * cov_circ[4]);

   // std::cout <<"Circular fit done successfully"<< std::endl;

   return 2;
}
