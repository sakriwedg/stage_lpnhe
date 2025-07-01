//
// Created by Wiliam SAENZ AREVALO on 20/05/2025.
//
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMath.h>
#include <TFile.h>
#include <fstream>
#include <TPaveText.h>
#include <cmath>
#include <iostream>
#include <iomanip>  // for std::setw
#include <chrono>
#include <thread>   // for std::this_thread::sleep_for
#include <filesystem>

int bord = 40;
int N_hits = 1000000;  // nombre max de hits à traiter
int ADC = 10; //Amplitude en dessous de laquelle on a du bruit
int X_charge = 18000;
int X_amplitude = 2500;
int bin_amplitude = 200;
int bin_charge = 100;
#define LENGTH_MIN 10.
#define WIDTH_OVER_LENGTH_MAX .05


void plot_wf(int argc, char** argv);





bool save_plots = false; // prints example plots of track fit
bool use_integral = true; // if false it uses amplitude for barycenter calculation

void printProgressBar(float progress) {
    const int barWidth = 50;
    std::cout << "[";
    int pos = static_cast<int>(barWidth * progress);

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }

    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

struct TClusterFeatures {
    std::vector<double> time, phi, x, z, y, dir_x, dir_y, dir_z, err, q;
    std::vector<double> y_eram, z_eram;
    std::vector<int> ori, size;

    void clear() {
        time.clear(); phi.clear(); x.clear(); z.clear(); y.clear();
        dir_x.clear(); dir_y.clear(); dir_z.clear(); err.clear(); q.clear();
        y_eram.clear(); z_eram.clear();
        ori.clear(); size.clear();
    }
};

double GetBarycenter_Z(double z0, double z1, double z2,
                       double q0, double q1, double q2,
                       bool h0, bool h1, bool h2) {
    double sum = 0.0, total = 0.0;

    if (h0) { sum += z0 * q0; total += q0; }
    if (h1) { sum += z1 * q1; total += q1; }
    if (h2) { sum += z2 * q2; total += q2; }

    return (total == 0.0) ? 0.0 : sum / total;
}

double GetBarycenter_Z_Error(double z0, double z1, double z2,
                              double q0, double q1, double q2,
                              bool h0, bool h1, bool h2) {
    double sum = 0.0, weight = 0.0;

    if (h0) { sum += z0 * q0; weight += q0; }
    if (h1) { sum += z1 * q1; weight += q1; }
    if (h2) { sum += z2 * q2; weight += q2; }

    if (weight == 0.0) return 0.0;

    double zbar = sum / weight;
    double var = 0.0;

    if (h0) { double dz = z0 - zbar; var += q0 * dz * dz; }
    if (h1) { double dz = z1 - zbar; var += q1 * dz * dz; }
    if (h2) { double dz = z2 - zbar; var += q2 * dz * dz; }

    return std::sqrt(var / weight);
}

Double_t GetIntegral(Double_t thr, TGraph* wf) {
    double integral = 0.;

    for (int j = 0; j < wf->GetN(); ++j)
        if (wf->GetPointY(j) > thr) integral += wf->GetPointY(j);

    return integral;
}

void FillMinMaxValues(TClusterFeatures &clus_features, double cenZ, double cenY,
                      double &phi_min, double &phi_max) {
    for (unsigned int i = 0; i < clus_features.z.size(); ++i) {
        double dZ = clus_features.z[i] - cenZ;
        double dY = clus_features.y[i] - cenY;
        double phi = atan2(dY, dZ);

        clus_features.phi.push_back(phi);

        if (clus_features.size[i] > 0) {
            if (phi_min > phi) phi_min = phi;
            if (phi_max < phi) phi_max = phi;
        }
    }
}


void ProcessEvent(
    TClusterFeatures &clu_feat,
    TTree *tree,
    Long64_t event,
    TH1F *q0_amp_hist, TH1F *q1_amp_hist, TH1F *q2_amp_hist,
    TH1F *q0_int_hist, TH1F *q1_int_hist, TH1F *q2_int_hist,
    TH1F *q1_0_amp_hist, TH1F *q1_0_int_hist,
    TH1F *q2_1_amp_hist, TH1F *q2_1_int_hist
) {
    tree->GetEntry(event);

    // Structures temporaires pour regrouper les hits par y
    std::map<double, std::vector<int Hit>> hits_by_y;

    for (int i = 0; i < nhits; i++) {
        if (hit_chip_id[i] == 999) continue; // hit invalide ?

        // Sélection des hits valides
        if (use_integral && hit_integral[i] < min_integral) continue;
        if (!use_integral && hit_amplitude[i] < min_amplitude) continue;

        // Sélection temporelle
        if (hit_time[i] < t_min || hit_time[i] > t_max) continue;

        double y = hit_pos_y[i];
        hits_by_y[y].push_back({hit_pos_z[i], hit_amplitude[i], hit_integral[i], true});
    }

    for (const auto &[y0, hits] : hits_by_y) {
        if (hits.size() < 3) continue;

        // Tri des hits par position z croissante
        std::vector<Hit> sorted_hits = hits;
        std::sort(sorted_hits.begin(), sorted_hits.end(), [](const Hit &a, const Hit &b) {
            return a.z < b.z;
        });

        // Balayage des hits trois par trois pour former des "clusters"
        for (size_t i = 0; i + 2 < sorted_hits.size(); i++) {
            const Hit &h0 = sorted_hits[i];
            const Hit &h1 = sorted_hits[i + 1];
            const Hit &h2 = sorted_hits[i + 2];

            double dz1 = h1.z - h0.z;
            double dz2 = h2.z - h1.z;

            // On impose une proximité entre les hits en z
            if (dz1 > max_delta_z || dz2 > max_delta_z) continue;

            ProcessCluster(
                clu_feat,
                h0.amplitude, h1.amplitude, h2.amplitude,
                h0.integral, h1.integral, h2.integral,
                y0,
                h0.z, h1.z, h2.z,
                h0.valid, h1.valid, h2.valid,
                3,
                q0_amp_hist, q1_amp_hist, q2_amp_hist,
                q0_int_hist, q1_int_hist, q2_int_hist,
                q1_0_amp_hist, q1_0_int_hist,
                q2_1_amp_hist, q2_1_int_hist
            );
        }
    }
}


void ProcessFile(
    const std::string &file_name,
    TClusterFeatures &clu_feat,
    TH1F *q0_amp_hist, TH1F *q1_amp_hist, TH1F *q2_amp_hist,
    TH1F *q0_int_hist, TH1F *q1_int_hist, TH1F *q2_int_hist,
    TH1F *q1_0_amp_hist, TH1F *q1_0_int_hist,
    TH1F *q2_1_amp_hist, TH1F *q2_1_int_hist
) {
    std::unique_ptr<TFile> file(TFile::Open(file_name.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Erreur à l'ouverture du fichier : " << file_name << std::endl;
        return;
    }

    TTree *tree = dynamic_cast<TTree *>(file->Get("hit_tree"));
    if (!tree) {
        std::cerr << "Arbre 'hit_tree' non trouvé dans le fichier : " << file_name << std::endl;
        return;
    }

    InitTree(tree);

    Long64_t nentries = tree->GetEntries();
    std::cout << "Traitement du fichier : " << file_name << " (" << nentries << " événements)" << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        ProcessEvent(
            clu_feat, tree, i,
            q0_amp_hist, q1_amp_hist, q2_amp_hist,
            q0_int_hist, q1_int_hist, q2_int_hist,
            q1_0_amp_hist, q1_0_int_hist,
            q2_1_amp_hist, q2_1_int_hist
        );
    }

    file->Close();
}




void ProcessCluster(
    TClusterFeatures &clu_feat,
    double q0_amp, double q1_amp, double q2_amp,
    double q0_int, double q1_int, double q2_int,
    double y0,
    double z0, double z1, double z2,
    bool h0, bool h1, bool h2,
    int clu_size,
    TH1F *q0_amp_hist, TH1F *q1_amp_hist, TH1F *q2_amp_hist,
    TH1F *q0_int_hist, TH1F *q1_int_hist, TH1F *q2_int_hist,
    TH1F *q1_0_amp_hist, TH1F *q1_0_int_hist,
    TH1F *q2_1_amp_hist, TH1F *q2_1_int_hist
) {
    if (h0) {
        q0_amp_hist->Fill(q0_amp);
        q0_int_hist->Fill(q0_int);
    }

    if (h1) {
        q1_amp_hist->Fill(q1_amp);
        q1_int_hist->Fill(q1_int);
    }

    if (h2) {
        q2_amp_hist->Fill(q2_amp);
        q2_int_hist->Fill(q2_int);
    }

    if (h0 && h1) {
        q1_0_amp_hist->Fill(q1_amp / q0_amp);
        q1_0_int_hist->Fill(q1_int / q0_int);
    }

    if (h2 && h1) {
        q2_1_amp_hist->Fill(q2_amp / q1_amp);
        q2_1_int_hist->Fill(q2_int / q1_int);
    }

    // Stockage dans les vecteurs de la structure clu_feat
    clu_feat.q.push_back(TMath::Log(q1_amp / q0_amp));
    clu_feat.ori.push_back(0);
    clu_feat.y.push_back(y0);

    if (use_integral) {
        clu_feat.z.push_back(GetBarycenter_Z(z0, z1, z2, q0_int, q1_int, q2_int, h0, h1, h2));
        clu_feat.err.push_back(GetBarycenter_Z_Error(z0, z1, z2, q0_int, q1_int, q2_int, h0, h1, h2));
    } else {
        clu_feat.z.push_back(GetBarycenter_Z(z0, z1, z2, q0_amp, q1_amp, q2_amp, h0, h1, h2));
        clu_feat.err.push_back(GetBarycenter_Z_Error(z0, z1, z2, q0_amp, q1_amp, q2_amp, h0, h1, h2));
    }

    clu_feat.size.push_back(clu_size);
}





int matinv3(double *a, double *b){
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
int FitCirc(int np, double *z, double *y, double *sig, double *pari, double *parf, double &chi2,
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
int FitParab(int nz, double *z, double *y, double *sigy, double *parf, double &chi2, double *C)
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
void Shape(int np, double *z, double *y, double &phi, double *zp, double *up, double &length,
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
int FitHelix(int np, double *z, double *y, int *dirm, double *sig, double *zp, double *up,
                                 double &zfit, double &yfit, double &phifit, double &Rfit, double *param, double &chi2,
                                 double *cov, bool forceCircFit)

{
   if (np == 0) {
      std::cout << "No points in the track" << std::endl;
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
   std::cout << "Width of pattern: " << length << " <- " << LENGTH_MIN << std::endl;

   if (length < LENGTH_MIN)
   return 0;
   cphi = cos(phi);
   sphi = sin(phi);

   // if(width/length<WIDTH_OVER_LENGTH_MAX and not forceCircFit) {

  std::cout << "Parabolic Fit" << std::endl; 
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
      std::cout << "Parabolic Fit Failed" << std::endl;
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

   std::cout << "parabola output: p1 = " << param[0]
          << ", p2 = " << param[1]
          << ", p3 = " << param[2] << std::endl;


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

   if (std::abs(param[2]) < 1e-3 and parab_fit_res == 1) {
      std::cout << "INFO: Parabolic fit chosen because of small curvature" << std::endl;
      // return 1;
   }

   std::cout << "INFO: Circular fit" << std::endl;
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
      std::cout << "INFO: Circular fit failed because of det=0" << std::endl;
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
      std::cerr << "WARNING: Circular fit failed because of calculation" << std::endl;
      if (parab_fit_res == 1) {
         std::cout << "Parabolic fit was successful; using it instead" << std::endl;
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

   std::cout << "Circular fit done successfully" << std::endl;


   return 2;
}











//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double integrale(double val_lim, TGraph *graph) {
    double integral = 0.0;
    int n = graph->GetN();

    for (int i = 0; i < n - 1; ++i) {
        double t1, v1, t2, v2;
        graph -> GetPoint(i,     t1, v1);
        graph -> GetPoint(i + 1, t2, v2);

        double dt = t2-t1;
        if (v1 > val_lim) {
            integral += v1*dt;
        }
    }
    return integral;
}



double coordxpic(TGraph *graph){
    double max_val = -1000;
    double t_at_max = 0;

    int n_points = graph -> GetN();
        for ( int i = 0; i < n_points; ++i ) {
            double t;
            double val;
            graph -> GetPoint(i, t, val);
            if (val > max_val) {
                        max_val = val;
                        t_at_max = t;
                    }
                }
                return t_at_max;  // Temps du pic maximal
            }



double coordypic(TGraph *graph){
    double max_val = -1000;

    int n_points = graph -> GetN();
        for ( int i = 0; i < n_points; ++i ) {
            double t;
            double val;
            graph -> GetPoint(i, t, val);
            if (val > max_val) {
                        max_val = val;
            }
        }

    return max_val;  // Hauteur du pic maximale
    }



void histo(TTree *tree, int hit_filter,TH1F *h_charge, TH1F *h_peak_time, TH1F *h_amplitude) {
    std::vector<double> charges;
    std::vector<double> peak_times;
    std::vector<double> max_amplitudes;

    Int_t cluster_per_trace;
    tree -> SetBranchAddress("trk_ncl", &cluster_per_trace);

    int hit;
    TGraph *wf = nullptr;

    tree -> SetBranchAddress("hit", &hit);
    tree -> SetBranchAddress("WF", &wf);

    int nentries = tree -> GetEntries();
    for (int i = 0; i < nentries && i < N_hits; ++i) {
        tree->GetEntry(i);
        if (cluster_per_trace > bord){

            if (hit == hit_filter) {
                double charge = integrale(10, wf);
                double v_max = coordypic(wf);
                double t_at_max = coordxpic(wf);

                charges.push_back(charge);
                peak_times.push_back(t_at_max);
                max_amplitudes.push_back(v_max);
            }
        }

    }


    for (size_t i = 0; i < charges.size(); ++i) {
        h_charge -> Fill(charges[i]);
        h_peak_time -> Fill(peak_times[i]);
        h_amplitude -> Fill(max_amplitudes[i]);
    }

    //On affiche dans des canva séparés
    TCanvas *c_charge = new TCanvas(Form("c_charge_hit%d", hit_filter), Form("Charge (hit=%d)", hit_filter));
    h_charge -> Draw();

    TCanvas *c_peak_time = new TCanvas(Form("c_peak_time_hit%d", hit_filter), Form("Temps du pic max (hit=%d)", hit_filter));
    h_peak_time -> Draw();

    TCanvas *c_amplitude = new TCanvas(Form("c_amplitude_hit%d", hit_filter), Form("Amplitude max (hit=%d)", hit_filter));
    h_amplitude -> Draw();

    c_charge -> Update();
    c_peak_time ->  Update();
    c_amplitude -> Update();
}





    // Curve fit pour approximer les courbes à une landau + gauss

    double gaussPlusLandau(Double_t *x, Double_t *par) {
        //[0] = moy gaussienne, [1] = amp gaussienne,[2] = pic Landau,[3] = largeur Landau, [4] = amp Landau
        double gauss = par[1] * exp(- pow(x[0] - par[0], 2));
        double lambda = (x[0] - par[2]) / par[3];
        double landau = par[4] * (1.0 / (2 * TMath::Pi())) * exp(-0.5 * (lambda + exp(-lambda)));
        return (gauss + landau);
    }

    double Gauss(Double_t *x, Double_t *par) {
        //[0] = moy, [1]= amplitude
        double arg = x[0] - par[0];
        return (par[1] * exp(- arg * arg));
    }

    double Landau(Double_t *x, Double_t *par) {
        //[0] =abscisse du pic, [1]=largeur à mi hauteur ,[2] = amp
        double lambda = (x[0] - par[0]) / par[1];
        double val = (1.0 / sqrt(2 * TMath::Pi())) * exp(-0.5 * (lambda + exp(-lambda)));
        return (par[2] * val);
    }



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_wf(int argc, char** argv) {


    std::vector<double> Q0_s_Q1_c;
    std::vector<double> Q1_c;
    std::vector<double> Q1_s_Q2_c;
    std::vector<double> Q2_c;


    //On utilise ces mêmes vecteurs avec un a (amplitude : objectif de comparer charger et amplitude
    std::vector<double> Q0_s_Q1_a;
    std::vector<double> Q1_a;
    std::vector<double> Q1_s_Q2_a;
    std::vector<double> Q2_a;

    std::vector<double> z_hit;
    std::vector<double> z_clu;
    std::vector<double> ln_Ratio_ZERO_c;
    std::vector<double> ln_Ratio_ZERO_a;

    //Vceteurs à remplir pour comparer l'amplitue et la charge dues aux cluster avec un seul hit (0)
    std::vector<double> AMPLITUDE_1;
    std::vector<double> AMPLITUDE_PLUS_QUE_1;
    std::vector<double> CHARGE_1;
    std::vector<double> CHARGE_PLUS_QUE_1;


    //Vecteur à remplir pour histo des clusters par trace
    std::vector<int> cluster_par_trace;

    //Vecteur à remplir pour avoir un histogramme de hits par cluster
    std::vector<int> hit_average;

    const std::string inputFile = argv[1];
    const std::string outputFile = argv[2];

    // Reading input
    TFile *myf = TFile::Open(inputFile.c_str());
    if(!myf)
    {
        std::cout << " ERROR: input file not found" << std::endl;
        throw;
    }

    TTree *tree = (TTree *) myf->Get("outTree");
    if(!tree)
    {
        std::cout << " ERROR: input TTree not found" << std::endl;
        throw;
    }

    // Declaration of leaf types
    Int_t event;
    Int_t pattern;
    Int_t cluster;
    Int_t hit;
    Int_t cluster_per_trace;
    Double_t clu_pos[3];


    Int_t    trk_ncl;
    Double_t trk_pos[3];
    Double_t trk_dir[3];
    Int_t    clu_size;
    TGraph *wf = 0;

    // List of branches
    TBranch *b_event;   //!
    TBranch *b_pattern;   //!
    TBranch *b_cluster;   //!
    TBranch *b_hit;   //!
    TBranch *b_trk_pos;   //!
    TBranch *b_clu_size;   //!
    TBranch *b_wf;   //!
    TBranch *b_cluster_per_trace;  //!
    TBranch *b_trk_dir; //!
    TBranch *b_clu_pos; //!

    // Linking branches to variables
    tree->SetBranchAddress("event",    &event, &b_event);
    tree->SetBranchAddress("pattern",  &pattern, &b_pattern);
    tree->SetBranchAddress("cluster",  &cluster, &b_cluster);
    tree->SetBranchAddress("hit",      &hit, &b_hit);
    tree->SetBranchAddress("trk_pos",  &trk_pos, &b_trk_pos);
    tree->SetBranchAddress("clu_size", &clu_size, &b_clu_size);
    tree->SetBranchAddress("WF",       &wf);
    tree->SetBranchAddress("trk_ncl",  &cluster_per_trace, &b_cluster_per_trace);
    tree->SetBranchAddress("trk_dir",  &trk_dir, &b_trk_dir);
    tree->SetBranchAddress("clu_pos",  &clu_pos, &b_clu_pos);


    Int_t nentries = tree->GetEntries();
    std::cout << "entries = " << nentries << std::endl;


    TGraph * wf_1 = new TGraph();
    TGraph * wf_2 = new TGraph();
    TGraph * wf_3 = new TGraph();
    TGraph * wf_4 = new TGraph();

    Int_t look_event   = 0;
    Int_t look_pattern = 0;
    Int_t look_cluster = 7;

    std::vector<int> ZERO;
    std::vector<int> UNO;

    std::vector<double> av_charge_tot;
    std::vector<double> av_amp_tot;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int temp_cluster_per_trace;

    tree -> GetEntry(0);
    temp_cluster_per_trace = cluster_per_trace;
    std::vector<int> hit_count;
    std::vector<int> hit_count_temp;

    int x = 0;
    int indice = 0;

    std::vector<std::vector<double>> zp_all;
    std::vector<std::vector<double>> up_all;
    int indice_t = 0;
    int indice_cl = 0;

    double z_av_c = 0;
    double z_av_a = 0;
    double temp_c_av = 0;
    double temp_a_av = 0;


    // Reading TTree and filling data of interest
    for (int i = 0; i < nentries && i < N_hits ; ++i) {

        tree -> GetEntry(i);



        if (trk_dir[1] > 0.9){  //Pour la direction verticqle des particules

        if (fabs (trk_pos[0]) > 900){

        if (cluster_per_trace > bord){  //On regarde si on a bien un nombre de cluster assez grand qui signifie que l'on a bien un trace est pas du bruit

        if (temp_cluster_per_trace != cluster_per_trace){  //on regarde si on est dans la même trace
            cluster_par_trace.push_back(temp_cluster_per_trace);
        }
        temp_cluster_per_trace = cluster_per_trace;


        //pour faire Q1/Q0 en amplitude et en charge
        if (hit == 0 && clu_size > 1){
            Q0_s_Q1_c.push_back(integrale(10, wf));
            Q0_s_Q1_a.push_back(coordypic(wf));
            z_clu.push_back(clu_pos[2]);
        }

        if (hit == 1){
            Q1_c.push_back(integrale(10, wf));
            Q1_a.push_back(coordypic(wf));
        }

        //Pour faire Q2/Q1 en amplitude et en charge
        if (hit == 1 && clu_size > 2){
            Q1_s_Q2_c.push_back(integrale(10, wf));
            Q1_s_Q2_a.push_back(coordypic(wf));
        }

        if (hit == 2){
            Q2_c.push_back(integrale(10, wf));
            Q2_a.push_back(coordypic(wf));
        }

        //Pour "split" les différentes composantes des hits 0 :celles oùil n'y a qu'un hit par cluster ou plus
        if (clu_size == 1 && hit == 0){
           AMPLITUDE_1.push_back(coordypic(wf));
           CHARGE_1.push_back(integrale(10,wf));
        }
        if (clu_size > 1 && hit == 0){
           AMPLITUDE_PLUS_QUE_1.push_back(coordypic(wf));
           CHARGE_PLUS_QUE_1.push_back(integrale(10,wf));
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Looking at a single track

        //if (event == look_event && pattern == look_pattern && cluster == look_cluster) {
            if (hit == 0){
                hit_count.push_back(hit_count_temp.size());
            	hit_count_temp.clear();
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_1   -> SetPoint(j, t, q);
		}
            }

            if (hit == 1){
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_2   -> SetPoint(j, t, q);
                }
            }
            if (hit == 2){
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_3   -> SetPoint(j, t, q);
                }
            }
            if (hit == 3){
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_4   -> SetPoint(j, t, q);
                }
            }
	    hit_count_temp.push_back(hit);

        //}
    }
    }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Plotting

    //Pour les ratios de Q1/Q0 et Q2/Q1
    std::vector<double> ratio_ZERO_c;
    std::vector<double> ratio_UNO_c;

    for (size_t i = 0; i < Q0_s_Q1_c.size(); ++i) {
            ratio_ZERO_c.push_back((double)Q1_c[i] / Q0_s_Q1_c[i]);
        }

    for (size_t i = 0; i < Q1_s_Q2_c.size(); ++i) {
            ratio_UNO_c.push_back((double)Q2_c[i] / Q1_s_Q2_c[i]);
        }


    std::vector<double> ratio_ZERO_a;
    std::vector<double> ratio_UNO_a;

    for (size_t i = 0; i < Q0_s_Q1_a.size(); ++i) {
        ratio_ZERO_a.push_back((double)Q1_a[i] / Q0_s_Q1_a[i]);
    }

    for (size_t i = 0; i < Q1_s_Q2_c.size(); ++i) {
        ratio_UNO_a.push_back((double)Q2_a[i] / Q1_s_Q2_a[i]);
    }


    TCanvas * myc = new TCanvas();
    TGraph * wf_0 = new TGraph();

    wf_0 -> SetPoint(0,0,-100);
    wf_0 -> SetPoint(1,511,2e3);
    wf_0 -> SetTitle(Form("event %i, pattern %i, cluster %i;t (time stamps); ADC", look_event, look_pattern, look_cluster));
    wf_0 -> Draw("AP");

    if( wf_1 -> GetN() > 0 ){
        wf_1 -> SetLineColor(1);
        wf_1 -> Draw("L same");
    }
    if( wf_2 -> GetN() > 0 ){
        wf_2-> SetLineColor(2);
        wf_2 -> Draw("L same");
    }
    if( wf_3 -> GetN() > 0 ){
        wf_3 -> SetLineColor(3);
        wf_3 -> Draw("L same");
    }
    if( wf_4 -> GetN() > 0 ){
        wf_4 -> SetLineColor(4);
        wf_4 -> Draw("L same");
    }


    TLegend* myl = new TLegend(); 


    if( wf_1 -> GetN() > 0 )
        myl -> AddEntry(wf_1 , "hit_0");
    if( wf_2 -> GetN() > 0 )
        myl -> AddEntry(wf_2 , "hit_1");
    if( wf_3 -> GetN() > 0 )
        myl -> AddEntry(wf_3 , "hit_2");
    if( wf_4 -> GetN() > 0 )
        myl -> AddEntry(wf_4 , "hit_3");

    myl -> Draw("same");


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Histogrammes des charges
    TFile * outfile = new TFile(outputFile.c_str(),"recreate");
    outfile -> cd();

    for (int hit_filter = 0; hit_filter <=3; ++hit_filter) {
        // Création des histogrammes
        TH1F *h_charge    = new TH1F(Form("h_charge_hit%d",    hit_filter), Form("Distribution des charges (hit=%d);Charge;Nombre de hits", hit_filter), bin_charge, 0, X_charge);
        TH1F *h_peak_time = new TH1F(Form("h_peak_time_hit%d", hit_filter), Form("Temps du pic max (hit=%d);Temps (ns);Nombre de hits",     hit_filter), 50, 0, 511);
        TH1F *h_amplitude = new TH1F(Form("h_amplitude_hit%d", hit_filter), Form("Amplitude max (hit=%d);Amplitude (ADC);Nombre de hits",   hit_filter), bin_amplitude, 0, X_amplitude);

        //on remplit les histogrammes
        histo(tree, hit_filter, h_charge, h_peak_time, h_amplitude);

        h_charge ->   Write();
        h_peak_time -> Write();
        h_amplitude -> Write();

        //On fait 1 canva par histogramme
        TCanvas *c_charge = new TCanvas(Form("c_charge_hit%d", hit_filter), Form("Charge hit %d", hit_filter)) ;
        h_charge -> Draw();
        c_charge -> Update();

        TCanvas *c_peak_time = new TCanvas(Form("c_peak_time_hit%d", hit_filter), Form("Peak time hit %d", hit_filter)) ;
        h_peak_time -> Draw();
        c_peak_time -> Update();

        TCanvas *c_amplitude = new TCanvas(Form("c_amplitude_hit%d", hit_filter), Form("Amplitude hit %d", hit_filter)) ;
        h_amplitude -> Draw();
        c_amplitude -> Update();

        }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Superposition des histogrammes de charge
    TCanvas *c_all_charges = new TCanvas("c_all_charges");

    TH1F *h0 = (TH1F*) gDirectory -> Get("h_charge_hit0");
    TH1F *h1 = (TH1F*) gDirectory -> Get("h_charge_hit1");
    TH1F *h2 = (TH1F*) gDirectory -> Get("h_charge_hit2");
    TH1F *h3 = (TH1F*) gDirectory -> Get("h_charge_hit3");

    h0 -> SetLineColor(1);
    h1 -> SetLineColor(2);
    h2 -> SetLineColor(3);
    h3 -> SetLineColor(4);

    float max_val = std::max({h0 -> GetMaximum(), h1 -> GetMaximum(), h2 -> GetMaximum(), h3 -> GetMaximum()});
    h0 -> SetMaximum(1.1 * max_val);


    gStyle->SetOptStat(0);


    h0 -> Draw();
    h1 -> Draw("same");
    h2 -> Draw("same");
    h3 -> Draw("same");

   TIter next(c_all_charges->GetListOfPrimitives());
	TObject *obj;
	while ((obj = next())) {
    		if (obj->InheritsFrom("TPaveStats")) {
        	obj->Delete();
    }
}

///////////////////////////////////////////////////////////////////

    //Fit gauss + landau

    TF1 *fitFunc0 = new TF1("fitFunc0", gaussPlusLandau, 0, 18000, 5);
    TF1 *fitFunc1 = new TF1("fitFunc1", gaussPlusLandau, 0, 18000, 5);
    TF1 *fitFunc2 = new TF1("fitFunc2", gaussPlusLandau, 0, 18000, 5);
    TF1 *fitFunc3 = new TF1("fitFunc3", gaussPlusLandau, 0, 18000, 5);

    TF1 *fitFunc0B = new TF1("fitFunc0B", Landau, 0, 18000, 3);
    TF1 *fitFunc1B = new TF1("fitFunc1B", Landau, 0, 18000, 3);
    TF1 *fitFunc2B = new TF1("fitFunc2B", Landau, 0, 18000, 3);
    TF1 *fitFunc3B = new TF1("fitFunc3B", Landau, 0, 18000, 3);


    fitFunc0 -> SetParameters(6000, 85, 6000, 2000, 765);
    fitFunc1 -> SetParameters(3000, 150, 3000, 3000, 1350);
    fitFunc2 -> SetParameters(900, 700, 2600, 800, 100);
    fitFunc3 -> SetParameters(900, 120, 900, 1000, 100);

    //fitFunc0B -> SetParameters(6200,6000, 850);
    //fitFunc1B -> SetParameters(3000,3000, 1500);
    //fitFunc2B -> SetParameters(2100,1000, 1200);
    //fitFunc3B -> SetParameters();


    fitFunc0 -> SetLineColor(5);
    fitFunc1 -> SetLineColor(6);
    fitFunc2 -> SetLineColor(7);
    fitFunc3 -> SetLineColor(8);

    //fitFunc0B -> SetLineColor(9);
    //fitFunc1B -> SetLineColor(1);
    //fitFunc2B -> SetLineColor(2);
    //fitFunc3B -> SetLineColor(3);

    fitFunc2 -> FixParameter(0, 900);
    //fitFunc2 -> FixParameter(1, 1200);
    //fitFunc2 -> FixParameter(2, 2000);


    h0 -> Fit(fitFunc0, "R");
    h1 -> Fit(fitFunc1, "R");
    h2 -> Fit(fitFunc2, "R");
    h3 -> Fit(fitFunc3, "R");

    //h0 -> Fit(fitFunc0B, "R");
    //h1 -> Fit(fitFunc1B, "R");
    //h2 -> Fit(fitFunc2B, "R");
    //h3 -> Fit(fitFunc3B, "R");


    fitFunc0 -> Draw("same");
    fitFunc1 -> Draw("same");
    fitFunc2 -> Draw("same");
    fitFunc3 -> Draw("same");

    //fitFunc0B -> Draw("same");
    //fitFunc1B -> Draw("same");
    //fitFunc2B -> Draw("same");
    //fitFunc3B -> Draw("same");


    fitFunc0 -> Write();
    fitFunc1 -> Write();
    fitFunc2 -> Write();
    fitFunc3 -> Write();

    //fitFunc0B -> Write();
    //fitFunc1B -> Write();
    //fitFunc2B -> Write();
    //fitFunc3B -> Write();



    TLegend *leg = new TLegend(0.8, 0.7, 0.9, 0.9);
    leg -> AddEntry(h0, "Hit 0", "l");
    leg -> AddEntry(h1, "Hit 1", "l");
    leg -> AddEntry(h2, "Hit 2", "l");
    leg -> AddEntry(h3, "Hit 3", "l");
    leg -> Draw();


    //On veut afficher les entries sur ce même histo
    TPaveText *infos_entries = new TPaveText(0.35, 0.75, 0.6, 0.9, "NDC");

    infos_entries -> AddText(Form("Entries hit 0 : %d", (int)h0 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 1 : %d", (int)h1 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 2 : %d", (int)h2 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 3 : %d", (int)h3 -> GetEntries()));


    infos_entries -> Draw();
    c_all_charges -> Modified();
    c_all_charges -> Update();
    c_all_charges -> Write();

    outfile -> cd();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //On travaille en amplitude
    //Histogramme QO_A sachant qu'il existe Q1
    TH1F *h_Q0_s_Q1B = new TH1F("Q0_sachant_Q1_A", "Q0_sachant_Q1_A", bin_amplitude, 0, X_amplitude);
    for (double r : Q0_s_Q1_a) h_Q0_s_Q1B -> Fill(r);
    h_Q0_s_Q1B -> Write();

    //Histogramme de Q1
    TH1F *h_Q1B = new TH1F("Q1_A", "Q1_A", bin_amplitude, 0, X_amplitude);
    for (double r : Q1_a) h_Q1B -> Fill(r);
    h_Q1B -> Write();

    //Histogramme Q1_A sachant qu'il existe Q2
    TH1F *h_Q1_s_Q2B = new TH1F("Q1_sachant_Q2_A", "Q1_sachant_Q2_A", bin_amplitude, 0, X_amplitude);
    for (double r : Q1_s_Q2_a) h_Q1_s_Q2B -> Fill(r);
    h_Q1_s_Q2B-> Write();

    //Histogramme de Q2
    TH1F *h_Q2B = new TH1F("Q2_A", "Q2_A", bin_amplitude, 0, X_amplitude);
    for (double r : Q2_a) h_Q2B -> Fill(r);
    h_Q2B -> Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //On travaille en charge
    //Histogramme QO_c sachant qu'il existe Q1
    TH1F *h_Q0_s_Q1 = new TH1F("Q0_sachant_Q1_C", "Q0_sachant_Q1_C", bin_charge, 0, X_charge);
    for (double r : Q0_s_Q1_c) h_Q0_s_Q1 -> Fill(r);
    h_Q0_s_Q1 -> Write();

    //Histogramme de Q1
    TH1F *h_Q1 = new TH1F("Q1_C", "Q1_C", bin_charge, 0, X_charge);
    for (double r : Q1_c) h_Q1 -> Fill(r);
    h_Q1 -> Write();

    //Histogramme Q1_A sachant qu'il existe Q2
    TH1F *h_Q1_s_Q2 = new TH1F("Q1_sachant_Q2_C", "Q1_sachant_Q2_C", bin_charge, 0, X_charge);
    for (double r : Q1_s_Q2_c) h_Q1_s_Q2 -> Fill(r);
    h_Q1_s_Q2-> Write();

    //Histogramme de Q2
    TH1F *h_Q2 = new TH1F("Q2_C", "Q2_C", bin_charge, 0, X_charge);
    for (double r : Q2_c) h_Q2 -> Fill(r);
    h_Q2 -> Write();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Q0 en fonction de Q0 pour amplitude et charge sachant qu'il existe Q1

    TGraph *gr0 = new TGraph(Q0_s_Q1_a.size());
    for (size_t i = 0; i < Q0_s_Q1_a.size(); ++i) {
        gr0 -> SetPoint(i, Q0_s_Q1_a[i], Q0_s_Q1_c[i]); //charge en fonction de l'amplitude
    }

    gr0 -> SetTitle("Q0_c en fonction de Q0_a sachant qu'il existe Q1; Amplitude;Charge");
    gr0 -> SetMarkerStyle(20);
    gr0 -> SetLineWidth(0);
    gr0 -> SetMarkerSize(0.4);
    gr0 -> SetName("graph_Q0_vs_Q1");
    gr0 -> Draw("AP");
    gr0 -> Write();
////////////////////////////////////////////////////////////

    //Q1 en fonction de Q1 pour amplitude et charge sachant qu'il existe Q2

    TGraph *gr1 = new TGraph(Q1_s_Q2_a.size());
    for (size_t i = 0; i < Q1_s_Q2_a.size(); ++i) {
        gr1 -> SetPoint(i, Q1_s_Q2_a[i], Q1_s_Q2_c[i]); //charge en fonction de l'amplitude
    }

    gr1 -> SetTitle("Q1_c en fonction de Q1_a sachant qu'il existe Q2; Amplitude;Charge");
    gr1 -> SetMarkerStyle(20);
    gr1 -> SetLineWidth(0);
    gr1 -> SetMarkerSize(0.4);
    gr1 -> SetName("graph_Q1_c_vs_Q1_a_BIS");
    gr1 -> Draw("AP");
    gr1 ->Write();

////////////////////////////////////////////////////////////

    //Q1 en fonction de Q1 pour amplitude et charge

    TGraph *gr2 = new TGraph(Q1_a.size());
    for (size_t i = 0; i < Q1_a.size(); ++i) {
        gr2 -> SetPoint(i, Q1_a[i], Q1_c[i]); //charge en fonction de l'amplitude
    }

    gr2 -> SetTitle("Q1_c en fonction de Q1_a; Amplitude;Charge");
    gr2 -> SetMarkerStyle(20);
    gr2 -> SetLineWidth(0);
    gr2 -> SetMarkerSize(0.4);
    gr2 -> SetName("graph_Q1_c_vs_Q1_a");
    gr2 -> Draw("AP");
    gr2 ->Write();
////////////////////////////////////////////////////////////

    //Q2 en fonction de Q2 pour amplitude et charge

    TGraph *gr3 = new TGraph(Q2_a.size());
    for (size_t i = 0; i < Q2_a.size(); ++i) {
        gr3 -> SetPoint(i, Q2_a[i], Q2_c[i]); //charge en fonction de l'amplitude
    }

    gr3 -> SetTitle("Q2_c en fonction de Q2_a; Amplitude;Charge");
    gr3 -> SetMarkerStyle(20);
    gr3 -> SetLineWidth(0);
    gr3 -> SetMarkerSize(0.4);
    gr3 -> SetName("graph_Q2_c_vs_Q2_a");
    gr3 -> Draw("AP");
    gr3 ->Write();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Histo du nombre de hits par cluster
    TH1F *hist = new TH1F("hit_per_cluster", "Histo des hit/cluster;Hit N;Nombre de cluster avec N hits", 6, 0, 6);
    for (int val : hit_count) {
        hist -> Fill(val);
    }

    hist->Draw("hist");

    TLegend *legici = new TLegend(0.6, 0.7, 0.9, 0.85);
    legici -> AddEntry(hist, Form("N_{elements} = %d", (int)hit_count.size()), "f");
    legici -> Draw();

    hist->Write();

    //Histo du nombre de cluster par trace
    TH1F *h_cluster_par_trace = new TH1F("cluster_per_trace", "Clusters per trace", 100, 0, 100);
    for (double r : cluster_par_trace) h_cluster_par_trace -> Fill(r);
    h_cluster_par_trace -> Write();



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Histo des ratios Q1/Q0 et Q2/Q1 en amplitude et charge
    TH1F *h_ratio_ZERO = new TH1F("Q1_sur_Q0_charge", "Ratio Q1/Q0;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_ZERO_c) h_ratio_ZERO -> Fill(r);
    h_ratio_ZERO -> Write();



//////////////////////////////////////////////////////

    TH1F *h_ratio_UNO = new TH1F("Q2_sur_Q1_charge", "Ratio Q2/Q1;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_UNO_c) h_ratio_UNO -> Fill(r);
    h_ratio_UNO -> Write();



////////////////////////////////////////////////////////

    //Histo du ratio des amplitudes
    TH1F *h_ratio_ZEROB = new TH1F("Q1_sur_Q0_amplitude", "Ratio Q1/Q0;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_ZERO_a) h_ratio_ZEROB -> Fill(r);
    h_ratio_ZEROB -> Write();



////////////////////////////////////////////////////////

    TH1F *h_ratio_UNOB = new TH1F("Q2_sur_Q1_amplitude", "Ratio Q2/Q1;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_UNO_a) h_ratio_UNOB -> Fill(r);
    h_ratio_UNOB -> Write();


////////////////////////////////////////////////////////

    int nPoints = z_clu.size();
    for (double i : ratio_ZERO_c) {
            ln_Ratio_ZERO_c.push_back(std::log(i));

    }


    std::cout << "z_clu.size() = " << z_clu.size() << std::endl;
    std::cout << "ln_Ratio_ZERO_c.size() = " << ln_Ratio_ZERO_c.size() << std::endl;


    TGraph *gr4 = new TGraph(nPoints, &z_clu[0], &ln_Ratio_ZERO_c[0]);

    gr4 -> SetTitle("Ratio ZERO_c vs Z_clu;Z;lnQ1/Q0_c");
    gr4 -> SetMarkerStyle(20);
    gr4 -> SetMarkerColor(1);
    gr4 -> SetLineWidth(0);
    gr4 -> SetMarkerSize(0.4);
    gr4 -> SetName("ln_ratio_vs_z_c");

    gr4 -> Write();

    // Création d’un canvas et affichage
    TCanvas *patatoide_c = new TCanvas("patatoide_c", "patatoide_c");

    gr4 -> Draw("AP");
    patatoide_c -> Update();
    patatoide_c -> Write();

///////////////////////////////////////////////////////


    for (double i : ratio_ZERO_a) {
            ln_Ratio_ZERO_a.push_back(std::log(i));

    }


    std::cout << "z_clu.size() = " << z_clu.size() << std::endl;
    std::cout << "ln_Ratio_ZERO_c.size() = " << ln_Ratio_ZERO_a.size() << std::endl;


    TGraph *gr5 = new TGraph(nPoints, &z_clu[0], &ln_Ratio_ZERO_a[0]);

    gr5 -> SetTitle("Ratio ZERO_a vs Z_clu;Z;lnQ1/Q0_a");
    gr5 -> SetMarkerStyle(20);
    gr5 -> SetMarkerColor(1);
    gr5 -> SetLineWidth(0);
    gr5 -> SetMarkerSize(0.4);
    gr5 -> SetName("ln_ratio_vs_z_a");

    gr5 -> Write();

    // Création d’un canvas et affichage
    TCanvas *patatoide_a = new TCanvas("patatoide_a", "patatoide_a");

    gr5 -> Draw("AP");
    patatoide_a -> Update();
    patatoide_a -> Write();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Superposition des histogrammes dans un Canvas
    TCanvas *c_ratio_01 = new TCanvas("c_ratio_01");
    h_ratio_ZEROB -> SetLineColor(1);
    h_ratio_ZERO ->  SetLineColor(2);
    double max1 = h_ratio_ZEROB -> GetMaximum();
    double max2 = h_ratio_ZERO -> GetMaximum();
    double ymax = std::max( max1, max2) * 1.1;
    h_ratio_ZEROB  -> SetMaximum(ymax);
    h_ratio_ZERO -> SetMaximum(ymax);

    h_ratio_ZEROB -> Draw();
    h_ratio_ZERO -> Draw("SAME");
    TLegend *leg1 = new TLegend(0.65, 0.7, 0.9, 0.9);

    TString label1 = Form("Q1 / Q0 Amplitude (%.0f entries)", h_ratio_ZEROB->GetEntries());
    TString label2 = Form("Q1 / Q0 Charge (%.0f entries)", h_ratio_ZERO->GetEntries());

    leg1->AddEntry(h_ratio_ZEROB, label1, "l");
    leg1->AddEntry(h_ratio_ZERO,  label2, "l");

    leg1 -> Draw();
    c_ratio_01 -> Write();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c_ratio_02 = new TCanvas("c_ratio_02");
    h_ratio_UNOB -> SetLineColor(1);
    h_ratio_UNO -> SetLineColor(2);
    double max3 = h_ratio_UNOB -> GetMaximum();
    double max4 = h_ratio_UNO -> GetMaximum();
    double ymax2 = std::max( max3, max4) * 1.1;
    h_ratio_UNOB -> SetMaximum(ymax2);
    h_ratio_UNO -> SetMaximum(ymax2);

    h_ratio_UNOB -> Draw();
    h_ratio_UNO -> Draw("SAME");
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2 -> AddEntry(h_ratio_UNOB, " Q2 / Q1 Amplitude", "l");
    leg2 -> AddEntry(h_ratio_UNO, " Q2 / Q1 Charge", "l");
    leg2 -> Draw();
    c_ratio_02 -> Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH1F* h_charge_1 = new TH1F("h_charge_1", "Charge comparaison;Charge;Nombre de hits", bin_charge, 0, X_charge);
    for (double val : CHARGE_1) h_charge_1->Fill(val);

    TH1F* h_charge_plus = new TH1F("h_charge_plus", "Charge comparaison;Charge;Nombre de hits", bin_charge, 0, X_charge);
    for (double val : CHARGE_PLUS_QUE_1) h_charge_plus->Fill(val);

    TCanvas* c_charge = new TCanvas("c_charge_split");

    h_charge_1 -> SetLineColor(1);
    h_charge_plus -> SetLineColor(2);

    double maxc1 = h_charge_1 -> GetMaximum();
    double maxc2 = h_charge_plus -> GetMaximum();
    double ymax3 = std::max( maxc1, maxc2) * 1.1;
    h_charge_1 -> SetMaximum(ymax3);
    h_charge_plus -> SetMaximum(ymax3);
    h_charge_1->GetXaxis()->SetRangeUser(0, X_charge);
    h_charge_1 -> Draw();
    h_charge_plus -> Draw("SAME");

    int n1_c = h_charge_1->GetEntries();
    int nplus_c = h_charge_plus->GetEntries();

    TString label1_c = Form("clusters avec 1 hit (n = %d)", n1_c);
    TString label2_c = Form("clusters avec >1 hit (n = %d)", nplus_c);


    TLegend *leg_charge = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg_charge->AddEntry(h_charge_1, label1_c, "l");
    leg_charge->AddEntry(h_charge_plus, label2_c, "l");
    leg_charge -> Draw();
    c_charge -> Draw();
    c_charge -> Write();
    c_charge -> Update();


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    TH1F * h_amp_1 = new TH1F("h_amp_1", "Amplitude comparison;Amplitude;Entries", bin_amplitude, 0, X_amplitude);
    TH1F * h_amp_plus = new TH1F("h_amp_plus", "Amplitude comparison;Amplitude;Entries", bin_amplitude, 0, X_amplitude);

    for (double val : AMPLITUDE_1)
        h_amp_1 -> Fill(val);

    for (double val : AMPLITUDE_PLUS_QUE_1)
        h_amp_plus -> Fill(val);

    TCanvas * c_amp = new TCanvas("c_amp_split", "Amplitude per cluster size");
    h_amp_1 -> SetLineColor(3);
    h_amp_plus -> SetLineColor(4);
    gStyle -> SetOptStat(0);


    double maxa1 = h_amp_1 -> GetMaximum();
    double maxa2 = h_amp_plus -> GetMaximum();
    double ymax4 = std::max( maxa1, maxa2) * 1.1;
    h_amp_1 -> SetMaximum(ymax4);
    h_amp_plus -> SetMaximum(ymax4);

    h_amp_1->GetXaxis()->SetRangeUser(0, X_amplitude);
    h_amp_1 -> Draw();
    h_amp_plus -> Draw("SAME");

    int n1_a = h_amp_1->GetEntries();
    int nplus_a = h_amp_plus->GetEntries();

    TString label1_a = Form("clusters avec 1 hit (n = %d)", n1_a);
    TString label2_a = Form("clusters avec >1 hit (n = %d)", nplus_a);

    TLegend * leg_amp = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg_amp->AddEntry(h_amp_1, label1_a, "l");
    leg_amp->AddEntry(h_amp_plus, label2_a, "l");
    leg_amp -> Draw();
    c_amp -> Draw();
    c_amp -> Write();
    c_amp -> Update();





    //data export
    outfile -> cd();
    myc -> Write();
    outfile -> Close();



}
// Utiliser la charge et amplitude pour calculer la position du cluster
// en log Q on utilise l'amplitude
// cool de faire de même mais que avec la charge


//on va comparer les 2 méthodes
//on peut avoir la position du cluster
//on calculer larésolution = ecart type position du cluster et position trace reconstruite

//quelle est la meilleure méthode ?



//puis on calculeera l'écart type


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_root_file>\n";
        return 1; }
    plot_wf(argc, argv);
   // analyse_data(argc, argv);
}

