//
// Created by Wiliam SAENZ AREVALO on 20/05/2025.
//
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TTree.h>
#include <TGraph.h>
#include <TFile.h>
#include <TArc.h>
#include "THATFitMethods.hxx"
#include <iostream>
#include <iomanip>  // for std::setw
#include <chrono>
#include <thread>   // for std::this_thread::sleep_for
#include <filesystem>


bool save_plots   = false; // prints example plots of track fit
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

double GetBarycenter_Z(double z0, double z1, double z2, double q0, double q1, double q2, bool h0, bool h1, bool h2) {
    double sum = 0.0, total = 0.0;
    if (h0) { sum += z0 * q0; total += q0; }
    if (h1) { sum += z1 * q1; total += q1; }
    if (h2) { sum += z2 * q2; total += q2; }
    return (total == 0.0) ? 0.0 : sum / total;
}

double GetBarycenter_Z_Error(double z0, double z1, double z2, double q0, double q1, double q2, bool h0, bool h1, bool h2) {
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

void FillMinMaxValues(TClusterFeatures &clus_features, double cenZ, double cenY, double &phi_min, double &phi_max) {
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

void ProcessCluster(TClusterFeatures &clu_feat,
                    double q0_amp, double q1_amp, double q2_amp,
                    double q0_int, double q1_int, double q2_int,
                    double y0, double z0, double z1, double z2,
                    bool h0, bool h1, bool h2, int clu_size,
                    TH1F *q0_amp_hist, TH1F *q1_amp_hist, TH1F *q2_amp_hist,
                    TH1F *q0_int_hist, TH1F *q1_int_hist, TH1F *q2_int_hist,
                    TH1F *q1_0_amp_hist, TH1F *q1_0_int_hist,
                    TH1F *q2_1_amp_hist, TH1F *q2_1_int_hist) {

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

    clu_feat.q.push_back(TMath::Log(q1_amp / q0_amp));
    clu_feat.ori.push_back(0);
    clu_feat.y.push_back(y0);
    if(use_integral) {
        clu_feat.z  .push_back(GetBarycenter_Z(z0, z1, z2, q0_int, q1_int, q2_int, h0, h1, h2));
        clu_feat.err.push_back(GetBarycenter_Z_Error(z0, z1, z2, q0_int, q1_int, q2_int, h0, h1, h2));
    }
    else{
        clu_feat.z  .push_back(GetBarycenter_Z(z0, z1, z2, q0_amp, q1_amp, q2_amp, h0, h1, h2));
        clu_feat.err.push_back(GetBarycenter_Z_Error(z0, z1, z2, q0_amp, q1_amp, q2_amp, h0, h1, h2));
    }

    clu_feat.size.push_back(clu_size);
}

void ProcessEvent(int event, int pattern, TClusterFeatures &clu_feat, TH1F * residuals) {

    std::vector<double> zp(clu_feat.z.size()), up(clu_feat.z.size());
    double param_circ[3], chi2_circ, zfit, yfit, phic, Rfit, cov_circ[6];

    ND::THATFitMethods::FitHelix((int)clu_feat.z.size(), clu_feat.z.data(), clu_feat.y.data(), clu_feat.ori.data(),
                                 clu_feat.err.data(), zp.data(), up.data(), zfit, yfit, phic, Rfit,
                                 param_circ, chi2_circ, cov_circ, true);

    auto *track = new TGraphErrors();
    for (size_t i = 0; i < clu_feat.z.size(); ++i) {
        track->SetPoint     (i, clu_feat.z[i], clu_feat.y[i]);
        track->SetPointError(i, clu_feat.err[i], 0);
    }

    double centerZ = zfit - Rfit * sin(phic);
    double centerY = yfit + Rfit * cos(phic);
    if (Rfit < 0) {
        Rfit = -Rfit;
        centerZ = zfit - Rfit * sin(phic + M_PI);
        centerY = yfit + Rfit * cos(phic + M_PI);
    }


    // Here we fill the residuals histogram with Delta Z =

    for (size_t i = 0; i < clu_feat.z.size(); ++i) {

        double res = pow(pow((clu_feat.z[i] - centerZ), 2.) + pow((clu_feat.y[i]  - centerY), 2.),0.5) - Rfit;
        residuals ->Fill(res);

    }


    if (save_plots && event < 200) {

        double phi_min = 1e9, phi_max = -1e9;
        FillMinMaxValues(clu_feat, centerZ, centerY, phi_min, phi_max);
        track->GetXaxis()->SetTitle("z (mm)");
        track->GetYaxis()->SetTitle("y (mm)");
        track->SetMarkerStyle(20); track->SetMarkerColor(kBlack); track->Draw("AP");
        auto *arc = new TArc(centerZ, centerY, std::abs(Rfit), phi_min * 180 / TMath::Pi(), phi_max * 180 / TMath::Pi());
        arc->SetLineWidth(2); arc->SetLineColor(kBlue); arc->SetFillStyle(0); arc->Draw("sameL");

        auto *canvas = new TCanvas("c", "", 800, 600);
        canvas->SetLeftMargin(0.15);  // Increase if needed (default is ~0.1)
        if (! std::filesystem::exists("plots")) {
            std::filesystem::create_directory("plots");
        }
        track->Draw("AP");
        arc->Draw("same");
        canvas->SaveAs(Form("plots/pattern_gr_%i_%i.pdf", event, pattern));
        delete canvas;

    }
}

void analyse_data(int argc, char** argv) {
    std::string inputFile = argv[1], outputFile = argv[2];
    TFile *myf = TFile::Open(inputFile.c_str());
    if (!myf) throw std::runtime_error("Input file not found.");

    TTree *tree = (TTree*) myf->Get("outTree");
    if (!tree) throw std::runtime_error("TTree not found in input file.");

    Int_t event, pattern, cluster, hit, clu_size, trk_ncl;
    Double_t hit_pos[3], trk_pos[3], trk_dir[3];
    TGraph *wf = nullptr;

    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("pattern", &pattern);
    tree->SetBranchAddress("cluster", &cluster);
    tree->SetBranchAddress("hit", &hit);
    tree->SetBranchAddress("hit_pos", &hit_pos);
    tree->SetBranchAddress("trk_pos", &trk_pos);
    tree->SetBranchAddress("trk_dir", &trk_dir);
    tree->SetBranchAddress("clu_size", &clu_size);
    tree->SetBranchAddress("WF", &wf);
    tree->SetBranchAddress("trk_ncl", &trk_ncl);

    auto *q0_amp = new TH1F("q0_amp", "q0_amp", 1000, 0, 2000);
    auto *q1_amp = new TH1F("q1_amp", "q1_amp", 1000, 0, 2000);
    auto *q2_amp = new TH1F("q2_amp", "q2_amp", 1000, 0, 2000);
    auto *q0_int = new TH1F("q0_int", "q0_int", 1000, 0, 50000);
    auto *q1_int = new TH1F("q1_int", "q1_int", 1000, 0, 50000);
    auto *q2_int = new TH1F("q2_int", "q2_int", 1000, 0, 50000);
    auto *q1_0_amp = new TH1F("q1_0_amp", "q1_0_amp", 200, 0, 2);
    auto *q2_1_amp = new TH1F("q2_1_amp", "q2_1_amp", 200, 0, 2);
    auto *q1_0_int = new TH1F("q1_0_int", "q1_0_int", 200, 0, 2);
    auto *q2_1_int = new TH1F("q2_1_int", "q2_1_int", 200, 0, 2);
    auto *residuals = new TH1F("residuals", "residuals", 200, -10, 10);

    double q0a, q1a, q2a, q0i, q1i, q2i, y0, y1, y2, z0, z1, z2;
    bool h0 = false, h1 = false, h2 = false;
    TClusterFeatures clu_feat;
    int prev_clu = -999, prev_evt = -999;
    double thr = 10.;

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; ++i) {

        tree->GetEntry(i);

        // Print progress
        float progress = static_cast<float>(i + 1) / nentries;
        printProgressBar(progress);

        if (pattern != 0 )  // only event with single pattern
            continue;

        if (!(clu_size > 1 && trk_ncl > 50 && fabs(trk_pos[0]) > 900 && fabs(trk_dir[1]) > 0.9))
            continue;

        if (prev_clu != cluster) {
            ProcessCluster(clu_feat, q0a, q1a, q2a, q0i, q1i, q2i, y0, z0, z1, z2, h0, h1, h2, clu_size, q0_amp, q1_amp,
                           q2_amp, q0_int, q1_int, q2_int, q1_0_amp, q1_0_int, q2_1_amp,
                           q2_1_int);
                           h0 = h1 = h2 = false;
                           prev_clu = cluster;
        }

        if (prev_evt != event) {
            ProcessEvent(event, pattern, clu_feat, residuals);
            clu_feat.clear();
            prev_evt = event;
        }

        double amp = -9999.9;
        for (int j = 0; j < wf->GetN(); ++j)
            amp = std::max(amp, wf->GetPointY(j));

        if (hit == 0) {
            h0 = true;
            q0a = amp;
            q0i = GetIntegral(thr, wf);
            y0 = hit_pos[1];
            z0 = hit_pos[2]; }

        else if (hit == 1) {
            h1 = true;
            q1a = amp;
            q1i = GetIntegral(thr, wf);
            y1 = hit_pos[1];
            z1 = hit_pos[2]; }

        else if (hit == 2) {
            h2 = true;
            q2a = amp;
            q2i = GetIntegral(thr, wf);
            y2 = hit_pos[1];
            z2 = hit_pos[2]; }
    }

    TCanvas *myc = new TCanvas(); myc->Divide(2,2);
    myc->cd(1); q1_amp->Draw(); q0_amp->Draw("same");
    myc->cd(2); q1_int->Draw(); q0_int->Draw("same");
    myc->cd(3); q1_0_amp->Draw(); q1_0_int->Draw("same");
    myc->cd(4); q2_1_amp->Draw(); q2_1_int->Draw("same");

    TFile *outfile = new TFile(outputFile.c_str(), "RECREATE");
    myc->Write(); q0_amp->Write(); q1_amp->Write(); q2_amp->Write();
    q0_int->Write(); q1_int->Write(); q2_int->Write();
    q1_0_amp->Write(); q1_0_int->Write(); q2_1_amp->Write(); q2_1_int->Write();
    residuals->Write();
    outfile->Close();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_root_file>\n";
        return 1;
    }
    analyse_data(argc, argv);
}
