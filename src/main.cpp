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

struct   TClusterFeatures{
    std::vector<double>  time, phi, x, z, y, dir_x, dir_y, dir_z, err, q;
    std::vector<double>  y_eram, z_eram;
    std::vector<int>     ori, size;

    void clear() {
        time.clear(); phi.clear(); x.clear(); z.clear(); y.clear();
        dir_x.clear(); dir_y.clear(); dir_z.clear(); err.clear(); q.clear();
        y_eram.clear(); z_eram.clear();
        ori.clear(); size.clear();
    }

};

double   GetBarycenter_Z(double z_h0, double z_h1, double z_h2, double q0_amp, double q1_amp, double q2_amp, bool hit_0, bool hit_1, bool hit_2) {
    double weighted_sum = 0.0;
    double total_charge = 0.0;

    if (hit_0) {
        weighted_sum += z_h0 * q0_amp;
        total_charge += q0_amp;
    }

    if (hit_1) {
        weighted_sum += z_h1 * q1_amp;
        total_charge += q1_amp;
    }

    if (hit_2) {
        weighted_sum += z_h2 * q2_amp;
        total_charge += q2_amp;
    }

    // Avoid division by zero
    if (total_charge == 0.0) {
        return 0.0; // or return NAN, or throw, depending on use case
    }

    return weighted_sum / total_charge;
}

double   GetBarycenter_Z_Error(double z_h0, double z_h1, double z_h2, double q0_amp, double q1_amp, double q2_amp, bool hit_0, bool hit_1, bool hit_2) {
    double weighted_sum = 0.0;
    double total_weight = 0.0;

    if (hit_0) {
        weighted_sum += z_h0 * q0_amp;
        total_weight += q0_amp;
    }
    if (hit_1) {
        weighted_sum += z_h1 * q1_amp;
        total_weight += q1_amp;
    }
    if (hit_2) {
        weighted_sum += z_h2 * q2_amp;
        total_weight += q2_amp;
    }

    if (total_weight == 0.0) return 0.0;

    double z_bar = weighted_sum / total_weight;

    // Now compute variance
    double variance_sum = 0.0;
    if (hit_0) {
        double dz = z_h0 - z_bar;
        variance_sum += q0_amp * dz * dz;
    }
    if (hit_1) {
        double dz = z_h1 - z_bar;
        variance_sum += q1_amp * dz * dz;
    }
    if (hit_2) {
        double dz = z_h2 - z_bar;
        variance_sum += q2_amp * dz * dz;
    }

    return std::sqrt(variance_sum / total_weight);
}

void     analyse_data(int argc, char** argv);

Double_t GetIntegral(Double_t thr, TGraph * my_wf){
    double integral = 0.;
    for (int j = 0; j < my_wf->GetN(); ++j) {
        auto amp_j = my_wf->GetPointY(j);
        if ( amp_j > thr)
            integral += amp_j;
    }
    return integral;
}

void     FillMinMaxValues(TClusterFeatures & clus_features, double cen_cir_Z, double cen_cir_Y, double & phi_min, double & phi_max){

    for (unsigned int i = 0; i < clus_features.z.size(); i++) {

        double delta_z = clus_features.z[i] - cen_cir_Z;
        double delta_y = clus_features.y[i] - cen_cir_Y;
        double phi_i   = atan2(delta_y, delta_z);

        clus_features.phi.push_back(phi_i);
        if (clus_features.size[i] > 0) { // Only consider clusters with size > 1 for defining track and node states
            if (phi_min > atan2(delta_y, delta_z))
                phi_min = atan2(delta_y, delta_z);

            if (phi_max < atan2(delta_y, delta_z))
                phi_max = atan2(delta_y, delta_z);

        }
    }
}


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_root_file>\n";
        return 1;
    }

    analyse_data(argc, argv);

}

void analyse_data(int argc, char** argv) {

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
    Int_t trk_ncl;

    Double_t hit_pos[3];
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
    TBranch *b_hit_pos;   //!
    TBranch *b_trk_dir;   //!
    TBranch *b_clu_size;   //!
    TBranch *b_wf;   //!
    TBranch *b_trk_ncl;   //!

    // Linking branches to variables
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("pattern", &pattern, &b_pattern);
    tree->SetBranchAddress("cluster", &cluster, &b_cluster);
    tree->SetBranchAddress("hit", &hit, &b_hit);
    tree->SetBranchAddress("hit_pos", &hit_pos, &b_hit_pos);
    tree->SetBranchAddress("trk_pos", &trk_pos, &b_trk_pos);
    tree->SetBranchAddress("trk_dir", &trk_dir, &b_trk_dir);
    tree->SetBranchAddress("clu_size", &clu_size, &b_clu_size);
    tree->SetBranchAddress("WF", &wf);
    tree->SetBranchAddress("trk_ncl", &trk_ncl, &b_trk_ncl);


    Int_t nentries = tree->GetEntries();
    std::cout << "N entries = " << nentries << std::endl;

    // histogram declaration
    TH1F * q0_amp = new TH1F("q0_amp","q0_amp",1000,0,2000);
    TH1F * q1_amp = new TH1F("q1_amp","q1_amp",1000,0,2000);
    TH1F * q2_amp = new TH1F("q2_amp","q2_amp",1000,0,2000);

    TH1F * q1_0_amp = new TH1F("q1_0_amp","q1_0_amp",200,0,2);
    TH1F * q2_1_amp = new TH1F("q2_1_amp","q2_1_amp",200,0,2);

    TH1F * q0_int = new TH1F("q0_int","q0_int",1000,0,50000);
    TH1F * q1_int = new TH1F("q1_int","q1_int",1000,0,50000);
    TH1F * q2_int = new TH1F("q2_int","q2_int",1000,0,50000);

    TH1F * q1_0_int = new TH1F("q1_0_int","q1_0_int",200,0,2);
    TH1F * q2_1_int = new TH1F("q2_1_int","q2_1_int",200,0,2);

    double q_h0_amp, q_h1_amp, q_h2_amp;
    double q_h0_int, q_h1_int, q_h2_int;
    double y_h0, y_h1, y_h2;
    double z_h0, z_h1, z_h2;

    bool hit_0 = false, hit_1 = false, hit_2 = false;

    Int_t previous_clu = - 9;
    Int_t previous_eve = - 9;
    Double_t my_threshold = 10.;

    TClusterFeatures clu_feat;

    // Reading TTree and filling data of interest
    for (int i = 0; i < nentries; ++i) {

        tree->GetEntry(i);

        // Just events with single pattern
        if (pattern != 0)
            continue;

        if (event > 200)
            break;

    	if( clu_size > 1 && trk_ncl > 50 && fabs(trk_pos[0]) > 900 &&  fabs(trk_dir[1]) > 0.9 ){
            //std::cout << " - passed filter" << std::endl;

            // This if determines when going from one cluster to the next
            if (previous_clu != cluster ) {

                //std::cout << " - new cluster" << std::endl;
                if (hit_0) {
                    q0_amp->Fill(q_h0_amp);
                    q0_int->Fill(q_h0_int);
                }
                if (hit_1) {
                    q1_amp->Fill(q_h1_amp);
                    q1_int->Fill(q_h1_int);
                }
                if (hit_2) {
                    q2_amp->Fill(q_h2_amp);
                    q2_int->Fill(q_h2_int);
                }
                if (hit_0 && hit_1) {
                    q1_0_amp->Fill(q_h1_amp / q_h0_amp);
                    q1_0_int->Fill(q_h1_int / q_h0_int);
                }
                if (hit_2 && hit_1) {
                    q2_1_amp->Fill(q_h2_amp / q_h1_amp);
                    q2_1_int->Fill(q_h2_int / q_h1_int);
                }

                // Filling cluster position q, y, z and orientation (0)... z computed from barycenter, q from Log(Q1/Q0)
                clu_feat.q   .push_back(TMath::Log(q_h1_amp/q_h0_amp));
                clu_feat.ori .push_back(0);
                clu_feat.y   .push_back(y_h0);
                clu_feat.z   .push_back(GetBarycenter_Z      (z_h0,z_h1,z_h2,q_h0_amp,q_h1_amp,q_h2_amp,hit_0,hit_1,hit_2));
                clu_feat.err .push_back(GetBarycenter_Z_Error(z_h0,z_h1,z_h2,q_h0_amp,q_h1_amp,q_h2_amp,hit_0,hit_1,hit_2));
                clu_feat.size.push_back(clu_size);


                // Resetting pad booleans
                previous_clu = cluster;
                hit_0 = false;
                hit_1 = false;
                hit_2 = false;

            }

            // This if determines when going from one event to the next
            if (previous_eve != event) {
                std::cout << "new event" << std::endl;


                std::vector<double> zp, up;
                zp.resize(clu_feat.z.size());
                up.resize(clu_feat.z.size());
                double param_circ[3], chi2_circ;
                double zfit, yfit, phic, Rfit;
                double cov_circ[6];

                int fit_result = ND::THATFitMethods::FitHelix(int(clu_feat.z.size()), clu_feat.z.data(), clu_feat.y.data(),
                                                              clu_feat.ori.data(), clu_feat.err.data(), zp.data(), up.data(), zfit,
                                                              yfit, phic, Rfit, param_circ, chi2_circ, cov_circ, true);


                auto * can = new TCanvas("c", "", 800, 600);
                can->cd();

                auto * track_ZY = new TGraphErrors();
                for (int j = 0; j < clu_feat.size.size(); ++j) {
                    track_ZY->SetPoint(track_ZY->GetN(), clu_feat.z[j], clu_feat.y[j]);
                    track_ZY->SetPointError(j,clu_feat.err[j],0);
                }

                track_ZY->SetMarkerStyle(20);
                track_ZY->SetMarkerSize(0.8);
                track_ZY->SetMarkerColor(kBlack);
                track_ZY->SetLineColor(kBlack);
                track_ZY->Draw("AP");


                double center_circle_Z, center_circle_Y;
                if (Rfit < 0) {
                    Rfit = -Rfit;
                    center_circle_Z = zfit - Rfit * sin(phic + M_PI);
                    center_circle_Y = yfit + Rfit * cos(phic + M_PI);
                } else {
                    center_circle_Z = zfit - Rfit * sin(phic);
                    center_circle_Y = yfit + Rfit * cos(phic);
                }

                double phi_min, phi_max;
                FillMinMaxValues(clu_feat, center_circle_Z, center_circle_Y, phi_min, phi_max);
                TArc *arc = new TArc(center_circle_Z, center_circle_Y, TMath::Abs(Rfit), phi_min * 180 / TMath::Pi(),
                                     phi_max * 180 / TMath::Pi());
                arc->SetLineWidth(2);
                arc->SetFillStyle(0);
                arc->SetLineColor(kBlue);
                arc->Draw("sameL");
                can->SaveAs(Form("plots/pattern_gr_%i_%i.pdf", event, pattern));

                // Add Delta Z vs Log Q point to 2D Histo

                previous_eve = event;
                track_ZY -> Clear();
                clu_feat.clear();
            }

            double amp = -9999.9;
            for (int j = 0; j < wf->GetN(); ++j) {
                auto amp_j = wf->GetPointY(j);
                if ( amp_j > amp)
                    amp = amp_j;
            }

            if      (hit == 0) {
                hit_0    = true;
                q_h0_amp = amp;
                q_h0_int = GetIntegral(my_threshold,wf);
                y_h0     = hit_pos[1];
                z_h0     = hit_pos[2];
            }
            else if (hit == 1){
                hit_1    = true;
                q_h1_amp = amp;
                q_h1_int = GetIntegral(my_threshold,wf);
                y_h1     = hit_pos[1];
                z_h1     = hit_pos[2];
            }
            else if (hit == 2){
                hit_2    = true;
                q_h2_amp = amp;
                q_h2_int = GetIntegral(my_threshold,wf);
                y_h2     = hit_pos[1];
                z_h2     = hit_pos[2];
            }


	    }
    }







    // Plotting and storing information

    q0_amp->SetLineColor(1);
    q1_amp->SetLineColor(2);
    q2_amp->SetLineColor(1);
    q1_0_amp->SetLineColor(1);
    q2_1_amp->SetLineColor(1);


    q2_int->SetLineColor(3);
    q0_int->SetLineColor(4);
    q1_int->SetLineColor(5);
    q1_0_int->SetLineColor(4);
    q2_1_int->SetLineColor(4);

    TCanvas * myc = new TCanvas();

    myc->Divide(2,2);


    myc->cd(1);
    q1_amp->Draw();
    q0_amp->Draw("same");

    myc->cd(2);
    q1_int->Draw();
    q0_int->Draw("same");

    myc->cd(3);
    q1_0_amp->Draw();
    q1_0_int->Draw("same");

    myc->cd(4);
    q2_1_amp->Draw();
    q2_1_int->Draw("same");

    // Data export
    TFile * outfile = new TFile(outputFile.c_str(),"recreate");

    myc ->Write();
    q0_amp->Write();
    q0_int->Write();
    q1_amp->Write();
    q1_int->Write();
    q2_amp->Write();
    q2_int->Write();

    q1_0_amp->Write();
    q2_1_amp->Write();
    q1_0_int->Write();
    q2_1_int->Write();
    outfile->Close();

}
