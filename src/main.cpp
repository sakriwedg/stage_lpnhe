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
#include <fstream>
#include <TPaveText.h>
#include <cmath>

double R = 100;
double C = 1;
int bord = 40;
int N_hits = 100000;  // nombre max de hits à traiter
int ADC = 100; //Amplitude en dessous de laquelle on a du bruit

void plot_wf(int argc, char** argv);
int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_root_file>\n";
        return 1;
    }

    plot_wf(argc, argv);

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


    for (int i = 0; i < nentries && charges.size() < N_hits; ++i) {
        tree->GetEntry(i);
        //if (cluster_per_trace > bord){

            if (hit == hit_filter) {
                double charge = integrale(10, wf);
                double v_max = coordypic(wf);
                double t_at_max = coordxpic(wf);

                charges.push_back(charge);
                peak_times.push_back(t_at_max);
                max_amplitudes.push_back(v_max);
                //}
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void plot_wf(int argc, char** argv) {
    int count_hit0 = 0;
    int count_hit1 = 0;
    int count_hit2 = 0;
    int count_hit3 = 0;
    int test_hit = 0;

    //Vecteur test pour hit/cluster
    std::vector<int> vtest_hit;

    //Vecteur à remplir por histo des clusters par trace
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

    // Linking branches to variables
    tree->SetBranchAddress("event",    &event, &b_event);
    tree->SetBranchAddress("pattern",  &pattern, &b_pattern);
    tree->SetBranchAddress("cluster",  &cluster, &b_cluster);
    tree->SetBranchAddress("hit",      &hit, &b_hit);
    tree->SetBranchAddress("trk_pos",  &trk_pos, &b_trk_pos);
    tree->SetBranchAddress("clu_size", &clu_size, &b_clu_size);
    tree->SetBranchAddress("WF",       &wf);
    tree-> SetBranchAddress("trk_ncl",  &cluster_per_trace, &b_cluster_per_trace);


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
    std::vector<int> DOS;
    std::vector<int> TRES;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HISTO DES CLUSTERS PAR TRACE ET COUPURE

    // Reading TTree and filling data of interest
    for (int i = 0; i < nentries && i < N_hits ; ++i) {
        tree -> GetEntry(i);
        if (cluster_per_trace > bord){  //On regarde si on a bien un nombre de cluster assez grand qui signifie que l'on a bien un trace est pas du bruit


        int test_trk = 0;  //On initialise les test pour avoir les clusters par trace
        int test_trk1 = 0;


        test_trk = cluster_per_trace;
        tree -> GetEntry(i+1);
        test_trk1 = cluster_per_trace;

        if (test_trk != test_trk1){  //on regarde si on est dans la même trace
            cluster_par_trace.push_back(test_trk);
            }



        tree -> GetEntry(i);
        // Looking at a single track

        //if (event == look_event && pattern == look_pattern && cluster == look_cluster) {
            if (hit == 0){

            vtest_hit.push_back(0);

            double charge = integrale(10, wf);
            double val_lim = coordypic(wf);


             //std::cout << "t_at_max ="  <<coordxpic(wf) << std::endl;
             //std::cout << "max_value =" <<coordypic(wf) << std::endl;
             //std::cout << "charge ="    <<integrale(0.05 * val_lim , wf) << std::endl;

                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_1   -> SetPoint(j, t, q);
                }
            }

            if (hit == 1){
                vtest_hit.push_back(1);
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_2   -> SetPoint(j, t, q);
                }
            }
            if (hit == 2){
                vtest_hit.push_back(2);
                count_hit2++;
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_3   -> SetPoint(j, t, q);
                }
            }
            if (hit == 3){
                vtest_hit.push_back(3);
                for (int j = 0; j < wf->GetN(); ++j) {
                    double t, q;
                    wf     -> GetPoint(j, t, q);
                    wf_4   -> SetPoint(j, t, q);
                }
            }
        //}
        }
    }
    //Traitement des hits

    //Les listes des ZERO serviront à calculer le premier rartio les UNO à l'autre
    std::vector<int> ZERO0;
    std::vector<int> ZERO1;
    std::vector<int> UNO0;
    std::vector<int> UNO1;
    std::vector<int> temp;


    //On utilise ces mêmes vecteurs avec un B à la fin pour dire Bis : objectif de comparer charger et amplitude
    std::vector<int> ZERO0B;
    std::vector<int> ZERO1B;
    std::vector<int> UNO0B;
    std::vector<int> UNO1B;


    for (int h = 0; h < vtest_hit.size(); ++h){
        tree -> GetEntry(h);
        temp.push_back(vtest_hit[h]);



        if (h + 1 < vtest_hit.size() && vtest_hit[h] + 1 != vtest_hit[h+1]){
            hit_average.push_back(temp.size());
            temp.clear();
            }


        if (h + 1 < vtest_hit.size() && vtest_hit[h] + 1 == vtest_hit[h+1] && vtest_hit[h] == 0){
               ZERO0.push_back(integrale(10, wf));
               ZERO0B.push_back(coordypic(wf));

               tree -> GetEntry(h+1);

               ZERO1.push_back(integrale(10, wf));
               ZERO1B.push_back(coordypic(wf));
            }


        if (h + 1 < vtest_hit.size() && vtest_hit[h] + 1 == vtest_hit[h+1] && vtest_hit[h] == 1){
                UNO0.push_back(integrale(10, wf));
                UNO0B.push_back(coordypic(wf));

                tree -> GetEntry(h+1);

                UNO1.push_back(integrale(10, wf));
                UNO1B.push_back(coordypic(wf));
            }
    }


    // Plotting

    //Pour les ratios de Q1/Q0 et Q2/Q1
    std::vector<double> ratio_ZERO;
    std::vector<double> ratio_UNO;

    for (size_t i = 0; i < ZERO0.size(); ++i) {
            ratio_ZERO.push_back((double)ZERO1[i] / ZERO0[i]);
        }

    for (size_t i = 0; i < UNO0.size(); ++i) {
            ratio_UNO.push_back((double)UNO1[i] / UNO0[i]);
        }


    std::vector<double> ratio_ZEROB;
    std::vector<double> ratio_UNOB;

    for (size_t i = 0; i < ZERO0B.size(); ++i) {
            ratio_ZEROB.push_back((double)ZERO1B[i] / ZERO0B[i]);
        }

    for (size_t i = 0; i < UNO0B.size(); ++i) {
            ratio_UNOB.push_back((double)UNO1B[i] / UNO0B[i]);
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

    TLegend * myl = new TLegend();


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
        TH1F *h_charge    = new TH1F(Form("h_charge_hit%d",    hit_filter), Form("Distribution des charges (hit=%d);Charge;Nombre de hits", hit_filter), 50, 0, 15000);
        TH1F *h_peak_time = new TH1F(Form("h_peak_time_hit%d", hit_filter), Form("Temps du pic max (hit=%d);Temps (ns);Nombre de hits",     hit_filter), 50, 0, 511);
        TH1F *h_amplitude = new TH1F(Form("h_amplitude_hit%d", hit_filter), Form("Amplitude max (hit=%d);Amplitude (ADC);Nombre de hits",   hit_filter), 100, 0, 1000);

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

    h0 -> Draw();
    h1 -> Draw("same");
    h2 -> Draw("same");
    h3 -> Draw("same");


    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h0, "Hit 0", "l");
    leg->AddEntry(h1, "Hit 1", "l");
    leg->AddEntry(h2, "Hit 2", "l");
    leg->AddEntry(h3, "Hit 3", "l");
    leg->Draw();


    //On veut afficher les entries sur ce même histo
    TPaveText *infos_entries = new TPaveText(0.15, 0.75, 0.4, 0.9, "NDC");

    infos_entries -> AddText(Form("Entries hit 0 : %d", (int)h0 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 1 : %d", (int)h1 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 2 : %d", (int)h2 -> GetEntries()));
    infos_entries -> AddText(Form("Entries hit 3 : %d", (int)h3 -> GetEntries()));


    infos_entries -> Draw();
    c_all_charges -> Modified();
    c_all_charges -> Update();
    c_all_charges -> Write();

    outfile -> cd();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Histogramme QO_A sachant qu'il existe Q1
    TH1F *h_Q0_s_Q1 = new TH1F("Q0_sachant_Q1_A", "Q0_sachant_Q1_A", 100, 0, 1400);
    for (double r : ZERO0B) h_Q0_s_Q1 -> Fill(r);
    h_Q0_s_Q1 -> Write();

    //Histogramme de Q1
    TH1F *h_Q1 = new TH1F("Q1_A", "Q1_A", 100, 0, 1400);
    for (double r : UNO0B) h_Q1 -> Fill(r);
    h_Q1 -> Write();

    //Histogramme Q1_A sachant qu'il existe Q2
    TH1F *h_Q1_s_Q2= new TH1F("Q1_sachant_Q2_A", "Q1_sachant_Q2_A", 100, 0, 1400);
    for (double r : ZERO1B) h_Q1_s_Q2 -> Fill(r);
    h_Q1_s_Q2-> Write();

    //Histogramme de Q2
    TH1F *h_Q2 = new TH1F("Q2_A", "Q2_A", 100, 0, 1400);
    for (double r : UNO1B) h_Q2 -> Fill(r);
    h_Q2 -> Write();




    //Histo du nombre de hits par cluster
    TH1F *hist = new TH1F("hit_per_cluster", "Histo des hit/cluster;Hit N;Nombre de cluster avec N hits", 6, 0, 6);
    for (int val : hit_average) {
        hist -> Fill(val);
    }
    hist->Write();



    //Histo du nombre de cluster par trace
    TH1F *h_cluster_par_trace = new TH1F("cluster_per_trace", "Clusters per trace", 100, 0, 100);
    for (double r : cluster_par_trace) h_cluster_par_trace -> Fill(r);
    h_cluster_par_trace -> Write();



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Histo du ratio des charges
    TH1F *h_ratio_ZERO = new TH1F("Q1_sur_Q0_charge", "Ratio Q1/Q0;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_ZERO) h_ratio_ZERO -> Fill(r);
    h_ratio_ZERO -> Write();

    TH1F *h_ratio_UNO = new TH1F("Q2_sur_Q1_charge", "Ratio Q2/Q1;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_UNO) h_ratio_UNO -> Fill(r);
    h_ratio_UNO -> Write();


    //Histo du ln du ratio des charges
    TH1F *h_lnratio_ZERO = new TH1F("log_Q1_sur_Q0_charge", "lnQ1/Q0;Ratio;Nombre d'entrees", 100, -5, 5);
    for (double r : ratio_ZERO) h_lnratio_ZERO -> Fill(std::log(r));
    h_lnratio_ZERO -> Write();

    TH1F *h_lnratio_UNO = new TH1F("log_Q2_sur_Q1_charge", "lnQ2/Q1 Charge;Ratio;Nombre d'entrees", 100, -5, 5);
    for (double r : ratio_UNO) h_lnratio_UNO -> Fill(std::log(r));
    h_lnratio_UNO -> Write();


    //Histo du ratio des amplitudes
    TH1F *h_ratio_ZEROB = new TH1F("Q1_sur_Q0_amplitude", "Ratio Q1/Q0;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_ZEROB) h_ratio_ZEROB -> Fill(r);
    h_ratio_ZEROB -> Write();

    TH1F *h_ratio_UNOB = new TH1F("Q2_sur_Q1_amplitude", "Ratio Q2/Q1;Ratio;Nombre d'entrees", 50, 0, 2);
    for (double r : ratio_UNOB) h_ratio_UNOB -> Fill(r);
    h_ratio_UNOB -> Write();


    //Histo du ln du ratio des amplitudes
    TH1F *h_lnratio_ZEROB = new TH1F("log_Q1_sur_Q0_amplitude", "lnQ1/Q0;Ratio;Nombre d'entrees", 100, -5, 5);
    for (double r : ratio_ZEROB) h_lnratio_ZEROB -> Fill(std::exp(log(r)));
    h_lnratio_ZEROB -> Write();

    TH1F *h_lnratio_UNOB = new TH1F("log_Q2_sur_Q1_amplitude", "lnQ2/Q1 ;Ratio;Nombre d'entrees", 100, -5, 5);
    for (double r : ratio_UNOB) h_lnratio_UNOB -> Fill(std::log(r));
    h_lnratio_UNOB -> Write();





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
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1 -> AddEntry(h_ratio_ZEROB, " Q1 / Q0 Amplitude", "l");
    leg1 -> AddEntry(h_ratio_ZERO, "Q1 / Q0 Charge", "l");
    leg1 -> Draw();
    c_ratio_01 -> Write();




    TCanvas *c_ratio_02 = new TCanvas("c_ratio_02");
    h_ratio_UNOB -> SetLineColor(1);
    h_ratio_UNO -> SetLineColor(2);

    double max3 = h_ratio_UNOB -> GetMaximum();
    double max4 = h_ratio_UNO -> GetMaximum();
    double ymax2 = std::max( max3, max4) * 1.1;
    h_ratio_UNOB -> SetMaximum(ymax);
    h_ratio_UNO -> SetMaximum(ymax);

    h_ratio_UNOB -> Draw();
    h_ratio_UNO -> Draw("SAME");
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2 -> AddEntry(h_ratio_UNOB, " Q2 / Q1 Amplitude", "l");
    leg2 -> AddEntry(h_ratio_UNO, " Q2 / Q1 Charge", "l");
    leg2 -> Draw();
    c_ratio_02 -> Write();




    TCanvas *c_lnratio_01 = new TCanvas("c_lnratio_01");
    h_lnratio_ZEROB -> SetLineColor(1);
    h_lnratio_ZERO->SetLineColor(2);

    double max5 = h_lnratio_ZEROB -> GetMaximum();
    double max6 = h_lnratio_ZERO -> GetMaximum();
    double ymax3 = std::max( max5, max6) * 1.1;
    h_lnratio_ZEROB  -> SetMaximum(ymax);
    h_lnratio_ZERO -> SetMaximum(ymax);

    h_lnratio_ZEROB -> Draw();
    h_lnratio_ZERO -> Draw("SAME");
    TLegend *legln1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legln1 -> AddEntry(h_lnratio_ZEROB, "ln (Q1 / Q0) Amplitude", "l");
    legln1 -> AddEntry(h_lnratio_ZERO, "ln (Q1 / Q0) Charge", "l");
    legln1 -> Draw();
    c_lnratio_01 -> Write();




    TCanvas *c_lnratio_02 = new TCanvas("c_lnratio_02");
    h_lnratio_UNOB -> SetLineColor(1);
    h_lnratio_UNO->SetLineColor(2);

    double max8 = h_lnratio_UNOB -> GetMaximum();
    double max7 = h_lnratio_UNO -> GetMaximum();
    double ymax4 = std::max( max7, max8) * 1.1;
    h_lnratio_UNOB  -> SetMaximum(ymax);
    h_lnratio_UNO -> SetMaximum(ymax);

    h_lnratio_UNOB -> Draw();
    h_lnratio_UNO -> Draw("SAME");
    TLegend *legln2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legln2 -> AddEntry(h_lnratio_UNOB, "ln (Q2 / Q1) Amplitude", "l");
    legln2 -> AddEntry(h_lnratio_UNO, "ln (Q2 / Q1) Charge", "l");
    legln2 -> Draw();
    c_lnratio_02 -> Write();


// CURVE FIT A FAIRE


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

//FAIRE 2 HISTOS/RATION q0/q1   //q1/q2                 OK

//ON Va calculer q avec amplitude / intégrale  et faire ces 2 histos,             OK

//puis on calculeera l'écart type





//lancer le code  avec des vraies données

//faire un fit sur les derniers histos obtenus
//comment utiliser git
// faire git push()
