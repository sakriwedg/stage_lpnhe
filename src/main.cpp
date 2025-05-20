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

void plot_wf(int argc, char** argv);

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_data_file> <output_root_file>\n";
        return 1;
    }

    plot_wf(argc, argv);

}

void plot_wf(int argc, char** argv) {

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

    // Linking branches to variables
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("pattern", &pattern, &b_pattern);
    tree->SetBranchAddress("cluster", &cluster, &b_cluster);
    tree->SetBranchAddress("hit", &hit, &b_hit);
    tree->SetBranchAddress("trk_pos", &trk_pos, &b_trk_pos);
    tree->SetBranchAddress("clu_size", &clu_size, &b_clu_size);
    tree->SetBranchAddress("WF", &wf);


    Int_t nentries = tree->GetEntries();
    std::cout << "entries = " << nentries << std::endl;


    TGraph * wf_1 = new TGraph();
    TGraph * wf_2 = new TGraph();
    TGraph * wf_3 = new TGraph();
    TGraph * wf_4 = new TGraph();

    Int_t look_event   = 0;
    Int_t look_pattern = 0;
    Int_t look_cluster = 7;

    // Reading TTree and filling data of interest
    for (int i = 0; i < nentries; ++i) {

        tree->GetEntry(i);

        // Looking at a single track

        if (event == look_event && pattern == look_pattern && cluster == look_cluster) {

            if (hit == 0){

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
        }
    }


    // Plotting
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




    // Data export
    TFile * outfile = new TFile(outputFile.c_str(),"recreate");
    myc->Write();
    outfile->Close();

}
