{

    TFile * f1 = TFile::Open("/Users/wsaenzar/post/nd/stage_lpnhe/build/int.root");

    TH1F * inte = (TH1F*) f1 -> Get("residuals");

    TFile * f2 = TFile::Open("/Users/wsaenzar/post/nd/stage_lpnhe/build/amp.root");

    TH1F * ampl = (TH1F*) f2 -> Get("residuals");

    inte->SetLineColor(1);
    ampl->SetLineColor(2);

    inte->Scale(1./inte->GetEntries());
    ampl->Scale(1./ampl->GetEntries());

    // Need to implement a fit for each distribution

    inte -> Draw ();
    ampl -> Draw ("same");

}