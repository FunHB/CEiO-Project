TH1F* CreateHist(const char* id, const char* labels, int bins, float min, float max) {
    TH1F *hist = new TH1F(id, labels, bins, min, max);
    hist->SetMarkerColor(9);
    hist->SetMarkerSize(1.1);
    hist->SetMarkerStyle(20);
    return hist;
}

void Project() {
    TFile *file = new TFile("V0_data.root");
    TTree *tree = (TTree*)file->Get("Ks2pipi_dd/DecayTree");

    Double_t h1_PX;
    Double_t h1_PY;
    Double_t h1_PZ;
    Double_t h1_M;

    Double_t h2_PX;
    Double_t h2_PY;
    Double_t h2_PZ;
    Double_t h2_M;

    Double_t KS0_OWNPV_X;
    Double_t KS0_OWNPV_Y;
    Double_t KS0_OWNPV_Z;

    Double_t KS0_ENDVERTEX_X;
    Double_t KS0_ENDVERTEX_Y;
    Double_t KS0_ENDVERTEX_Z;

    tree->SetBranchAddress("h1_PX", &h1_PX);
    tree->SetBranchAddress("h1_PY", &h1_PY);
    tree->SetBranchAddress("h1_PZ", &h1_PZ);
    tree->SetBranchAddress("h1_M",  &h1_M);

    tree->SetBranchAddress("h2_PX", &h2_PX);
    tree->SetBranchAddress("h2_PY", &h2_PY);
    tree->SetBranchAddress("h2_PZ", &h2_PZ);
    tree->SetBranchAddress("h2_M",  &h2_M);

    tree->SetBranchAddress("KS0_OWNPV_X",  &KS0_OWNPV_X);
    tree->SetBranchAddress("KS0_OWNPV_Y",  &KS0_OWNPV_Y);
    tree->SetBranchAddress("KS0_OWNPV_Z",  &KS0_OWNPV_Z);

    tree->SetBranchAddress("KS0_ENDVERTEX_X",  &KS0_ENDVERTEX_X);
    tree->SetBranchAddress("KS0_ENDVERTEX_Y",  &KS0_ENDVERTEX_Y);
    tree->SetBranchAddress("KS0_ENDVERTEX_Z",  &KS0_ENDVERTEX_Z);

    int nev = tree->GetEntries();

    auto hist_h1_M = CreateHist("hist_h1_M", " h1; Mass [MeV]", 40, 100, 200);
    auto hist_h2_M = CreateHist("hist_h2_M", " h2; Mass [MeV]", 40, 100, 200);
    auto hist_sum  = CreateHist("hist_sum",  " Sum of h1 i h2; Mass [MeV]", 100, 450, 550);
    auto hist_life  = CreateHist("hist_life",  " Life Time; Time [s]", 50, 0, 0.6e-9);

    for (int i = 0; i < nev; ++i) {
        tree->GetEntry(i);

        Double_t h1_E = sqrt(pow(h1_PX, 2) + pow(h1_PY, 2) + pow(h1_PZ, 2) + pow(h1_M, 2));
        TLorentzVector h1_P(h1_PX, h1_PY, h1_PZ, h1_E);

        Double_t h2_E = sqrt(pow(h2_PX, 2) + pow(h2_PY, 2) + pow(h2_PZ, 2) + pow(h2_M, 2));
        TLorentzVector h2_P(h2_PX, h2_PY, h2_PZ, h2_E);

        TLorentzVector sum_P = h1_P + h2_P;

        hist_h1_M->Fill(h1_P.M());
        hist_h2_M->Fill(h2_P.M());
        hist_sum->Fill(sum_P.M());

        ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> Own(1e-3 * KS0_OWNPV_X, 1e-3 * KS0_OWNPV_Y, 1e-3 * KS0_OWNPV_Z);
        ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> End(1e-2 * KS0_ENDVERTEX_X, 1e-2 * KS0_ENDVERTEX_Y, 1e-2 * KS0_ENDVERTEX_Z);

        auto diff = End - Own;

        if (abs(sum_P.M() - 500) < 10) {
            Double_t lifeTime = diff.R() * sum_P.M() / (sum_P.P() * 3.0e8);
            hist_life->Fill(lifeTime);
        }
    }

    hist_h1_M->Sumw2();
    hist_h2_M->Sumw2();
    hist_sum->Sumw2();
    hist_life->Sumw2();

    TCanvas *c1 = new TCanvas("c1", "", 100, 10, 640, 480);
    c1->SetLeftMargin(0.13);


    hist_h1_M->Draw();
    c1->Print("hist/hist_h1_M.png");


    hist_h2_M->Draw();
    c1->Print("hist/hist_h2_M.png");


    hist_sum->Draw();

    auto gauss = new TF1("hist_sum_gauss_fit", "gaus(0)", 490, 510);
    hist_sum->Fit("hist_sum_gauss_fit", "R");
    // gauss->Draw("same");

    auto gauss_pol2 = new TF1("hist_sum_fit", "gaus(0) + pol2(3)", 450, 550);
    gauss_pol2->SetParameters(0, gauss->GetParameter(1), gauss->GetParameter(2), 0, 0, 0);
    hist_sum->Fit("hist_sum_fit", "R");
    gauss_pol2->Draw("same");

    auto filtered_gauss = new TF1("filtered_gauss", "gaus(0)", 490, 508);
    filtered_gauss->SetParameters(gauss_pol2->GetParameter(0), gauss_pol2->GetParameter(1), gauss_pol2->GetParameter(2));

    auto bg = new TF1("bg", "pol2(0)", 490, 508);
    bg->SetParameters(gauss_pol2->GetParameter(3), gauss_pol2->GetParameter(4), gauss_pol2->GetParameter(5));

    double S = filtered_gauss->Integral(490, 508) / hist_sum->GetBinWidth(0);
    double B = bg->Integral(490, 508) / hist_sum->GetBinWidth(0);

    double u_S = S * gauss_pol2->GetParError(0) / gauss_pol2->GetParameter(0);
    double u_B = B * gauss_pol2->GetParError(3) / abs(gauss_pol2->GetParameter(3));

    cout << "\nS: " << S << "\n";
    cout << "B: " << B << "\n";
    cout << "S/B: " << S / B << "\n";
    cout << "u(S): " << u_S << "\n";
    cout << "u(B): " << u_B << "\n";

    c1->Print("hist/hist_sum.png");

    cout << "Particle Mass: " << hist_sum->GetMean() << " +- " << hist_sum->GetMeanError() << " \n";

    hist_life->Draw();

    auto hist_life_fit = new TF1("hist_life_fit", "expo", 0, 0.6e-9);
    hist_life->Fit("hist_life_fit", "R");
    hist_life_fit->Draw("same");

    c1->Print("hist/hist_life.png");
}