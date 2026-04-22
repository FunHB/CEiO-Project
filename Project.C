void Project() {
    TFile *file = new TFile("V0_data.root");
    TTree *tree = (TTree*)f1_s->Get("Ks2pipi_dd/DecayTree");

    Double_t h1_PX;
    Double_t h1_PY;
    Double_t h1_PZ;
    Double_t h1_M;

    Double_t h2_PX;
    Double_t h2_PY;
    Double_t h2_PZ;
    Double_t h2_M;

    tree->SetBranchAddress("h1_PX", &h1_PX);
    tree->SetBranchAddress("h1_PY", &h1_PY);
    tree->SetBranchAddress("h1_PZ", &h1_PZ);
    tree->SetBranchAddress("h1_M",  &h1_M);

    tree->SetBranchAddress("h2_PX", &h2_PX);
    tree->SetBranchAddress("h2_PY", &h2_PY);
    tree->SetBranchAddress("h2_PZ", &h2_PZ);
    tree->SetBranchAddress("h2_M",  &h2_M);

    int nev = tree->GetEntries();

    for (int i = 0; i < nev; ++i) {
        Double_t h1_E = sqrt(pow(h1_PX, 2) + pow(h1_PY, 2) + pow(h1_PZ, 2) + pow(h1_M, 2));
        TLorentzVector h1_P(h1_PX, h1_PY, h1_PZ, h1_E);

        Double_t h2_E = sqrt(pow(h2_PX, 2) + pow(h2_PY, 2) + pow(h2_PZ, 2) + pow(h2_M, 2));
        TLorentzVector h2_P(h2_PX, h2_PY, h2_PZ, h2_E);
    }
}