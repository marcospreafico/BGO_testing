
gErrorIgnoreLevel = kWarning;

  const double deltat_evt = 500e3; // delay between trg and plastic
   const double deltat_coinc = 30e3;

 struct hit{
    ULong64_t t; 
    UShort_t ch; 
    double A; 
    double Q; 
  };


  struct timeslice{
    vector<hit> hit;  
  };

struct event{
    ULong64_t tseed;
    UShort_t chseed; 
    double Aseed; 
    vector<hit> hit;  
  };

int find_evt(ULong64_t t, const vector<event> &evts){
            int id = -1; 
            for(int ii = 0; ii < evts.size(); ii++){
                ULong64_t deltat = (t > evts[ii].tseed) ? t-evts[ii].tseed : evts[ii].tseed-t;
                if(deltat<deltat_evt) id = ii;
            }
            return id;
        }

double langau(double *x, double *par) {

    //Fit parameters:
       //par[0]=Width (scale) parameter of Landau density
       //par[1]=Most Probable (MP, location) parameter of Landau density
       //par[2]=Total area (integral -inf to inf, normalization constant)
       //par[3]=Width (sigma) of convoluted Gaussian function
      
    // Numeric constants
          Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
          Double_t mpshift  = -0.22278298;       // Landau maximum location

          // Control constants
          Double_t np = 1000.0;      // number of convolution steps
          Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

          // Variables
          Double_t xx;
          Double_t mpc;
          Double_t fland;
          Double_t sum = 0.0;
          Double_t xlow,xupp;
          Double_t step;
          Double_t i;

          // MP shift correction
          mpc = par[1] - mpshift * par[0];

          // Range of convolution integral
          xlow = x[0] - sc * par[3];
          xupp = x[0] + sc * par[3];

          step = (xupp-xlow) / np;

          // Convolution integral of Landau and Gaussian by sum
          for(i=1.0; i<=np/2; i++) {
             xx = xlow + (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);

             xx = xupp - (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);
          }

          return (par[2]*step*sum*invsq2pi / par[3]);
}

void event_builder(int runN, int batch){
    bool dofit = true; 

    string fnamedata = Form("./data/%i.root", runN);
    
    ifstream crs_map_file(Form("data/batch_%i_map.dat", batch));
    map<int, int> crs_map; 

    int crs_id, crs_counter = 0; 
    while(crs_map_file >> crs_id){
        crs_counter ++; 
        crs_map[56-crs_counter] = crs_id; 
    }
    for(int ii = 56; ii < 60; ii++){
        crs_map[ii] = 8000+ii; 
    }
    for(int ii = 60; ii < 64; ii++){
        crs_map[ii] = 9000+ii; 
    }

    const size_t start = fnamedata.find_last_of("/\\") + 1;
    const size_t end   = fnamedata.find_last_of('.');

    const string basename = fnamedata.substr(start, end - start);

    string fname = "./out/" + basename + "_preprocessed.root";
    string outname  = "./pdf/" + basename + "_ana.pdf";
    string outroot = "./out/" + basename + "_evt.root"; 

    TFile *f = new TFile(fname.c_str(), "read"); 
    TTree *t = (TTree*) f->Get("tout"); 

    ULong64_t ti; 
    UShort_t ch; 
    double A, Q; 
    int ids; 
    ULong64_t ts; 

    t->SetBranchAddress("ch", &ch); 
    t->SetBranchAddress("t", &ti); 
    t->SetBranchAddress("A", &A); 
    t->SetBranchAddress("Q", &Q);
    t->SetBranchAddress("ids", &ids); 
    t->SetBranchAddress("ts", &ts); 

  // divide data in timeslices
  map<int, timeslice> data_stream; 
  // you need to fix this with the name of the tree you define when you read
  cout << "Divide data into timeslices" << endl;
  for(int ii = 0; ii < t->GetEntries(); ii++){
    t->GetEntry(ii); 
    
    hit thishit; 
    thishit.t = ts; 
    thishit.ch = ch; 
    thishit.Q = Q; 
    thishit.A = A; 

    if(data_stream.find(ids) != data_stream.end()){
      data_stream[ids].hit.push_back(thishit); 
    }else{
      timeslice new_slice; 
      data_stream[ids] = new_slice; 
      data_stream[ids].hit.push_back(thishit); 
    }
  }

  vector<event> evts; 

  double maxdt = deltat_evt; 
  TH1D* htdiff = new TH1D("htdiff", "t diff from seed", 1000, -maxdt, maxdt);
  TH1D* htdiff_crs[55];
  for(int crs = 0; crs < 55; crs++){
    htdiff_crs[crs] = new TH1D(Form("htdiff_crs_%i", crs), 
    Form("t diff from seed crs %i", crs), 1000, -maxdt, maxdt);
  }
  TH1D* htdiff_crs_all = new TH1D("htdiff_crs_all", 
    "t diff from seed all crs", 1000, -maxdt, maxdt);
  TH1D* htdiff_pl[8]; 
    for(int pl = 0; pl < 8; pl++){
        htdiff_pl[pl] = new TH1D(Form("htdiff_pl_%i", pl), 
        Form("t diff from seed pl %i", pl), 1000, -maxdt, maxdt);
    }
    TH1D* htdiff_pl_all = new TH1D("htdiff_pl_all", 
        "t diff from seed all pl", 1000, -maxdt, maxdt);
    TH1D* htdiff_trg = new TH1D("htdiff_trg", 
        "t diff from seed trg", 1000, -maxdt, maxdt);


    TH1D* hAseed = new TH1D("hAseed", "Leading seed A", 1000, 0, 400);
    TH1D* hAseed_second = new TH1D("hAseed_second", "Second seed A", 1000, 0, 400);
    TH1D* hAdiff = new TH1D("hAdiff", "Leading seed A - second seed A", 1000, 0, 400);
    TH1D* htseeddiff = new TH1D("htseeddiff", "t diff between seeds", 1000, -maxdt, maxdt);

    TH2D* hAdiff_vs_tdiff = new TH2D("hAdiff_vs_tdiff", 
      "Leading seed A - second seed A vs t diff", 1000, -maxdt, maxdt, 1000, 0, 400); 

    TH1D* hseedlessch = new TH1D("hseedlessch", "Channels without seed", 64, 0, 64);

  // process timeslices
  cout << "Process timeslices" << endl; 
  int counter = 0; 
  for(auto &slice_pair : data_stream){
    timeslice &slice = slice_pair.second; 

    cout << counter << endl; 
    counter ++; 
    vector<event> slice_evts;

    // search for trigger
    cout << "search for trigger" << endl; 
    for(auto hit : slice.hit){
      if(hit.ch == 0){
        event thisevt; 
        if(slice_evts.size() > 0) {
          ULong64_t lastseed = slice_evts[slice_evts.size()-1].tseed; 
          ULong64_t deltat = (hit.t > lastseed) ? hit.t-lastseed : lastseed-hit.t;
          if(deltat<deltat_evt) cout << "Two trg within the same event" << endl; 
        }
        thisevt.tseed = hit.t; 
        thisevt.chseed = 0;
        slice_evts.push_back(thisevt); 
      }
    }

    // search for weird slices 
    if(slice_evts.size() == 0) {
      continue; 
    }

    bool u[4] = {false, false, false, false}, d[4] = {false, false, false, false}; 
    for(auto hit : slice.hit){
      if(hit.ch > 55 && hit.A > 0.4){
        if(hit.ch < 60) d[hit.ch - 56] = true; 
        else u[hit.ch-60] = true; 
      }
    }

    int uon = 0, don = 0; 
    for(int ii = 0; ii < 4; ii++){
      if(u[ii]) uon ++; 
      if(d[ii]) don ++; 
    }

    if(uon > 1) cout << "More than one up on" << endl; 
    else if(don > 1) cout << "More than one down on " << endl; 
    if(don > 1 && uon > 1) cout << "More than a pair on" << endl; 


    // now match to plastic hit
    for(auto hit : slice.hit){
      if(hit.ch > 59) {
        // search for closest event
        int idhit = find_evt(hit.t, slice_evts); 
        if(idhit < 0) continue; 

      auto& thisevt = slice_evts[idhit];
      if(thisevt.chseed != 0){
            if( (((hit.t > thisevt.tseed) ? hit.t-thisevt.tseed : thisevt.tseed-hit.t) < deltat_evt)){
                htseeddiff->Fill(static_cast<Long64_t>(hit.t) - static_cast<Long64_t>(thisevt.tseed));
            if(thisevt.Aseed > hit.A) {
                hAseed_second->Fill(hit.A); 
                hAseed->Fill(thisevt.Aseed);
            } else {
                hAseed_second->Fill(thisevt.Aseed);
                hAseed->Fill(hit.A);
            }
            hAdiff->Fill(fabs(thisevt.Aseed - hit.A));  
            hAdiff_vs_tdiff->Fill(static_cast<Long64_t>(hit.t) - static_cast<Long64_t>(thisevt.tseed), 
              fabs(thisevt.Aseed - hit.A));
        }
        if(fabs(thisevt.Aseed - hit.A) > 0.8) {
                        cout << "Discarding second seed because A diff > 800 mV" << endl; 
                        continue; 
                    }else{
                        if(hit.A>thisevt.Aseed) {
        thisevt.tseed = hit.t; 
        thisevt.chseed = hit.ch;
        thisevt.Aseed = hit.A; }
                    }
    }

        // update the event seed 
        thisevt.tseed = hit.t; 
        thisevt.chseed = hit.ch;
        thisevt.Aseed = hit.A; 
      }
    }

    bool near_start = (slice_evts[0].tseed < deltat_evt && data_stream[slice_pair.first-1].hit.size() > 0);
    auto prev_it = data_stream.find(slice_pair.first - 1);
    if(prev_it != data_stream.end()){
        cout << "Error: I should look in the previous slice but i did not find it!" << endl; 
        continue; 
    }

    if(near_start){
        cout << "Seed near the start of the slice - tseed = " << slice_evts[0].tseed << endl; 
        cout << "Look at slice before" << endl; 
        auto& prev_slice = data_stream[slice_pair.first-1];
        
        auto& thisevt = slice_evts[0];
        thisevt.tseed = thisevt.tseed + 1e9; // add 1 us to avoid double counting with the previous slice


        for(auto hit : prev_slice.hit){
            if(hit.ch > 59) {

                ULong64_t thistseed = thisevt.tseed;
                ULong64_t deltat = (thistseed > hit.t) ? thistseed-hit.t : hit.t-thistseed;
                if(deltat < deltat_evt) {
                    if(thisevt.chseed != 0){
                    if( (((hit.t > thistseed) ? hit.t-thistseed : thistseed-hit.t) < deltat_evt)){
                        htseeddiff->Fill(static_cast<Long64_t>(hit.t) - static_cast<Long64_t>(thisevt.tseed));
                     }
                    if(thisevt.Aseed > hit.A) {
                        hAseed_second->Fill(hit.A); 
                        hAseed->Fill(thisevt.Aseed);
                    } else {
                    hAseed_second->Fill(thisevt.Aseed);
                    hAseed->Fill(hit.A);
                    }
                hAdiff->Fill(fabs(thisevt.Aseed - hit.A));  
                hAdiff_vs_tdiff->Fill(static_cast<Long64_t>(hit.t) - static_cast<Long64_t>(thistseed), 
                fabs(thisevt.Aseed - hit.A));

                // discard second seed if A diff > 800 mV
                    if(fabs(thisevt.Aseed - hit.A) > 0.8) {
                        cout << "Discarding second seed because A diff > 800 mV" << endl; 
                        continue; 
                    }else{
                        if(hit.A>thisevt.Aseed) {
        thisevt.tseed = hit.t; 
        thisevt.chseed = hit.ch;
        thisevt.Aseed = hit.A; }
                    }
                }   
            }else{
                // update the event seed 
                thisevt.tseed = hit.t; 
                thisevt.chseed = hit.ch;
                thisevt.Aseed = hit.A; 
            }

        }
    }
}
cout << "Check for plastic hit" << endl; 
   auto it = slice_evts.begin();
   while(it != slice_evts.end()) {
    if (it->chseed == 0) {
        std::cout << "No plastic in the event\n";
        for(auto hit : slice.hit){
            ULong64_t deltat = (hit.t > it->tseed) ? hit.t-it->tseed : it->tseed-hit.t;
          if(deltat<deltat_evt) hseedlessch->Fill(hit.ch); 
        }

        it = slice_evts.erase(it);
    } else {
        ++it;
    }
}

    // search for hits related to the event
    for(auto hit : slice.hit){
      int idhit = find_evt(hit.t, slice_evts);
        if(idhit < 0) continue; 

      auto& thisevt = slice_evts[idhit];

        ULong64_t thist = hit.t; 
        if(near_start) thist=thist+1e9; 

        ULong64_t deltat = (thist > thisevt.tseed) ? thist-thisevt.tseed : thisevt.tseed-thist; 
        Long64_t deltat_sign = static_cast<Long64_t>(thist) - static_cast<Long64_t>(thisevt.tseed);
        if(deltat<deltat_coinc) thisevt.hit.push_back(hit); 

        //if(deltat_sign < deltat_coinc) cout << thisevt.chseed << endl; 

        htdiff->Fill(deltat_sign);
        if(hit.ch > 55 && hit.ch < 64) {
            htdiff_pl_all->Fill(deltat_sign);
            htdiff_pl[hit.ch-56]->Fill(deltat_sign);
        }else if(hit.ch == 0){
            htdiff_trg->Fill(deltat_sign);
        }else if(hit.ch > 0 && hit.ch < 56){
            htdiff_crs_all->Fill(deltat_sign);
            htdiff_crs[hit.ch-1]->Fill(deltat_sign);
        }
    }

    if(near_start){
        auto prev_slice = data_stream.find(slice_pair.first-1);
        auto& thisevt = slice_evts[0]; 
        for(auto hit : prev_slice->second.hit){
            ULong64_t deltat = (hit.t > thisevt.tseed) ? hit.t-thisevt.tseed : thisevt.tseed-hit.t; 
            Long64_t deltat_sign = static_cast<Long64_t>(hit.t) - static_cast<Long64_t>(thisevt.tseed);
            if(deltat<deltat_coinc) thisevt.hit.push_back(hit); 

            //if(deltat_sign < deltat_coinc) cout << thisevt.chseed << endl; 

            htdiff->Fill(deltat_sign);
            if(hit.ch > 55 && hit.ch < 64) {
                htdiff_pl_all->Fill(deltat_sign);
                htdiff_pl[hit.ch-56]->Fill(deltat_sign);
            }else if(hit.ch == 0){
                htdiff_trg->Fill(deltat_sign);
            }else if(hit.ch > 0 && hit.ch < 56){
                htdiff_crs_all->Fill(deltat_sign);
                htdiff_crs[hit.ch-1]->Fill(deltat_sign);
            }
        }
    }

    // add events to the event list 
    for(auto thisevt : slice_evts){
        if(thisevt.hit.size() == 1 || thisevt.chseed == 0) {
          cout << "No hits in the event" << endl;
          continue;
        }
        evts.push_back(thisevt); 
    }
  }


TCanvas* cAseed = new TCanvas("cAseed", "cAseed", 1600, 600);
cAseed->Divide(2,2);
cAseed->cd(1);
cAseed->cd(1)->SetLogy();
hAseed->Draw("histo");
hAseed_second->SetLineColor(2);
hAseed_second->Draw("samehisto");
cAseed->cd(2);
cAseed->cd(2)->SetLogy();
hAdiff->Draw("histo");
cAseed->cd(3);
cAseed->cd(3)->SetLogz();
hAdiff_vs_tdiff->Draw("colz");
cAseed->cd(4);
cAseed->cd(4)->SetLogy();
htseeddiff->Draw("histo");
cAseed->Print((outname + "(").c_str());


    TCanvas* cseedless = new TCanvas("cseedless", "cseedless", 1600, 600);
    hseedlessch->Draw("histo");
    cseedless->Print(outname.c_str());

  TCanvas* ctdiff = new TCanvas("ctdiff", "ctdiff", 1600, 1600);
  ctdiff->Divide(2,2);
  ctdiff->cd(1); ctdiff->cd(1)->SetLogy(); 
  htdiff->SetLineColor(1);
  htdiff->Draw();
  htdiff_crs_all->SetLineColor(2); 
  htdiff_crs_all->Draw("same");
  htdiff_pl_all->SetLineColor(3);
  htdiff_pl_all->Draw("same");
  htdiff_trg->SetLineColor(4);
  htdiff_trg->Draw("same");
  ctdiff->cd(1)->BuildLegend();

  ctdiff->cd(2); ctdiff->cd(2)->SetLogy();
  htdiff_crs_all->Draw();
  for(int ii = 0; ii < 8; ii++){
    htdiff_crs[ii]->SetLineColor(ii+3); 
    htdiff_crs[ii]->Draw("same");
  }
  ctdiff->cd(2)->BuildLegend();

  ctdiff->cd(3); ctdiff->cd(3)->SetLogy();
    htdiff_pl_all->Draw();
    for(int ii = 0; ii < 4; ii++){
        htdiff_pl[ii]->SetLineColor(ii+4); 
        htdiff_pl[ii]->Draw("same");
}
ctdiff->cd(3)->BuildLegend();

ctdiff->Print(outname.c_str());

  TH1D* hQ[55][4]; 
  for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        hQ[crs][pl] = new TH1D(Form("hQ_%i_%i", crs, pl), 
        Form("Charge crs %i - trg %i; Q (nWb); A.U.", crs, pl), 30, 0, 4000); 
    }
  }

TFile * fout = new TFile(outroot.c_str(), "recreate");
int chout, seedout; 
double qout; 
TTree* tout = new TTree("t", "t"); 
tout->Branch("ch", &chout); 
tout->Branch("seed", &seedout); 
tout->Branch("q", &qout);  

  for(auto evt : evts){
    for(auto hit : evt.hit){
        int crs = hit.ch - 1; 
        int pl = evt.chseed-60; 
        if(crs < 55 && crs > -1) {
          chout = crs; 
          seedout = pl; 
          qout = hit.Q; 
          tout->Fill(); 
          hQ[crs][pl]->Fill(hit.Q); 
        }
    }
  }

  fout->cd(); 
  tout->Write(); 
  fout->Write(); 


TCanvas* c = new TCanvas(); 
c->Divide(6, 10); 
for(int crs = 0; crs < 55; crs ++){
    c->cd(crs+1); 
    hQ[crs][0]->Scale(1./hQ[crs][0]->GetBinContent(hQ[crs][0]->GetMaximumBin()));
    hQ[crs][0]->Draw("histo"); 
    for(int pl = 1; pl < 4; pl++){
        hQ[crs][pl]->SetLineColor(pl); 
        hQ[crs][pl]->Scale(1./hQ[crs][pl]->GetBinContent(hQ[crs][pl]->GetMaximumBin()));
        hQ[crs][pl]->Draw("samehisto"); 
    }
}

c->Print(outname.c_str()); 

for(int crs = 0; crs < 55; crs ++){
TCanvas* cqi = new TCanvas(); 
    for(int pl = 0; pl < 4; pl++){
        cqi->cd(); 
        hQ[crs][pl]->Draw("samehisto");
    }
    cqi->Print(outname.c_str()); 
}

if(!dofit) {
    return;
    c->Print("out.pdf]");
}

TCanvas* cfitlandau = new TCanvas(); 
cfitlandau->Divide(6, 10); 
TCanvas* cfitlangau = new TCanvas(); 
cfitlangau->Divide(6, 10); 

TCanvas* cfit = new TCanvas(); 

TF1* flandau[55][4]; 
TF1* flangau[55][4]; 
for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        cfitlandau->cd(1+pl*8+crs);
        double probs[2] = {0.2, 0.90};
        double q[2];
        hQ[crs][pl]->GetQuantiles(2, q, probs);
        flandau[crs][pl] = new TF1(Form("flandau_%i_%i", crs, pl), "landau", q[0], q[1]); 
        hQ[crs][pl]->Fit(flandau[crs][pl], "rq"); 
        hQ[crs][pl]->Draw(); 
        flandau[crs][pl]->Draw("samel"); 

        cfit->cd(); 
        hQ[crs][pl]->Draw(); 
        flandau[crs][pl]->Draw("samel"); 
        cfit->Print(outname.c_str()); 


        cfitlangau->cd(1+pl*8+crs);
        double mu = flandau[crs][pl]->GetParameter(1); 
        double sigma = flandau[crs][pl]->GetParameter(2); 
        flangau[crs][pl] = new TF1(Form("flangau_%i_%i", crs, pl), langau, mu-1.5*sigma, mu+5*sigma, 4);
        flangau[crs][pl]->SetParameters(sigma, mu, 1000, 2); 
        hQ[crs][pl]->Fit(flangau[crs][pl], "rq"); 
        hQ[crs][pl]->Draw(); 
        flangau[crs][pl]->Draw("samel"); 


        cfit->cd(); 
        hQ[crs][pl]->Draw(); 
        flangau[crs][pl]->Draw("samel"); 
        cfit->Print(outname.c_str()); 

    }
  }


        cfitlandau->Print(outname.c_str()); 
        cfitlangau->Print(outname.c_str()); 

 
  double Edep_1[4] = {40, 38, 34, 32};
  TCanvas* cg = new TCanvas(); 

  TGraphErrors* gr_att[55]; 
  TCanvas* catt = new TCanvas(); 
  catt->Divide(10, 6); 
  for(int crs = 0; crs < 55; crs++){
    gr_att[crs] = new TGraphErrors(); 
    gr_att[crs]->SetTitle(Form("Crs %i; x (cm); phe/MeV", crs+1)); 
    for(int pl = 0; pl < 4; pl++){
        double x = (crs < 5 || (crs > 9 && crs < 15) || (crs > 19 && crs < 20) || (crs > 29 && crs < 30) || (crs > 39 && crs < 40) || (crs > 49 && crs < 55) || (crs > 59 && crs < 65)) ? 3+pl*6 : 21-pl*6;

        double Ediff = 2*7; 
        if(crs < 5 || (crs > 9 && crs < 15) || (crs > 19 && crs < 20) || (crs > 29 && crs < 30) || (crs > 39 && crs < 40) || (crs > 49 && crs < 55) || (crs > 59 && crs < 65)) Ediff = Edep_1[pl]; 
        else{ Ediff = Edep_1[3-pl];}

        gr_att[crs]->SetPoint(pl, x, flangau[crs][pl]->GetParameter(1)/(Ediff)); 
        gr_att[crs]->SetPointError(pl, 1, flangau[crs][pl]->GetParError(1)/(Ediff));
    }
    catt->cd(crs+1); 
    gr_att[crs]->SetMarkerColor(20); 
    gr_att[crs]->Draw("AP"); 

    cg->cd();
    gr_att[crs]->Draw("ap"); 
    cg->Print(outname.c_str()); 
  }

  catt->Print((outname + ")").c_str()); 

}