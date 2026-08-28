
gErrorIgnoreLevel = kWarning;

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

double langauExpo(double *x, double *par)
{
    // par[0] = Landau width
    // par[1] = MPV
    // par[2] = area
    // par[3] = Gaussian sigma
    // par[4] = exponential amplitude
    // par[5] = exponential decay constant

    double lang = langau(x, par);

    double expo = par[4] * TMath::Exp(-x[0] / par[5]);

    return lang + expo;
}

void ana_fit(int batch){

    gStyle->SetOptFit(1111); 

    string fname = Form("./out/batch%i.root", batch);
    string outname  = Form("./pdf/batch%i_fit.pdf", batch);

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

    TFile *f = new TFile(fname.c_str(), "read"); 
    TTree *t = (TTree*) f->Get("t"); 

    int ch, seed; 
    double q; 

    t->SetBranchAddress("ch", &ch); 
    t->SetBranchAddress("seed", &seed); 
    t->SetBranchAddress("q", &q); 
    
  TH1D* hQ_range[55][4]; 
  for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        hQ_range[crs][pl] = new TH1D(Form("hQ_range_%i_%i", crs, pl), 
        Form("Charge crs %i - trg %i; Q (phe); A.U.", crs, pl), 25,  200, 4000); 
    }
  }

  for(int ii = 0; ii < t->GetEntries(); ii++){
    t->GetEntry(ii); 
    hQ_range[ch][seed]->Fill(q); 
  }

 
TH1D* hQ[55][4]; 
  for(int crs = 0; crs < 55; crs++){
    double xmin = 4000, xmax = 0; 

    for(int pl = 0; pl < 4; pl++){
        double probs[2] = {0.05, 0.95};
        double q[2];

        hQ_range[crs][pl]->GetQuantiles(2, q, probs);

        if(q[0] < xmin) xmin = q[0];
        if(q[1] > xmax) xmax = q[1];
    }

    for(int pl = 0; pl < 4; pl++){
        hQ[crs][pl] = new TH1D(Form("hQ_%i_%i", crs, pl), 
        Form("Charge crs %i - trg %i; Q (phe); A.U.", crs_map[crs+1], pl), 25 , xmin, xmax); 
        hQ[crs][pl]->SetLineColor(pl+1); 
    }
  }

  TH2D* hcount[4]; 
  for(int ii = 0; ii < 4; ii++){
    hcount[ii] = new TH2D(Form("hcount_%i", ii), Form("Count trg %i", ii), 5, 0, 5, 12, 0, 12); 
  }

  for(int ii = 0; ii < t->GetEntries(); ii++){
    t->GetEntry(ii); 
    
    int x = 4 - ch % 5, y = (int) ch / 5; 
    hcount[seed]->Fill(x, y); 
    
    hQ[ch][seed]->Fill(q); 
  }

TCanvas* c = new TCanvas(); 
c->Divide(7, 8); 
for(int crs = 0; crs < 55; crs ++){
    c->cd(crs+1); 
    hQ_range[crs][0]->Scale(1./hQ_range[crs][0]->GetBinContent(hQ_range[crs][0]->GetMaximumBin()));
    hQ_range[crs][0]->Draw("histo"); 
    for(int pl = 1; pl < 4; pl++){
        hQ_range[crs][pl]->SetLineColor(pl); 
        hQ_range[crs][pl]->Scale(1./hQ_range[crs][pl]->GetBinContent(hQ_range[crs][pl]->GetMaximumBin()));
        hQ_range[crs][pl]->Draw("samehisto"); 
    }
}

c->Print(Form("%s(", outname.c_str())); 

for(int crs = 0; crs < 55; crs ++){
TCanvas* cqi = new TCanvas(); 
    for(int pl = 0; pl < 4; pl++){
        cqi->cd(); 
        hQ_range[crs][pl]->Draw("samehisto");
    }
    cqi->Print(outname.c_str()); 
}


TCanvas* cfitlandau[4]; 
for(int ii = 0; ii < 4; ii++){
    cfitlandau[ii] = new TCanvas(); 
    cfitlandau[ii]->Divide(7, 8); 
}
TCanvas* cfitlangau[4]; 
for(int ii = 0; ii < 4; ii++){
    cfitlangau[ii] = new TCanvas(); 
    cfitlangau[ii]->Divide(7, 8); 
}


TCanvas* cfit = new TCanvas(); 
cfit->Divide(2,1); 

TF1* flandau[55][4]; 
TF1* flangau[55][4]; 

for(int crs = 0; crs < 55; crs++){
    for(int pl = 0; pl < 4; pl++){
        cfitlandau[pl]->cd(1+crs);

TH1D *hist = hQ[crs][pl];
int peakBin = 2;
double maxSmooth = -1;
for(int ib = 2; ib < hist->GetNbinsX(); ib++){
    double smooth =
        hist->GetBinContent(ib-1) +
        hist->GetBinContent(ib)   +
        hist->GetBinContent(ib+1);
    if(smooth > maxSmooth){
        maxSmooth = smooth;
        peakBin = ib;
    }
}

double peak = 2000; 
for(int ii = 3; ii < hQ[crs][pl]->GetNbinsX(); ii++){
    double x = hQ_range[crs][pl]->GetBinCenter(ii);
    double y = hQ_range[crs][pl]->GetBinContent(ii);
    if(y > hQ_range[crs][pl]->GetBinContent(hQ_range[crs][pl]->FindBin(peak))) peak = x;
}

double rms = hQ[crs][pl]->GetRMS();

// Range più largo
double fitLow  = std::max(0.0, peak - 1.5*rms);
double fitHigh = peak + 2.5*rms;


// ============================================================
// LANDAU + EXPONENTIAL
// ============================================================

flandau[crs][pl] = new TF1(
    Form("flandau_%i_%i", crs, pl),
    "landau(0) + expo(3)",
    fitLow,
    fitHigh
);

// Landau parameters
flandau[crs][pl]->SetParameter(0, 0.5*hQ[crs][pl]->GetEntries());
flandau[crs][pl]->SetParameter(1, peak);
flandau[crs][pl]->SetParameter(2, 100);

flandau[crs][pl]->SetParLimits(0, 0, 1000);
flandau[crs][pl]->SetParLimits(1, peak-rms, peak+rms);
flandau[crs][pl]->SetParLimits(2, 10, 300);

// Exponential initial parameters
// expo(3) = exp(p3 + p4*x)
flandau[crs][pl]->SetParameter(3, TMath::Log(std::max(1.0, hQ[crs][pl]->GetBinContent(hQ[crs][pl]->FindBin(fitLow)))));

flandau[crs][pl]->SetParameter(4, -1.0/rms);
flandau[crs][pl]->SetParLimits(4, -0.1, -1e-6);

hQ[crs][pl]->Fit(flandau[crs][pl], "QBL");

hQ[crs][pl]->Draw();
flandau[crs][pl]->Draw("same");

cfit->cd(1);

hQ[crs][pl]->Draw();
flandau[crs][pl]->Draw("same");

cfitlangau[pl]->cd(1+crs);


double Afit = flandau[crs][pl]->GetParameter(0) * 100.0;

double mu =  flandau[crs][pl]->GetParameter(1);

double sigma = flandau[crs][pl]->GetParameter(2);


// ============================================================
// LANDAU-GAUSSIAN + EXPONENTIAL
// ============================================================

// Range ancora più largo rispetto a prima
double langauLow  = std::max(0.0, mu - 3.0*sigma);
double langauHigh = mu + 6.0*sigma;

flangau[crs][pl] = new TF1(
    Form("flangau_%i_%i", crs, pl),
    langauExpo,
    langauLow,
    langauHigh,
    6
);

flangau[crs][pl]->SetNpx(1000);


// Langau initial parameters
flangau[crs][pl]->SetParameter(0, sigma);
flangau[crs][pl]->SetParameter(1, mu);
flangau[crs][pl]->SetParLimits(1, mu - 200, mu + 200);
flangau[crs][pl]->SetParameter(2, Afit);
flangau[crs][pl]->SetParameter(3, 100);


// Exponential initial parameters
double expAmplitude =
    hQ[crs][pl]->GetBinContent(
        hQ[crs][pl]->FindBin(langauLow)
    );

flangau[crs][pl]->SetParameter(
    4,
    std::max(1.0, expAmplitude)
);

flangau[crs][pl]->SetParameter(
    5,
    std::max(100.0, rms)
);


// Parameter limits
flangau[crs][pl]->SetParLimits(0, 10, 300);

flangau[crs][pl]->SetParLimits(
    1,
    mu - 200,
    mu + 200
);

flangau[crs][pl]->SetParLimits(3, 20, 300);

// Exponential amplitude
flangau[crs][pl]->SetParLimits(
    4,
    0,
    1e6
);

// Exponential decay scale
flangau[crs][pl]->SetParLimits(
    5,
    20,
    5000
);


hQ[crs][pl]->Fit(flangau[crs][pl], "QBL");

hQ[crs][pl]->Draw();
flangau[crs][pl]->Draw("same");


cfit->cd(2);

hQ[crs][pl]->Draw();
flangau[crs][pl]->Draw("same");

cfit->Print(outname.c_str());

    }
  }

  for(int ii = 0; ii < 4; ii++){
    //cfitlandau[ii]->Print(outname.c_str()); 
    cfitlangau[ii]->Print(outname.c_str()); 
  }

 
  TCanvas* cg = new TCanvas(); 
  cg->Divide(2,1); 

  TH1D* hly_distrib = new TH1D("hly_distrib", "Distribution of light yield; phe; A.U.", 100, 0, 0);

  TGraphErrors* gr_att[55]; 
  TGraphErrors* gr_ly[55]; 
  TGraphErrors* gr_ly_norm[55]; 

  double H[8] = {4.96, 5.31, 5.27, 6.65, 6.28, 6.08, 5.80, 5.75};
  double h[8] = {2.44, 2.38, 2.13, 2.45, 2.08, 1.88, 1.74, 1.69};
  double Ediff = 2*7; 

  TCanvas* catt = new TCanvas(); 

  vector<double> ly_map(55); 

  catt->Divide(10, 6); 
  for(int crs = 0; crs < 55; crs++){
    gr_att[crs] = new TGraphErrors(); 
    gr_att[crs]->SetTitle(Form("Crs %i; x (cm); phe", crs_map[crs+1])); 

    gr_ly[crs] = new TGraphErrors();
    gr_ly[crs]->SetTitle(Form("Crs %i; x (cm); phe/MeV", crs_map[crs+1]));

    gr_ly_norm[crs] = new TGraphErrors();
    gr_ly_norm[crs]->SetTitle(Form("Crs %i; x (cm); phe/MeV (normalized to average)", crs_map[crs+1]));

    double ly_average = 0; 
    for(int pl = 0; pl < 4; pl++){
      
        if(hQ[crs][pl]->GetEntries() < 20) continue;
        double x = (crs < 5 || 
                    (crs > 9 && crs < 15) ||
                     (crs > 19 && crs < 25) || 
                     (crs > 29 && crs < 35) || 
                     (crs > 39 && crs < 45) || 
                     (crs > 49 && crs < 55) ) ? 3+pl*6 : 21-pl*6;

        int type = (int) crs_map[crs+1]/100;
        if(type > 8) type = 16-type; 

        double  y = H[type-1] + (h[type-1]-H[type-1])/24*x;

        gr_att[crs]->SetPoint(pl, x, flangau[crs][pl]->GetParameter(1)); 
        gr_att[crs]->SetPointError(pl, 1, flangau[crs][pl]->GetParError(1));

        gr_ly[crs]->SetPoint(pl, x, flangau[crs][pl]->GetParameter(1)/(Ediff*y));
        gr_ly[crs]->SetPointError(pl, 1, flangau[crs][pl]->GetParError(1)/(Ediff*y));

        gr_ly_norm[crs]->SetPoint(pl, x, flangau[crs][pl]->GetParameter(1)/(Ediff*y));
        gr_ly_norm[crs]->SetPointError(pl, 1, flangau[crs][pl]->GetParError(1)/(Ediff*y));

        ly_average += flangau[crs][pl]->GetParameter(1)/(Ediff*y);
    }

    hly_distrib->Fill(ly_average/4.);

    ly_map[crs] = ly_average/4.;


    catt->cd(crs+1); 
    gr_att[crs]->SetMarkerColor(20); 
    gr_att[crs]->Draw("AP"); 

    cg->cd(1);
    gr_att[crs]->Draw("ap"); 
    cg->cd(2);
    gr_ly[crs]->Draw("ap");
    cg->Print(outname.c_str()); 
  }


  TCanvas* cly_distrib = new TCanvas(); 
  hly_distrib->Draw("histo");
  cly_distrib->Print(outname.c_str());

  gStyle->SetOptStat(0); 

TCanvas* ccount = new TCanvas(); 
ccount->Divide(2,2); 
for(int ii = 0; ii < 4; ii++){
  ccount->cd(ii+1); 
  //ccount->cd(ii+1)->SetLogz(); 
  hcount[ii]->Draw("histo"); 
}
ccount->Print(outname.c_str()); 

TCanvas* cly_type[8]; 
for(int ii = 0; ii < 8; ii++){
  cly_type[ii] = new TCanvas(); 
  cly_type[ii]->Divide(2,1); 
}

TH1D* hly_type[8];
for(int ii = 0; ii < 8; ii++){
  hly_type[ii] = new TH1D(Form("hly_type_%i", ii), Form("Distribution of light yield type %i; phe/MeV; A.U.", ii+1), 100, 0, 0);
}

for(int ii = 0; ii < 55; ii++){
  int type = (int) crs_map[ii+1]/100;
  if(type > 8) type = 16-type; 

  hly_type[type-1]->Fill(ly_map[ii]);

  cly_type[type-1]->cd(1);
  gr_ly_norm[ii]->SetMarkerStyle(20);
  gr_ly_norm[ii]->SetMarkerColor(hly_type[type-1]->GetEntries());
  for(int jj = 0; jj < 4; jj++){
    gr_ly_norm[ii]->SetPoint(jj, jj, gr_ly_norm[ii]->GetPointY(jj)/ly_map[ii]);
    gr_ly_norm[ii]->SetPointError(jj, 0, gr_ly_norm[ii]->GetErrorY(jj)/ly_map[ii]);
  }

  if(hly_type[type-1]->GetEntries() == 1) gr_ly_norm[ii]->Draw("AP");
  else gr_ly_norm[ii]->Draw("P same");
}

for(int ii = 0; ii < 8; ii++){
  cly_type[ii]->cd(2);
  hly_type[ii]->Draw("histo");
  cly_type[ii]->Print(outname.c_str());
}


ofstream ofile_ly(Form("out/LY_batch_%i.dat", batch));

for(int ii = 0; ii < 55; ii++){
    ofile_ly << crs_map[ii+1] << "\t" ; 
  gr_ly[ii]->Sort();

    for(int jj = 0; jj < 4; jj ++){
      ofile_ly << gr_ly[ii]->GetPointY(jj) << "\t";
    } 
    ofile_ly << ly_map[ii] << endl; 
}


  catt->Print((outname + ")").c_str()); 

}