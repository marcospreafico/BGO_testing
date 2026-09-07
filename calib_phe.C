ULong64_t time_slice_length = 1e9; 

void calib_phe(){
    string fname = "./data/phe_260807.root";
    
    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');

    const string basename = fname.substr(start, end - start);

    string outname = "./out/" + basename + "_preprocessed.root";
    string outpdf  = "./pdf/" + basename + "_preprocess.pdf";
    string outtxt = "./calib/" + basename + ".dat"; 

    ofstream ofile(outtxt.c_str());

    TFile *f = new TFile(fname.c_str(), "read"); 
    TTree *t = (TTree*) f->Get("Data_R");  

    TArrayS* wf = nullptr; 
    ULong64_t ti; 
    UShort_t ch; 

    t->SetBranchAddress("Samples", &wf); 
    t->SetBranchAddress("Timestamp", &ti); 
    t->SetBranchAddress("Channel", &ch); 

    t->GetEntry(0); 
    int nsamples = wf->GetSize(); 
    int nbaseline = 5; 

    TH1D* hq_range[63]; 
    TH1D* hA_range[63]; 
    for(int ii = 0; ii < 63; ii++){
        hq_range[ii] = new TH1D(Form("hq_range_%i", ii), Form("Q crs %i; Q (pWb); A.U.", ii), 100, 0, 0); 
        hA_range[ii] = new TH1D(Form("hA_range_%i", ii), Form("A crs %i; A (V); A.U.", ii), 50, 0, 0); 
    }

    for(int ii = 0; ii < t->GetEntries(); ii++){
         t->GetEntry(ii);

        int idch = ch-1; 
        if(idch < 0) continue; 

        double baseline = 0; 
        for(int jj = 0; jj < nbaseline; jj++){
            baseline += wf->GetAt(jj)/nbaseline; 
        }

        double Q = 0; 
        double ymin = 1e9, ymax = -1e9; 
        double thr = 0.01; 
        for(int jj = 0; jj < nsamples; jj++){
            double y = (baseline-wf->GetAt(jj))/16384.;
            if(8*jj > 60 && 8*jj < 300) {
                Q += (y*1000*8);
                if(y>ymax) ymax = y; 
                if(y<ymin) ymin = y; 
            }
        }

        double A = ymax-ymin;
        if(ch > 55 && A < 0.001) continue;

        if(ch != 0){ 
        hA_range[idch]->Fill(A); 
        hq_range[idch]->Fill(Q); 
        }
    }


    TH2D* hwf[63]; 
    TH1D* hq[63]; 
    TH1D* hA[63]; 
    TF1* fgaus[63]; 
    for(int ii = 0; ii < 63; ii++){
        hwf[ii] = new TH2D(Form("hwf_%i", ii), Form("Crs %i; t (ns); A (V)", ii+1), 
                nsamples, 0, 8*nsamples, 1000, -0.001, 0.01);
        double qrange[2] = {0, 0}; 
        hq_range[ii]->GetQuantiles(2, qrange, (double[2]){0.05, 0.995});
        hq[ii] = new TH1D(Form("hq_%i", ii), Form("Q crs %i; Q (pWb); A.U.", ii), 100, qrange[0], qrange[1]); 
        double arange[2] = {0, 0}; 
        hA_range[ii]->GetQuantiles(2, arange, (double[2]){0.02, 0.9999});
        hA[ii] = new TH1D(Form("hA_%i", ii), Form("A crs %i; A (V); A.U.", ii), 50, arange[0], arange[1]); 
    }

    TFile* fout = new TFile(outname.c_str(), "recreate"); 
    TTree* tout = new TTree("tout", "tout"); 
    double A, Q; 
    int ids; 
    ULong64_t ts; 

    tout->Branch("ch", &ch); 
    tout->Branch("t", &ti); 
    tout->Branch("A", &A); 
    tout->Branch("Q", &Q);
    tout->Branch("ids", &ids); 
    tout->Branch("ts", &ts); 
    
    for(int ii = 0; ii < t->GetEntries(); ii++){
        t->GetEntry(ii);

        ids = (int)(ti/time_slice_length);
        ts = ti-ids*time_slice_length; 

        int idch = ch-1; 
        if(idch < 0) continue; 

        double baseline = 0; 
        for(int jj = 0; jj < nbaseline; jj++){
            baseline += wf->GetAt(jj)/nbaseline; 
        }

        Q = 0; 
        double ymin = 1e9, ymax = -1e9; 
        int tot = 0; 
        double thr = 0.01; 
        for(int jj = 0; jj < nsamples; jj++){
            double y = (baseline-wf->GetAt(jj))/16384.;
            if(jj>1 && (baseline-wf->GetAt(jj))/16384.<thr && y > thr) tot = 8; ; 
            if(y>thr) tot+=8; 
            if(ch!=0) hwf[idch]->Fill(8*jj, y);  
            if(8*jj > 60 && 8*jj < 300) {
                Q += (y*1000*8);
                if(y>ymax) ymax = y; 
                if(y<ymin) ymin = y; 
            }
        }

        A = ymax-ymin;

        if(ch != 0){ 
        hA[idch]->Fill(A); 
        hq[idch]->Fill(Q); 
        }
        tout->Fill(); 
    }

    TCanvas* cdummy = new TCanvas(); 
    cdummy->Print(Form("%s(", outpdf.c_str())); 
    
    TCanvas* c[63]; 

    TGraphErrors* glin[63];
    TF1* flin[63];
    TH1D* phe_distr = new TH1D("hdistr", "Phe value distribution; Q phe (nWb); counts", 50, 0, 0); 

    TGraph* gaphe = new TGraph(); 

    ifstream ifile_param("./calib/fit_param.dat");
    map <int, double*> fit_param;
    double idch, p; 
    while(ifile_param >> idch){
        fit_param[idch] = new double[9];
        for(int jj = 0; jj < 9; jj++){
            ifile_param >> p; 
            fit_param[idch][jj] = p; 
        }
    }

    for(int ii = 0; ii < 63; ii++){
        cout << ii << endl; 
        c[ii] = new TCanvas(); 
        c[ii]->Divide(2,2); 
        c[ii]->cd(1); 
        c[ii]->cd(1)->SetLogz(); 
        hwf[ii]->Draw("colz"); 
        c[ii]->cd(2); 
        c[ii]->cd(2)->SetLogy(); 
        hA[ii]->Draw("histo"); 
        if(ii > 54){
        double Q1 = hA[ii]->GetBinCenter(hA[ii]->GetMaximumBin());
        cout << Q1 << endl; 
        double A1 = hA[ii]->GetBinContent(hA[ii]->GetMaximumBin());
        fgaus[ii] = new TF1(Form("fgaus_%i", ii), "gaus(0)+gaus(3)+gaus(6)+[9]*exp(x/[10])", Q1-0.001, 3*Q1+0.002);
        fgaus[ii]->SetParameters(A1, Q1, 0.0002, A1/10, 2*Q1, 0.0002, A1/100, 3*Q1, 0.0003, 5000, -0.002);
        hA[ii]->Fit(fgaus[ii], "R"); 
        fgaus[ii]->Draw("samel");
        }
        c[ii]->cd(3); 
        c[ii]->cd(3)->SetLogy(); 
        cout << "Q" << endl;  
        hq[ii]->Draw("histo");
        if(ii < 55){
        double Q1 = hq[ii]->GetBinCenter(hq[ii]->GetMaximumBin());
        double A1 = hq[ii]->GetBinContent(hq[ii]->GetMaximumBin());
        fgaus[ii] = new TF1(Form("fgaus_%i", ii), "gaus(0)+gaus(3)+gaus(6)", 50, 3*Q1+50);
        // for(int jj = 0; jj < 9; jj++){
        //     fgaus[ii]->SetParameter(jj, fit_param[ii][jj]); 
        // }
        fgaus[ii]->SetParameters(A1, Q1, 24, A1/10, 2*Q1, 24, A1/100, 3*Q1, 24);
        // 1 phe
    fgaus[ii]->SetParLimits(1, 0.7*Q1, 1.3*Q1);
    fgaus[ii]->SetParLimits(2, 10, 50);

    // 2 phe
    fgaus[ii]->SetParLimits(4, 1.5*Q1, 2.5*Q1);
    fgaus[ii]->SetParLimits(5, 10, 50);
            // 3 phe
    fgaus[ii]->SetParLimits(7, 2.5*Q1, 3.5*Q1);
    fgaus[ii]->SetParLimits(8, 10, 50);

        hq[ii]->Fit(fgaus[ii], "R"); 
        fgaus[ii]->Draw("samel");
        }
        c[ii]->cd(4); 
        glin[ii] = new TGraphErrors(); 
        glin[ii]->SetTitle("Phe linear fit; n phe; Q phe (nWb)");
        for(int jj = 0; jj < 3; jj++){
            glin[ii]->SetPoint(jj, jj+1, fgaus[ii]->GetParameter(3*jj+1));
            glin[ii]->SetPointError(jj, 0, fgaus[ii]->GetParError(3*jj+1)); 
        }
        flin[ii] = new TF1(Form("flin_%i", ii), "pol1", 0, 4);
        glin[ii]->Fit(flin[ii]); 
        glin[ii]->SetMarkerStyle(20); 
        cout << "lin" << endl;  
        glin[ii]->Draw("ap"); 
        flin[ii]->Draw("samel");
        double Aphe = flin[ii]->GetParameter(1); 
        ofile << ii << " " << Aphe << endl; 
        if(ii < 55) gaphe->SetPoint(gaphe->GetN(), ii, Aphe);
        if(ii < 55) phe_distr->Fill(Aphe);  
        c[ii]->Print(outpdf.c_str()); 
    }

    ofstream ofile_param("./calib/fit_param.dat");
    for(int ii = 0; ii < 63; ii++){
        ofile_param << ii << " "; 
        for(int jj = 0; jj < 9; jj++){
            ofile_param << fgaus[ii]->GetParameter(jj) << " "; 
        }
        ofile_param << endl;
    }

    TCanvas* caphe = new TCanvas(); 
    gaphe->SetMarkerStyle(20); 
    gaphe->Draw("AP"); 
    double amean = phe_distr->GetMean(); 
    TLine* l1 = new TLine(0, amean, 60, amean); 
    TLine* l2 = new TLine(0, 1.1*amean, 60, 1.1*amean); 
    TLine* l3 = new TLine(0, 0.9*amean, 60, 0.9*amean); 

    l1->SetLineColor(2); 
    l1->Draw("samel"); 
    l2->SetLineStyle(3); 
    l2->Draw("samel"); 
    l3->SetLineStyle(3); 
    l3->Draw("samel"); 


    TCanvas* cequalization = new TCanvas();
    phe_distr->Draw("histo");
    cequalization->Print(outpdf.c_str()); 

    cdummy->Print(Form("%s)", outpdf.c_str())); 

    fout->cd(); 
    tout->Write(); 
    fout->Write();

}
