ULong64_t time_slice_length = 1e9; 

void preprocess(){
    string fname = "./data/260818.root";
    
    int batch = 2; 

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

    const size_t start = fname.find_last_of("/\\") + 1;
    const size_t end   = fname.find_last_of('.');

    const string basename = fname.substr(start, end - start);

    string outname = "./out/" + basename + "_preprocessed.root";
    string outpdf  = "./pdf/" + basename + "_preprocess.pdf";

    string calibdata = "./calib/phe_"+basename+".dat"; 

    map<int, double> cc; 
    ifstream calibfile(calibdata.c_str()); 
    if (calibfile.is_open()) {
    double ch, calib;
    while (calibfile >> ch >> calib) {
        cc[ch] = calib; 
    }
} else {
    std::cout << "Calibration file " << calibdata
              << " not found. Skipping calibration." << std::endl;
    for(int ii = 0; ii < 63; ii ++){
        cc[ii] = 0; 
    }
}


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

    TH2D* hwf[63]; 
    TH2D* hwf_weird[63]; 
    TH1D* hq[63]; 
    TH1D* hA[63]; 
    TH2D* hAQ[63];
    TH2D* hAtot[63]; 
    for(int ii = 0; ii < 63; ii++){
        hwf[ii] = new TH2D(Form("hwf_%i", ii), Form("Crs %i - %i; t (ns); A (V)", ii+1, crs_map[ii+1]), 
                nsamples, 0, 8*nsamples, 1000, -0.5, 2);
        hwf_weird[ii] = new TH2D(Form("hwf_weird_%i", ii), Form("Crs %i - %i; t (ns); A (V)", ii+1, crs_map[ii+1]), 
                nsamples, 0, 8*nsamples, 1000, -0.5, 2);
        hq[ii] = new TH1D(Form("hq_%i", ii), Form("Q crs %i - %i; Q (pWb); A.U.", ii+1, crs_map[ii+1]), 1000, 0, 0); 
        hA[ii] = new TH1D(Form("hA_%i", ii), Form("A crs %i - %i; A (V); A.U.", ii+1, crs_map[ii+1]), 100, 0, 0); 
        hAQ[ii] = new TH2D(Form("hAQ_%i", ii), Form("Q vs A Crs %i - %i; Q (phe); A (V)", ii+1, crs_map[ii+1]), 
                1000, 0, 0, 1000, 0, 0);
        hAtot[ii] = new TH2D(Form("hAtot_%i", ii), Form("tot vs A Crs %i - %i; tot (ns); A (V)", ii+1, crs_map[ii+1]), 
                1000, 0, 0, 1000, 0, 0);
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

        double cci = cc[ch]; 
        ids = (int)(ti/time_slice_length);
        ts = ti-ids*time_slice_length; 

        int idch = ch-1; 

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
            if(y>ymax) ymax = y; 
            if(y<ymin) ymin = y; 
            if(jj>1 && (baseline-wf->GetAt(jj))/16384.<thr && y > thr) tot = 8; ; 
            if(y>thr) tot+=8; 
            if(ch!=0) hwf[idch]->Fill(8*jj, y);  
            Q += (y*1000*8);
        }

        A = ymax-ymin;

        if(A>Q/1400){
            for(int jj = 0; jj < nsamples; jj++){
            double y = (baseline-wf->GetAt(jj))/16384.;
            if(ch!=0) hwf_weird[idch]->Fill(8*jj, y);  
        }
        }

        if(idch < 55) Q = Q/cc[idch]; 
        else if(idch > 54 && idch < 63) A = A/cc[idch]; 



        if(ch != 0){ 
        hA[idch]->Fill(A); 
        hq[idch]->Fill(Q); 
        hAQ[idch]->Fill(Q, A); 
        hAtot[idch]->Fill(tot, A); 
        }
        tout->Fill(); 
    }
    
    fout->cd(); 
    tout->Write(); 
    fout->Write();


    TCanvas* cdummy = new TCanvas(); 
    cdummy->Print(Form("%s(", outpdf.c_str())); 

    TCanvas* canvout[63]; 

    for(int ii = 0; ii < 63; ii++){
        canvout[ii] = new TCanvas(); 
        canvout[ii]->Divide(3, 2); 
        canvout[ii]->cd(1); 
        canvout[ii]->cd(1)->SetLogz(); 
        hwf[ii]->Draw("colz"); 
        canvout[ii]->cd(2); 
        canvout[ii]->cd(2)->SetLogz(); 
        hwf_weird[ii]->Draw("colz"); 
        canvout[ii]->cd(3); 
        canvout[ii]->cd(3)->SetLogy(); 
        hq[ii]->Draw("histo"); 
        if(ii < 55) hq[ii]->GetXaxis()->SetTitle("Q (phe)"); 
        canvout[ii]->cd(4); 
        canvout[ii]->cd(4)->SetLogy(); 
        hA[ii]->Draw("histo"); 
        if(ii > 54) hA[ii]->GetXaxis()->SetTitle("A (phe)"); 
        canvout[ii]->cd(5);
        canvout[ii]->cd(5)->SetLogz(); 
        hAQ[ii]->Draw("colz"); 
        canvout[ii]->cd(6); 
        canvout[ii]->cd(6)->SetLogz(); 
        hAtot[ii]->Draw("colz"); 
        canvout[ii]->Print(outpdf.c_str()); 
    }

    cdummy->Print(Form("%s)", outpdf.c_str())); 
}
