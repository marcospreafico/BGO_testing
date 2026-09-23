void plot_calib_trend(){
    ifstream if1("./phe_260717.dat"); 
    ifstream if2("./phe_260722.dat"); 
    ifstream if3("./phe_260727.dat");
    
    double x, y; 
    TGraph* g1 = new TGraph(); g1->SetTitle("260717; ch; cc (nWb/phe)");
    TGraph* g2 = new TGraph(); g2->SetTitle("260722; ch; cc (nWb/phe)");
    TGraph* g3 = new TGraph(); g3->SetTitle("260724; ch; cc (nWb/phe)");

    TH1D* h1 = new TH1D("h1", "260717; cc (nWb/phe); A.U.", 50, 0, 180); 
    TH1D* h2 = new TH1D("h2", "260722; cc (nWb/phe); A.U.", 50, 0, 180);
    TH1D* h3 = new TH1D("h3", "260727; cc (nWb/phe); A.U.", 50, 0, 180);

    while(if1 >> x >> y){
        if(x > 53) break; 
        g1->SetPoint(x, x+1, y); 
        h1->Fill(y); 
    }
    while(if2 >> x >> y){
        if(x > 53) break; 
        g2->SetPoint(x, x+1, y); 
        h2->Fill(y); 
    }
    while(if3 >> x >> y){
        if(x > 53) break; 
        g3->SetPoint(x, x+1, y); 
        h3->Fill(y); 
    }

    TCanvas* cg = new TCanvas(); 
    g1->SetMarkerColor(1); 
    g1->SetMarkerStyle(20); 
    g1->SetLineColor(1); 
    g1->Draw("alp"); 

    g2->SetMarkerColor(2); 
    g2->SetMarkerStyle(20); 
    g2->SetLineColor(2); 
    g2->Draw("samelp"); 

    g3->SetMarkerColor(3); 
    g3->SetMarkerStyle(20); 
    g3->SetLineColor(3); 
    g3->Draw("samelp"); 

    TCanvas* ch = new TCanvas(); 
    h1->SetLineColor(1); h1->Draw("histo"); 
    h2->SetLineColor(2); h2->Draw("samehisto"); 
    h3->SetLineColor(3); h3->Draw("samehisto"); 
}
