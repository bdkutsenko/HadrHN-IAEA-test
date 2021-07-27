// pictures XS
// 
// 06.07.2012


string sInt(Int_t number)
  {
    stringstream ss;
    ss << number;
    return ss.str();
  }

void plotXS()
{
  cout <<  "Started plotXS() " << endl;
  int verbose = 0;
  gROOT->SetStyle("Plain");
  gROOT->Reset();
  static const G4int NVAR=16;
  static const G4int NVAR1=48;
  G4String type[NVAR] = {"Total p + p","Total n + p","Total pi+ + p",
			 "Total pi- + p","Total p + n","Total n + n",
			 "Total pi+ + n","Total pi- + n",
			 "Elastic p + p","Elastic n + p","Elastic pi+ + p",
			 "Elastic pi- + p","Elastic p + n","Elastic n + n",
			 "Elastic pi+ + n","Elastic pi- + n"};
  G4String num[NVAR1] = {"0","1","2","3","4","5","6","7",
			"8","9","10","11","12","13","14","15",
			"16","17","18","19","20","21","22","23",
			"24","25","26","27","28","29","30","31",
			"32","33","34","35","36","37","38","39",
			"40","41","42","43","44","45","46","47"};


  int ee[5] = {0, 1, 10, 100, 250};
  double Xmin[5] = {0,-1,-1,-1,-1};
  double Xmax[5] = {8,1,1,1,1};
  double Ymin[5] = {0.00015, 0., 1.1, 1.1, 1.1};
  double Ymax[5] = {.015,1000,10000,10000,10000};
  double Y1max[101]={0,
  .1,.1,.15,.15,.5,.5,.5,.5,.5,.5,  // 1-10
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 11-20
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 21-30
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 31-40
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 41-50
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 51-60
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 61-70
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 71-80
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15,  // 81-90
  .15,.15,.15,.15,.15,.15,.15,.15,.15,.15   // 91-100
  };

  TLegend* leg;
 
  TString datafile="test/xsect.root";
  TString Name = "Total and elastic cross sections p + p in mb";
  TFile* f = new TFile(datafile);
  if(!f->IsOpen()) {
    cout <<  "File " << datafile  << " is opened" << endl;
    exit 1;
  }
  if(0 <= verbose) cout <<  "File " << datafile  << " is opened" << endl;

  TCanvas* c1 = new TCanvas("c1",Name,-1, 7, 900, 600);

  TH1F* hhh = gPad->DrawFrame(-1, 0.0001, 7, 10000. ,"p + p cross sections in mb");
  hhh->GetXaxis()->SetTitle("log10(E/MeV)");
  hhh->SetTitle(Name);
  hhh->Draw("SAME");
  gPad->SetLogy(1);
  gPad->Update();
  leg = new TLegend(.5, .2, 0.7, .4);
  leg->SetTextSize(.03);
  leg->SetHeader("Cross sections"); 

  TH1F* hh1 = (TH1F*)f->Get("0");
    hh1->SetLineColor(2);
    hh1->SetLineWidth(3);
    hh1->SetMarkerSize(.2);
    hh1->SetMarkerStyle(20);
    hh1->SetMarkerColor(2);
    leg->AddEntry(hh1, "Penelope", "pl");
    hh1->Draw("HIST SAME PL");

    TH1F* hh2 = (TH1F*)f->Get("h2");
    hh2->SetLineColor(3);
    hh2->SetLineWidth(3.);
    hh2->SetMarkerSize(.2);
    hh2->SetMarkerStyle(21);
    hh2->SetMarkerColor(3);
    leg->AddEntry(hh2, "Livermore", "pl");
    hh2->Draw("HIST SAME PL");
    
    TH1F* hh3 = (TH1F*)f->Get("h3");
    hh3->SetMarkerSize(.2);
    hh3->SetLineWidth(2);
    hh3->SetLineColor(4);
    hh3->SetMarkerStyle(22);
    hh3->SetMarkerColor(4);
    leg->AddEntry(hh3, "Bethe-Heitler-5D", "pl");
    leg->Draw("SAME");
    hh3->Draw("HIST SAME PL");

    TH1F* hh4 = (TH1F*)f->Get("h4");
    hh4->SetMarkerSize(.2);
    hh4->SetLineColor(6);
    hh4->SetLineWidth(2);
    hh4->SetMarkerStyle(22);
    hh4->SetMarkerColor(6);
    leg->AddEntry(hh4, "Bethe-Heitler", "pl");
    leg->Draw("SAME");
    hh4->Draw("HIST SAME PL");

    TH1F* hh5 = (TH1F*)f->Get("h5");
    hh5->SetMarkerSize(.2);
    hh5->SetLineWidth(1);
    hh5->SetLineColor(46);
    //    hh5->SetLineStyle(3);
    hh5->SetMarkerStyle(22);
    hh5->SetMarkerColor(46);
    leg->AddEntry(hh5, "Relativistic", "pl");
    leg->Draw("SAME");
    hh5->Draw("HIST SAME PL");
    /*
    TH1F* hh6 = (TH1F*)f.Get("h6");
    hh6->SetMarkerSize(.2);
    hh6->SetLineColor(7);
    hh6->SetLineWidth(1.);
    //    hh6->SetLineStyle(3);
    hh6->SetMarkerStyle(22);
    hh6->SetMarkerColor(7);
    leg->AddEntry(hh6, "Boldyshev", "pl");
    leg->Draw("SAME");
    hh6->Draw("HIST SAME PL");
    */

    cv->Update();
    name = "_" + s;
    cv->Print("plotXS/ApicXS" + name + ".png");
    //    gPad->Close(); 
    cv->Close();
  }
}
