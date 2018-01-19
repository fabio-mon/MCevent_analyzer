#include "EvAnalyz.hh"
#include "ConfigFile.hh"

#include <vector>
#include <string>

#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace std;

int GetBinNumber(TProfile* p,float x)
{
  TProfile p_aux(*p);
  return p_aux.Fill(x,0.);
}

int GetBinNumber2d(TProfile2D* p2, float x, float y)
{
  TProfile2D p2_aux(*p2);
  return p2_aux.Fill(x,y,0.);
}

void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity);

EvAnalyz::EvAnalyz(const ConfigFile & config)//:
//fconfig(config)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   vector<string> Filename;

   if(config.keyExists("Filename"))
   {
      config.readIntoVect(Filename,"Filename");
      if(Filename.size()==0)
      {
         cerr<<"[ERROR]: error while reading <Filename> parameter"<<endl;
         exit(EXIT_FAILURE);
      }
   }
   else
   {
      cerr<<"[ERROR]: <Filename> is a required parameter"<<endl;
      exit(EXIT_FAILURE);
   }

   cout<<"> Parsing config file"<<endl; 
   ParseConfigFile(config); 

   fDataTree = new TChain("digi","digi");
   int nfiles=0;
   for(std::vector<string>::iterator it = Filename.begin() ; it != Filename.end(); ++it)
      nfiles += fDataTree->Add(it->c_str());
   if(nfiles==0)
   {
      cerr<<"[ERROR]: empty tree"<<endl;
      exit(EXIT_FAILURE);
   }
   else
      cout<<"> "<<nfiles<<" file added to chain for a total of "<<fDataTree->GetEntries()<<" entries"<<endl;
   SetBranchTree();
   CreateProfile();
   CreateHisto();
   FillProfile();
   FillHisto();

}


//---------------------------------------------------------------------------------------------------------------
EvAnalyz::EvAnalyz(TChain* outtree, int Nthr, vector<float> thr, string DataLabel, float amp_min, float amp_max, float risetime_min, float risetime_max, float time_offset):
fDataTree(outtree),
fNthr(Nthr),
fthr(thr),
fDataLabel(DataLabel),
famp_min(amp_min),
famp_max(amp_max),
frisetime_min(risetime_min),
frisetime_max(risetime_max),
ftime_offset(time_offset)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   SetBranchTree();
   CreateProfile();
   CreateHisto();
   FillProfile();
   FillHisto();
}

//---------------------------------------------------------------------------------------------------------------
EvAnalyz::~EvAnalyz()
{
   cout<<"> Deleting profiles";
   for(int i=0; i<fNthr; i++)
   {
      delete fp_time_amp[fthr[i]];
      delete fp_time_risetime[fthr[i]];
      delete fp2_time_x_y[fthr[i]];
   }
   cout<<"OK"<<endl;

   cout<<"> Deleting histos";
   for(int i=0; i<fNthr; i++)
   {
      delete fh_time[fthr[i]];
   }
   cout<<"OK"<<endl;

   cout<<"> Deleting chain";
   delete fDataTree;
   cout<<"OK"<<endl;

   cout<<"> Deleting canvases";
   for(std::map<string,TCanvas*>::iterator it=fPlots.begin(); it!=fPlots.end(); ++it)
      if( it->second ) delete it->second;
   cout<<"OK"<<endl;
}


//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::ParseConfigFile(const ConfigFile & config)
{
   vector<string> Filename;
   config.readIntoVect(Filename,"Filename");   //check on filename has already been made

   if(config.keyExists("DataLabel"))
      fDataLabel = config.read<string>("DataLabel");
   else
   {
      TString Label_ToBeModified = Filename.at(0);
      Label_ToBeModified.ReplaceAll(".root","");
      Label_ToBeModified.ReplaceAll("*","");
      Label_ToBeModified.ReplaceAll(".","p");
      fDataLabel = Label_ToBeModified.Data();
   }

   if(config.keyExists("thr"))
      config.readIntoVect(fthr,"thr");
   else
   {
      cerr<<"[ERROR]: <thr> is a required parameter"<<endl;
      exit(EXIT_FAILURE);
   }
   fNthr = fthr.size();

   if(config.keyExists("amp_min"))
      famp_min = config.read<float>("amp_min");
   else
      famp_min = 0;

   if(config.keyExists("amp_max"))
      famp_max = config.read<float>("amp_max");
   else
      famp_max = 8000;

   if(config.keyExists("risetime_min"))
      frisetime_min = config.read<float>("risetime_min");
   else
      frisetime_min = 0;

   if(config.keyExists("risetime_max"))
      frisetime_max = config.read<float>("risetime_max");
   else
      frisetime_max = 10;

   if(config.keyExists("time_offset"))
      ftime_offset = config.read<float>("time_offset");
   else
      ftime_offset = 10;

}


//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::SetBranchTree()
{
   cout<<">> Branching the chain"<<endl;
   fDataTree->SetBranchStatus("*", 0);
   
   fDataTree -> SetBranchStatus("mu_x_hit",1); 
   fDataTree -> SetBranchStatus("mu_y_hit",1); 
   fDataTree -> SetBranchStatus("AMP_MAX",1); 
   fDataTree -> SetBranchAddress("mu_x_hit",&fmu_y_hit);
   fDataTree -> SetBranchAddress("mu_y_hit",&fmu_x_hit);
   fDataTree -> SetBranchAddress("AMP_MAX",&fAMP_MAX);

   for(int i=0; i<fNthr; i++)
   {
      fDataTree -> SetBranchStatus(Form("LDE%.0f",fthr[i]),1); 
      fDataTree -> SetBranchAddress(Form("LDE%.0f",fthr[i]),&ftime[fthr[i]]);
   }

   //fDataTree -> SetBranchStatus("PH2",1); fDataTree -> SetBranchAddress("PH2",&Phtime2);
   //fDataTree -> SetBranchStatus("PH5",1); fDataTree -> SetBranchAddress("PH5",&Phtime5);
   //fDataTree -> SetBranchStatus("PH10",1); fDataTree -> SetBranchAddress("PH10",&Phtime10);
   //fDataTree -> SetBranchStatus("PH20",1); fDataTree -> SetBranchAddress("PH20",&Phtime20);
   //fDataTree -> SetBranchStatus("PH50",1); fDataTree -> SetBranchAddress("PH50",&Phtime50);
   //fDataTree -> SetBranchStatus("PH100",1); fDataTree -> SetBranchAddress("PH100",&Phtime100);
}


//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::CreateProfile(bool mkamp, bool mkrisetime, bool mkpos)
{
   cout<<">> Creating time profiles"<<endl;
   for(int i=0; i<fNthr; i++)
   {
      if(mkamp)
         fp_time_amp[fthr[i]] = new TProfile(	Form("%s, time vs AMP_MAX, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						Form("%s, time vs AMP_MAX, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						100,famp_min,famp_max);
      if(mkrisetime)
         fp_time_risetime[fthr[i]] = new TProfile(Form("%s, time vs risetime(50-20), thr = %.0f ph",fDataLabel.c_str()/*,frisetime_min,frisetime_max*/,fthr[i]),
						Form("%s, time vs risetime(50-20), thr = %.0f ph",fDataLabel.c_str()/*,frisetime_min,frisetime_max*/,fthr[i]),
						100,frisetime_min,frisetime_max);
      if(mkpos)
      fp2_time_x_y[fthr[i]] = new TProfile2D(	Form("%s, time vs impact point, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						Form("%s, time vs impact point, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						22,-6.,6.,22,-6.,6.);
   }
}


//-------------------------------------------------------------------------------------------------------------
void EvAnalyz::CreateHisto()
{
   cout<<">> Creating time histogram"<<endl;
   for(int i=0; i<fNthr; i++)
      fh_time[fthr[i]] = new TH1F(	Form("%s, time distribution, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
					Form("%s, time distribution, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
					/*150*/200,-0.505,/*0.995*/1.495);
}

//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::FillProfile(bool mkamp, bool mkrisetime, bool mkpos)
{
   cout<<">> Filling profiles"<<endl;
   Long64_t nentries = fDataTree->GetEntries();
   for(Long64_t ientry=0; ientry<nentries; ientry++)
   {
      cout<<"\tReading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      for(int i=0; i<fNthr; i++)
      {
         if(mkamp)
            fp_time_amp[fthr[i]] -> Fill(fAMP_MAX,ftime[fthr[i]]-ftime_offset);
         if(mkrisetime)
            fp_time_risetime[fthr[i]] -> Fill(ftime[50]-ftime[20], ftime[fthr[i]]-ftime_offset);
         if(mkpos)
            fp2_time_x_y[fthr[i]] -> Fill(fmu_x_hit, fmu_y_hit, ftime[fthr[i]]-ftime_offset);
      }
   }
   cout<<"\n";
}


//---------------------------------------------------------------------------------------------------------------------
void EvAnalyz::FillHisto()
{
   cout<<">> Filling time histo"<<endl;;
   Long64_t nentries = fDataTree->GetEntries();
   for(Long64_t ientry=0; ientry<nentries; ientry++)
   {
      cout<<"\tReading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      for(int i=0; i<fNthr; i++)
          fh_time[fthr[i]]->Fill(ftime[fthr[i]]-ftime_offset);
   }
   cout<<"\n";
}
//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::DrawProfiles(float time_min, float time_max)
{
   cout<<"> Drawing time profiles"<<endl;
//creating canvas
   fPlots[fDataLabel+", time vs AMP_MAX"] = new TCanvas( Form("%s, time vs AMP_MAX",fDataLabel.c_str()) , Form("%s, time vs AMP_MAX",fDataLabel.c_str()) );
   fPlots[fDataLabel+", time vs risetime"] = new TCanvas( Form("%s, time vs risetime",fDataLabel.c_str()) , Form("%s, time vs risetime",fDataLabel.c_str()) );

//creating legends
   TLegend leg_time_amp(0.1,0.7,0.48,0.9);
   TLegend leg_time_risetime(0.1,0.7,0.48,0.9);

//creating title object
   TLatex title;
   title.SetNDC();

   for(int i=0; i<fNthr; i++)
   {
      fPlots[Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i])] = new TCanvas( Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(), fthr[i]) , Form("%s, time vs impact point",fDataLabel.c_str()) ,400,400);
   }

//esthetical set
   for(int i=0; i<fNthr; i++)
   {
      fp_time_amp[fthr[i]]->SetLineColor(i+1);  
      fp_time_risetime[fthr[i]]->SetLineColor(i+1);  
   }

//draw time vs amp_max & vs risetime
   fPlots[fDataLabel+", time vs AMP_MAX"]->cd();
   fp_time_amp[fthr[0]]->Draw("");  
   leg_time_amp.AddEntry(fp_time_amp[fthr[0]],Form("thr = %.0f ph",fthr[0]),"l");
   fp_time_amp[fthr[0]]->GetYaxis()->SetRangeUser(time_min,time_max);
   fp_time_amp[fthr[0]]->GetXaxis()->SetTitle("amp max (ph)"); 
   fp_time_amp[fthr[0]]->GetYaxis()->SetTitle("time (ns)");  
 
   fPlots[fDataLabel+", time vs risetime"]->cd();
   fp_time_risetime[fthr[0]]->Draw("");  
   leg_time_risetime.AddEntry(fp_time_risetime[fthr[0]],Form("thr = %.0f ph",fthr[0]),"l");
   fp_time_risetime[fthr[0]]->GetYaxis()->SetRangeUser(time_min,time_max);
   fp_time_risetime[fthr[0]]->GetXaxis()->SetTitle("risetime 20-50(ns)"); 
   fp_time_risetime[fthr[0]]->GetYaxis()->SetTitle("time (ns)"); 

   for(int i=1; i<fNthr; i++)
   {
      fPlots[fDataLabel+", time vs AMP_MAX"]->cd();
      fp_time_amp[fthr[i]]->Draw("same"); 
      leg_time_amp.AddEntry(fp_time_amp[fthr[i]],Form("thr = %.0f ph",fthr[i]),"l");
      fPlots[fDataLabel+", time vs risetime"]->cd();
      fp_time_risetime[fthr[i]]->Draw("same"); 
      leg_time_risetime.AddEntry(fp_time_risetime[fthr[i]],Form("thr = %.0f ph",fthr[i]),"l");
   } 
   fPlots[fDataLabel+", time vs AMP_MAX"]->cd();
   leg_time_amp.Draw();
   title.DrawLatex(0.1,0.93,(fDataLabel+", time vs AMP_MAX").c_str());

   fPlots[fDataLabel+", time vs risetime"]->cd();
   leg_time_risetime.Draw();
   title.DrawLatex(0.1,0.93,(fDataLabel+", time vs risetime").c_str());

//draw time vs impact point
   for(int i=0; i<fNthr; i++)
   {
      fPlots[Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i])] -> cd();
      fp2_time_x_y[fthr[i]] -> Draw("COLZ");
      fp2_time_x_y[fthr[i]] -> GetZaxis()->SetRangeUser(time_min,time_max);
      fp2_time_x_y[fthr[i]] -> GetXaxis()->SetTitle("x (mm)"); 
      fp2_time_x_y[fthr[i]] -> GetYaxis()->SetTitle("y (mm)"); 
      fp2_time_x_y[fthr[i]] -> GetZaxis()->SetTitle("time (ns)"); 
      title.DrawLatex(0.1,0.93,Form("%s, time vs impact point, thr = %.0f",fDataLabel.c_str(),fthr[i]));
   }

//print time profiles
   TString plotname = fDataLabel+", time vs AMP_MAX";
   plotname.ReplaceAll(" ","_");
   plotname.ReplaceAll("=","");
   plotname.ReplaceAll(",","_");
   plotname.ReplaceAll(".","p");
   plotname+= ".pdf";
   fPlots[fDataLabel+", time vs AMP_MAX"]->Print(plotname.Data());
   plotname.ReplaceAll("pdf","png");
   fPlots[fDataLabel+", time vs AMP_MAX"]->Print(plotname.Data());

   plotname = fDataLabel+", time vs risetime";
   plotname.ReplaceAll(" ","_");
   plotname.ReplaceAll("=","");
   plotname.ReplaceAll(",","_");
   plotname.ReplaceAll(".","p");
   plotname+= ".pdf";
   fPlots[fDataLabel+", time vs risetime"]->Print(plotname.Data());
   plotname.ReplaceAll("pdf","png");
   fPlots[fDataLabel+", time vs risetime"]->Print(plotname.Data());

   for(int i=0; i<fNthr; i++)
   {
      plotname = Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i]);
      plotname.ReplaceAll(" ","_");
      plotname.ReplaceAll("=","");
      plotname.ReplaceAll(",","_");
      plotname.ReplaceAll(".","p");
      plotname+= ".pdf";
      fPlots[Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i])]->Print(plotname.Data());
      plotname.ReplaceAll("pdf","png");
      fPlots[Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i])]->Print(plotname.Data());
   }

}


//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::DrawHistos()
{
   cout<<"> Drawing time histos"<<endl;

//creating title object
   TLatex title;
   title.SetNDC();
   //creating canvas
   for(int i=0; i<fNthr; i++)
   {
      fPlots[Form("%s, time histo, thr=%.0f",fDataLabel.c_str(),fthr[i])] = new TCanvas( Form("%s, time histo, thr=%.0f",fDataLabel.c_str(), fthr[i]) , Form("%s, time histo",fDataLabel.c_str()));
   }

//esthetical set
   for(int i=0; i<fNthr; i++)
      fh_time[fthr[i]]->SetLineColor(i+1);  

//draw time vs impact point
   for(int i=0; i<fNthr; i++)
   {
      fPlots[Form("%s, time histo, thr=%.0f",fDataLabel.c_str(),fthr[i])] -> cd();
      fh_time[fthr[i]] -> Draw();
      fh_time[fthr[i]] -> GetXaxis()->SetTitle("time (ns)"); 
      //fp2_time_x_y[fthr[i]] -> GetZaxis()->SetRangeUser(time_min,time_max);
      title.DrawLatex(0.1,0.93,Form("%s, time distribution, thr = %.0f",fDataLabel.c_str(),fthr[i]));
   }

//print time histos
   TString plotname;
   for(int i=0; i<fNthr; i++)
   {
      plotname = Form("%s, time histo, thr=%.0f",fDataLabel.c_str(),fthr[i]);
      plotname.ReplaceAll(" ","_");
      plotname.ReplaceAll("=","");
      plotname.ReplaceAll(",","_");
      plotname.ReplaceAll(".","p");
      plotname+= ".pdf";
      fPlots[Form("%s, time histo, thr=%.0f",fDataLabel.c_str(),fthr[i])]->Print(plotname.Data());
      plotname.ReplaceAll("pdf","png");
      fPlots[Form("%s, time histo, thr=%.0f",fDataLabel.c_str(),fthr[i])]->Print(plotname.Data());
   }

}




//---------------------------------------------------------------------------------------------------------------
EvAnalyz EvAnalyz::AmpCorrection()
{
   cout<<"> Amplitude walk correction"<<endl;
   std::map<int,TF1*> fitamw;
   for(int i=0;i<fNthr;i++)
   {
      fitamw[fthr[i]] = new TF1(Form("amplitude walk correction, thr = %.0f",fthr[i]),"[2]+[0]*exp(-[1]*x)",famp_min,famp_max);
      fitamw[fthr[i]]->SetLineWidth(1);
      fitamw[fthr[i]]->SetLineColor(1);
   }
   fitamw[fthr[fNthr-1]]->SetParameters(0.01,0.0006,1.);
   for(int i=fNthr-1;i>-1;i--)
   {
      fp_time_amp[fthr[i]]->Fit(fitamw[fthr[i]],"R");
      if(i>=1)
      {
         fitamw[fthr[i-1]]->SetParameter(0,fitamw[fthr[i]]->GetParameter(0));
         fitamw[fthr[i-1]]->SetParameter(1,fitamw[fthr[i]]->GetParameter(1));
         fitamw[fthr[i-1]]->SetParameter(2,fitamw[fthr[i]]->GetParameter(2));
      }
   }

   //branch the new tree
   cout<<">> Branching new tree"<<endl;
   float mu_y_hit, mu_x_hit, AMP_MAX;
   std::map<float,float> time;
   TFile outfile(("/tmp/"+fDataLabel+"_amw.root").c_str(),"RECREATE");
   TTree *outtree=new TTree("digi",(fDataLabel+" amplitude walk corrected").c_str());
   //outtree->Branch("eventNb",&eventNb,"eventNb/I");
   outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
   outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");
   outtree->Branch("AMP_MAX",&AMP_MAX,"AMP_MAX/F");
   for(int i=0; i<fNthr; i++)
      outtree->Branch( Form("LDE%.0f",fthr[i]) , &time[fthr[i]] , Form("LDE%.0f/F",fthr[i]) );

   
   //amplitude correction 
   cout<<">> Filling new tree"<<endl;
   Long64_t nentries = fDataTree->GetEntries();
   for (Long64_t ientry=0;ientry<nentries;ientry++)
   {
      cout<<"Reading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      mu_x_hit = fmu_x_hit;
      mu_y_hit = fmu_y_hit;         
      AMP_MAX = fAMP_MAX;
      for(int i=0;i<fNthr;i++)
         time[fthr[i]] = ftime[fthr[i]] - ftime_offset - fitamw[fthr[i]]->Eval(AMP_MAX);

      outtree->Fill();//Fill the output ntuple
   }
   outtree->Write();
   outfile.Close();
   TChain* outchain = new TChain("digi","digi amplitude walk corrected");
   outchain->Add(("/tmp/"+fDataLabel+"_amw.root").c_str());
   //create the new EvAnalyz
   cout<<">> Creating "<<fDataLabel<<"_amw"<<endl;
   EvAnalyz data_amw(outchain, fNthr, fthr, fDataLabel+"_amw", famp_min, famp_max, frisetime_min, frisetime_max, 0./*ftime_offset=0*/);
   return data_amw;

}


//---------------------------------------------------------------------------------------------------------------
EvAnalyz EvAnalyz::MitigatedAmpCorrection(float amp_min_fit, float amp_max_fit)
{
   cout<<"> Mitigated amplitude walk correction"<<endl;
   std::map<int,TF1*> fitamw;
   for(int i=0;i<fNthr;i++)
   {
      fitamw[fthr[i]] = new TF1(Form("mitigated amplitude walk correction, thr = %.0f",fthr[i]),"[0] + [1]*log([2]*x)",amp_min_fit,amp_max_fit);
      fitamw[fthr[i]]->SetLineWidth(1);
      fitamw[fthr[i]]->SetLineColor(1);
   }
   fitamw[fthr[fNthr-1]]->SetParameters(10.4,-0.3,0.00013);
   for(int i=fNthr-1;i>-1;i--)
   {
      fp_time_amp[fthr[i]]->Fit(fitamw[fthr[i]],"R");
      if(i>=1)
      {
         fitamw[fthr[i-1]]->SetParameter(0,fitamw[fthr[i]]->GetParameter(0));
         fitamw[fthr[i-1]]->SetParameter(1,fitamw[fthr[i]]->GetParameter(1));
         fitamw[fthr[i-1]]->SetParameter(2,fitamw[fthr[i]]->GetParameter(2));
      }
   }

   //branch the new tree
   cout<<">> Branching new tree"<<endl;
   float mu_y_hit, mu_x_hit, AMP_MAX;
   std::map<float,float> time;
   TFile outfile(("/tmp/"+fDataLabel+"_mitigatedamw.root").c_str(),"RECREATE");
   TTree *outtree=new TTree("digi",(fDataLabel+" mitigated amplitude walk corrected").c_str());
   //outtree->Branch("eventNb",&eventNb,"eventNb/I");
   outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
   outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");
   outtree->Branch("AMP_MAX",&AMP_MAX,"AMP_MAX/F");
   for(int i=0; i<fNthr; i++)
      outtree->Branch( Form("LDE%.0f",fthr[i]) , &time[fthr[i]] , Form("LDE%.0f/F",fthr[i]) );

   
   //amplitude correction 
   cout<<">> Filling new tree"<<endl;
   Long64_t nentries = fDataTree->GetEntries();
   for (Long64_t ientry=0;ientry<nentries;ientry++)
   {
      cout<<"Reading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      mu_x_hit = fmu_x_hit;
      mu_y_hit = fmu_y_hit;         
      AMP_MAX = fAMP_MAX;
      for(int i=0;i<fNthr;i++)
         time[fthr[i]] = ftime[fthr[i]] - ftime_offset - fitamw[fthr[i]]->Eval(AMP_MAX);

      outtree->Fill();//Fill the output ntuple
   }
   outtree->Write();
   outfile.Close();
   TChain* outchain = new TChain("digi","digi mitigated amplitude walk corrected");
   outchain->Add(("/tmp/"+fDataLabel+"_mitigatedamw.root").c_str());
   //create the new EvAnalyz
   cout<<">> Creating "<<fDataLabel<<"_mitigatedamw"<<endl;
   EvAnalyz data_amw(outchain, fNthr, fthr, fDataLabel+"_mitigatedamw", famp_min, famp_max, frisetime_min, frisetime_max, 0./*ftime_offset=0*/);
   return data_amw;

}


//---------------------------------------------------------------------------------------------------------------------------
EvAnalyz EvAnalyz::PosCorrection()
{
   cout<<"> Position correction"<<endl;

   //branch the new tree
   cout<<">> Branching new tree"<<endl;
   float mu_y_hit, mu_x_hit, AMP_MAX;
   std::map<float,float> time;
   TFile outfile(("/tmp/"+fDataLabel+"_poscorr.root").c_str(),"RECREATE");
   TTree *outtree=new TTree("digi",(fDataLabel+" impact point correction").c_str());
   //outtree->Branch("eventNb",&eventNb,"eventNb/I");
   outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
   outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");
   outtree->Branch("AMP_MAX",&AMP_MAX,"AMP_MAX/F");
   for(int i=0; i<fNthr; i++)
      outtree->Branch( Form("LDE%.0f",fthr[i]) , &time[fthr[i]] , Form("LDE%.0f/F",fthr[i]) );

   //position correction 
   cout<<">> Filling new tree"<<endl;
   Long64_t nentries = fDataTree->GetEntries();
   for (Long64_t ientry=0;ientry<nentries;ientry++)
   {
      cout<<"Reading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      mu_x_hit = fmu_x_hit;
      mu_y_hit = fmu_y_hit;         
      AMP_MAX = fAMP_MAX;
      for(int i=0;i<fNthr;i++)
         time[fthr[i]] = ftime[fthr[i]] - ftime_offset - fp2_time_x_y[fthr[i]]->GetBinContent(GetBinNumber2d(fp2_time_x_y[fthr[i]],fmu_x_hit,fmu_y_hit));

      outtree->Fill();//Fill the output ntuple
   }
   outtree->Write();
   outfile.Close();
   TChain* outchain = new TChain("digi","digi impact point corrected");
   outchain->Add(("/tmp/"+fDataLabel+"_poscorr.root").c_str());
   //create the new EvAnalyz
   cout<<">> Creating "<<fDataLabel<<"_poscorr"<<endl;
   EvAnalyz data_poscorr(outchain, fNthr, fthr, fDataLabel+"_poscorr", famp_min, famp_max, frisetime_min, frisetime_max, 0./*ftime_offset=0*/);
   return data_poscorr;

}

EvAnalyz EvAnalyz::RiseTimeCorrection()
{
   cout<<"> Risetime correction"<<endl;

   //branch the new tree
   cout<<">> Branching new tree"<<endl;
   float mu_y_hit, mu_x_hit, AMP_MAX;
   std::map<float,float> time;
   TFile outfile(("/tmp/"+fDataLabel+"_risetimecorr.root").c_str(),"RECREATE");
   TTree *outtree=new TTree("digi",(fDataLabel+" risetime correction").c_str());
   //outtree->Branch("eventNb",&eventNb,"eventNb/I");
   outtree->Branch("mu_x_hit",&mu_x_hit,"mu_x_hit/F");
   outtree->Branch("mu_y_hit",&mu_y_hit,"mu_y_hit/F");
   outtree->Branch("AMP_MAX",&AMP_MAX,"AMP_MAX/F");
   for(int i=0; i<fNthr; i++)
      outtree->Branch( Form("LDE%.0f",fthr[i]) , &time[fthr[i]] , Form("LDE%.0f/F",fthr[i]) );

   //risetime correction 
   cout<<">> Filling new tree"<<endl;
   Long64_t nentries = fDataTree->GetEntries();
   for (Long64_t ientry=0;ientry<nentries;ientry++)
   {
      cout<<"Reading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      mu_x_hit = fmu_x_hit;
      mu_y_hit = fmu_y_hit;         
      AMP_MAX = fAMP_MAX;
      for(int i=0;i<fNthr;i++)
         time[fthr[i]] = ftime[fthr[i]] - ftime_offset - fp_time_risetime[fthr[i]]->GetBinContent(GetBinNumber(fp_time_risetime[fthr[i]],ftime[50]-ftime[20]));

      outtree->Fill();//Fill the output ntuple
   }
   outtree->Write();
   outfile.Close();
   TChain* outchain = new TChain("digi","digi risetime corrected");
   outchain->Add(("/tmp/"+fDataLabel+"_risetimecorr.root").c_str());
   //create the new EvAnalyz
   cout<<">> Creating "<<fDataLabel<<"_risetimecorr"<<endl;
   EvAnalyz data_risetimecorr(outchain, fNthr, fthr, fDataLabel+"_risetimecorr", famp_min, famp_max, frisetime_min, frisetime_max, 0./*ftime_offset=0*/);
   return data_risetimecorr;

}

void EvAnalyz::SetAmpRange(float amp_min,float amp_max)
{
   famp_min=amp_min;
   famp_max=amp_max;
   cout<<"> Updating time vs amp profile"<<endl;
   for(int i=0; i<fNthr; i++)
   {
      delete fp_time_amp[fthr[i]];
   }
   CreateProfile(true,false,false);
   FillProfile(true,false,false);
}


void EvAnalyz::SetRiseTimeRange(float risetime_min,float risetime_max)
{
   frisetime_min=risetime_min;
   frisetime_max=risetime_max;
   cout<<"> Updating time vs risetime profile"<<endl;
   for(int i=0; i<fNthr; i++)
      delete fp_time_risetime[fthr[i]];

   CreateProfile(false,true,false);
   FillProfile(false,true,false);
}

TGraphErrors* EvAnalyz::ThrScan(std::string option)
{

   TGraphErrors* res_thr = new TGraphErrors();
   res_thr->SetName((fDataLabel+"_res_thr").c_str());
   res_thr->SetTitle((fDataLabel+"_res_thr").c_str());
   TF1* fitfunc = new TF1("time distribution fit","gaus(0)",-0.505,0.995);
   float* vals = new float[4];
   float min;
   float max;
   float SmallInt;
   for(int i=0; i<fNthr; i++)
   {
      if(option=="RMS" || option=="rms" || option=="Rms")
      {
         res_thr->SetPoint(i,fthr[i],fh_time[fthr[i]]->GetRMS());
         res_thr->SetPointError(i,0.,fh_time[fthr[i]]->GetRMSError());
      }
      else
         if(option=="FIT" || option=="fit" || option=="Fit")
         {
            fitfunc->SetParameter(1,fh_time[fthr[i]]->GetMean());
            fitfunc->SetParameter(2,fh_time[fthr[i]]->GetRMS());
            fh_time[fthr[i]]->Fit(fitfunc);
            res_thr->SetPoint(i,fthr[i],fitfunc->GetParameter(2));
            res_thr->SetPointError(i,0.,fitfunc->GetParError(2));
         }
         else
            if(option=="SMALLESTINTERVAL" || option=="smallestinterval" || option=="SmallestInterval")
            {
               FindSmallestInterval(vals,fh_time[fthr[i]],0.68,true); 
               min = vals[2];
               max = vals[3];
               SmallInt = 0.5*(max-min);
               res_thr->SetPoint(i,fthr[i],SmallInt);
               res_thr->SetPointError(i,0.,fh_time[fthr[i]]->GetRMSError());
            }
            else
            {
               cout<<"[ERROR]: Option "<<option<<" not valid"<<endl;
               break;
            } 
   }
   if(vals) delete[] vals;
   return res_thr;   
}


void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity)
{
  float integralMax = fraction * histo->Integral();
  
  int N = histo -> GetNbinsX();
  int M1 = 0;
  int M2 = 0;
  for(int bin1 = 0; bin1 < N; ++bin1)
  {
    if( histo->GetBinContent(bin1+1) > 0. && M1 == 0 ) M1 = bin1-1;
    if( histo->GetBinContent(bin1+1) > 0. ) M2 = bin1+2;
  }
  
  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    binCenters[bin1] = histo->GetBinCenter(bin1+1);
    binContents[bin1] = histo->GetBinContent(bin1+1);
    
    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += binContents[bin2];
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(int bin1 = M1; bin1 < M2; ++bin1)
  {
    for(int bin2 = bin1+1; bin2 < M2; ++bin2)
    {
      if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
      float tmpMin = histo -> GetBinCenter(bin1+1);
      float tmpMax = histo -> GetBinCenter(bin2+1);
      
      if( (tmpMax-tmpMin) < delta )
      {
        delta = (tmpMax - tmpMin);
        min = tmpMin;
        max = tmpMax;
      }
      
      break;
    }
  }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
  {
    if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);
    
    if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
  }
  smallHisto -> SetFillColor(kYellow);
  
  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();  
  
  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = min;
  ret[3] = max;
}





