#include "EvAnalyz.hh"
#include "ConfigFile.hh"

#include <vector>

#include "TString.h"
#include "TCanvas.h"
//#include <utility>
//using std::string;
using namespace std;

EvAnalyz::EvAnalyz(const ConfigFile & config)
{
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
   FillProfile();

}


//---------------------------------------------------------------------------------------------------------------
EvAnalyz::~EvAnalyz()
{
   cout<<"> Deleting profiles"<<endl;
   for(int i=0; i<fNthr; i++)
   {
      delete fp_time_amp[fthr[i]];
      delete fp_time_risetime[fthr[i]];
      delete fp2_time_x_y[fthr[i]];
   }

   cout<<"> Deleting chain"<<endl;
   delete fDataTree;

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
void EvAnalyz::CreateProfile()
{
   cout<<">> Creating time profiles"<<endl;
   for(int i=0; i<fNthr; i++)
   {
      fp_time_amp[fthr[i]] = new TProfile(	Form("%s, time vs AMP_MAX, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						Form("%s, time vs AMP_MAX, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						100,famp_min,famp_max);

      fp_time_risetime[fthr[i]] = new TProfile(	Form("%s, time vs risetime(%.0f-%.0f), thr = %.0f ph",fDataLabel.c_str(),frisetime_min,frisetime_max,fthr[i]),
						Form("%s, time vs risetime(%.0f-%.0f), thr = %.0f ph",fDataLabel.c_str(),frisetime_min,frisetime_max,fthr[i]),
						100,frisetime_min,frisetime_max);

      fp2_time_x_y[fthr[i]] = new TProfile2D(	Form("%s, time vs impact point, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						Form("%s, time vs impact point, thr = %.0f ph",fDataLabel.c_str(),fthr[i]),
						22,-6.,6.,22,-6.,6.);
   }
}


//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::FillProfile()
{
   Long64_t nentries = fDataTree->GetEntries();
   for(Long64_t ientry=0; ientry<nentries; ientry++)
   {
      cout<<"Reading entry "<<ientry<< "\r" << std::flush;
      fDataTree->GetEntry(ientry);
      for(int i=0; i<fNthr; i++)
      {
         fp_time_amp[fthr[i]] -> Fill(fAMP_MAX,ftime[fthr[i]]-ftime_offset);
         fp_time_risetime[fthr[i]] -> Fill(ftime[50]-ftime[20], ftime[fthr[i]]-ftime_offset);
         fp2_time_x_y[fthr[i]] -> Fill(fmu_x_hit, fmu_y_hit, ftime[fthr[i]]-ftime_offset);
      }
   }
}



//---------------------------------------------------------------------------------------------------------------
void EvAnalyz::DrawProfiles()
{
   cout<<"> Drawing time profiles"<<endl;
//creating canvas
   fPlots[fDataLabel+", time vs AMP_MAX"] = new TCanvas( Form("%s, time vs AMP_MAX",fDataLabel.c_str()) , Form("%s, time vs AMP_MAX",fDataLabel.c_str()) );
   fPlots[fDataLabel+", time vs risetime"] = new TCanvas( Form("%s, time vs risetime",fDataLabel.c_str()) , Form("%s, time vs risetime",fDataLabel.c_str()) );
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
   fPlots[fDataLabel+", time vs risetime"]->cd();
   fp_time_risetime[fthr[0]]->Draw("");  
   for(int i=1; i<fNthr; i++)
   {
      fp_time_amp[fthr[i]]->Draw("same"); 
      fp_time_risetime[fthr[i]]->Draw("same"); 
   } 

//draw time vs impact point
   for(int i=0; i<fNthr; i++)
   {
      fPlots[Form("%s, time vs impact point, thr=%.0f",fDataLabel.c_str(),fthr[i])] -> cd();
      fp2_time_x_y[fthr[i]] -> Draw("COLZ");
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


