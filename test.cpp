#include <iostream>
#include "ConfigFile.hh"
#include "EvAnalyz.hh"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"

int main (int argc, char **argv)
{
   if(argc!=2)
   {
      cout<<"ERROR: unvalid number of input parameters\n";
      exit(EXIT_FAILURE);
   }
   ConfigFile config(argv[1]);
   
   bool interactive;
   if(config.keyExists("interactive"))
      interactive = config.read<bool>("interactive");
   else
      interactive = false;

   EvAnalyz data(config);
   EvAnalyz data_amw = data.AmpCorrection();

   TGraphErrors *gr_rms = data_amw.ThrScan("rms");
   TGraphErrors *gr_fit = data_amw.ThrScan("fit");
   TGraphErrors *gr_smallint = data_amw.ThrScan("smallestinterval");

   gr_rms->SetMarkerStyle(20);
   gr_fit->SetMarkerStyle(20);
   gr_smallint->SetMarkerStyle(20);

   gr_rms->SetMarkerColor(1);
   gr_fit->SetMarkerColor(2);
   gr_smallint->SetMarkerColor(3);


   TApplication *myapp;
   if(interactive)
      myapp=new TApplication("myapp",0,0);

   TMultiGraph* mg = new TMultiGraph();
   mg->Add(gr_rms);
   mg->Add(gr_fit);
   mg->Add(gr_smallint);
   TCanvas *cc = new TCanvas(); 
   mg->Draw("APL");
   cc->Print("RMS.pdf");

   data_amw.DrawHistos();
   data_amw.DrawProfiles();
   if(interactive)
      myapp->Run();
   //EvAnalyz data_rtcorr = data.RiseTimeCorrection();
   //EvAnalyz data_amw = data.MitigatedAmpCorrection(2000,5000);
   //data_amw.SetRiseTimeRange(-0.2,0.2);
   //EvAnalyz data_amw_poscorr = data_amw.PosCorrection();
   //data.DrawProfiles();
   //data_rtcorr.SetRiseTimeRange(-0.1,0.1);
   //data_rtcorr.DrawProfiles(-0.3,0.3);
   //data_amw.DrawProfiles(-0.3,0.3);
   //data_amw_poscorr.DrawProfiles(-0.3,0.3);
   return 0;
   
}
