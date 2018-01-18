#include <iostream>
#include "ConfigFile.hh"
#include "EvAnalyz.hh"

int main (int argc, char **argv)
{
   if(argc!=2)
   {
      cout<<"ERROR: unvalid number of input parameters\n";
      exit(EXIT_FAILURE);
   }
   ConfigFile config(argv[1]);
   EvAnalyz data(config);
   EvAnalyz data_rtcorr = data.RiseTimeCorrection();
   //EvAnalyz data_amw = data.MitigatedAmpCorrection(2000,5000);
   //data_amw.SetRiseTimeRange(-0.2,0.2);
   //EvAnalyz data_amw_poscorr = data_amw.PosCorrection();
   data.DrawProfiles();
   data_rtcorr.SetRiseTimeRange(-0.1,0.1);
   data_rtcorr.DrawProfiles(-0.3,0.3);
   //data_amw.DrawProfiles(-0.3,0.3);
   //data_amw_poscorr.DrawProfiles(-0.3,0.3);
   
}
