#ifndef EVANALYZ_H
#define EVANALYZ_H

#include <iostream>
#include <string>
#include <map>

#include "TChain.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "ConfigFile.hh"
#include "TCanvas.h"
//#include "TH2.h"
//#include "TH2F.h"

//#include <fstream>
//#include <sstream>
//#include <utility>
//#include <vector>

//using std::string;
using namespace std;

class EvAnalyz 
{
   // Data
   protected:
      //ConfigFile fconfig;
      TChain* fDataTree;
      float fmu_y_hit, fmu_x_hit, fAMP_MAX;
      std::map<float,float> ftime;
      int fNthr;
      std::vector<float> fthr;
      std::string fDataLabel;
      float famp_min, famp_max;
      float frisetime_min, frisetime_max;
      float ftime_offset;
      std::map<float,TProfile*> fp_time_amp;
      std::map<float,TProfile*> fp_time_risetime;
      std::map<float,TProfile2D*> fp2_time_x_y;
      std::map<std::string,TCanvas*> fPlots;

   // Methods
   public:
      EvAnalyz(const ConfigFile & config);
      EvAnalyz(TChain* outtree, int Nthr, vector<float> thr, string DataLabel, float famp_min, float famp_max, float frisetime_min, float frisetime_max, float ftime_offset);
      ~EvAnalyz();
      void FillProfile(bool mkamp=true, bool mkrisetime=true, bool mkpos=true);
      EvAnalyz AmpCorrection();
      EvAnalyz MitigatedAmpCorrection(float amp_min_fit, float amp_max_fit);
      EvAnalyz PosCorrection();
      EvAnalyz RiseTimeCorrection();
      //void DrawTimeRes();
      void DrawProfiles(float time_min=0,float time_max=2);
      void SetAmpRange(float amp_min,float amp_max);
      void SetRiseTimeRange(float risetime_min,float risetime_max);
      //TChain& GetTree();
      //std::map<float,TProfile*>& Getp_time_amp();
      //std::map<float,TProfile*>& Getp_time_risetime();
      //std::map<float,TProfile2D*>& Getp2_time_x_y();

   protected:
      void SetBranchTree();
      void CreateProfile(bool mkamp=true, bool mkrisetime=true, bool mkpos=true);
      void ParseConfigFile(const ConfigFile & config);
};

#endif  // EVANALYZ_H
