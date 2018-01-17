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
   data.DrawProfiles();
}
