//
// macro to read a mceep kinemtaics ntuple into a root ntuple
//

// need the Mceep.C macro with the valiables defined

#include <iostream.h>
#include <fstream.h>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"

#include "m2root.h"

int main()
{
  init_var(); // initialize variable names

  TString KINEM("KINEM"); // header info for kinematics n-tuple
  TString TRANS("TRANS"); // header info for transport n-tuple

  // create output file
  TFile *f = new TFile("mceep.root", "RECREATE", "Mceep results");
  
  // open input file
  ifstream mc("mceep.ntu");

  Text_t header[80];
  
  Int_t nwt=0;
  Int_t nvar = 0;

  Int_t ivar[100];
  Int_t i;

  Float_t val[100];
  //
  //
  mc >> header;
  mc >> nwt >> nvar ;

  Bool_t is_kin   = KINEM.Contains(Strip(header));
  Bool_t is_trans = TRANS.Contains(Strip(header));

  if (is_trans){
    cout << " Cannot handle tranport N-tuples at this time . Sorry ! " << endl;
    mc.close();
    f->Close();
    return (-1);
  }
  // print the header information

  cout << "Ntuple type                 = " << header << endl;
  cout << "Weight index                = " << nwt << endl ;
  cout << "Number of Ntuple variables = " << nvar << endl ;

  // check the header information
  
  // read the idices of the variables in this n-tuple
  for (i=0; i<nvar; i++) mc >> ivar[i];

  // create the variable name string for the ntuple

  Int_t n;
  TString var_string; //

  // add the remaining variable names
  
  // the weight variables
  if (nwt < 0){
    var_string.Append("NUMER:");
    var_string.Append("DENOM:");
    nwt = 2;                       // used later
  }
  else{
    var_string.Append("WEIGHT:");
    nwt = 1;
  }

  // the ntuple variable 
  for (i=0; i<nvar-1;i++){
    n = ivar[i];
    var_string.Append(var[n].Strip()); // get rid of blanks
    var_string.Append(":");
  }
  n = ivar[i]; // the last variable does not need a colon
  var_string.Append(var[n].Strip());

  cout << "variable string = " << var_string.Data() << endl;

  // create the N-tuple
  
  TNtuple *mceep = new TNtuple("mceep", "Mceep N-Tuple", var_string.Data());

  n = 0;
  while(! mc.eof() ){
    n++;
    for (i=0; i<nvar+nwt; i++) { // nvar+1 elements because of the weight
      mc >> val[i];              // which is stored in val[0]
                                 // nwt has been adjusted above
    }
    // fill the ntuple with it
    mceep->Fill(val);
  }
  cout << "Read " << n << " events" << endl;

  mc.close();

  f->Write(); // write output file
  return(0);
}
