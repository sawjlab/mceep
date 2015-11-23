//
// read the ntuple from stdin 
//

#include <iostream.h>
#include <fstream.h>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"

// this file is created by prep_m2root.py
#include "m2root.h"

int main()
{
  init_var(); // initialize variable names

  // create output file
  TFile *f = new TFile("mceep.root", "RECREATE", "Mceep results");
  
  Text_t header[80];
  
  Int_t nwt=0;
  Int_t nvar = 0;

  Int_t ivar[100];
  Int_t i;

  Float_t val[100];
  //
  //
  cin >> header;
  cin >> nwt >> nvar ;

  cout << "header = " << header << endl;
  cout << "nwt = " << nwt << endl ;
  cout << "nvar = " << nvar << endl ;

  // read the idices of the variables in this n-tuple
  for (i=0; i<nvar; i++) cin >> ivar[i];

  // create the variable name string for the ntuple

  TString var_string("WEIGHT:"); //
  Int_t n;

  // add the remaining variable names

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
  while(! cin.eof() ){
    n++;
    for (i=0; i<nvar+1; i++) { // nvar+1 elements because of the weight
      cin >> val[i];            // which is stored in val[0]
    }
    // fill the ntuple with it
    mceep->Fill(val);
  }
  cout << "Read " << n << " events" << endl;

  f->Write(); // write output file
  return(0);
}
