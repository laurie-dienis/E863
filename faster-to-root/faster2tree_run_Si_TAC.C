
/*
 *  The file '/usr/share/fasterac/data/smallgrp.fast' is a faster data file
 *  containing four channels of qdc_tdc data. 
 *	The channels 1 and 3 are systematically
 *  grouped in this file (coinc).
 *  This program reads the data file, put the grouped channels 1 & 3 in a root tree
 *  and saves the given tree in 'group2tree_data.root'.
 *
 *  
 *
 */


//  std includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//  root includes
#include "TFile.h"
#include "TTree.h"

//  fasterac includes
#include "fasterac/fasterac.h"
#include "fasterac/fast_data.h"
#include "fasterac/adc.h"
#include "fasterac/group.h"

using namespace std;


#define DATAFILENAME "Si_TAC_0001.fast"
#define ROOTFILENAME "run.root"




int main (int argc, char** argv) {



  /************/
  /*  FASTER  */
  /************/
  //  file reader
  faster_file_reader_p   reader;
  //  data
  faster_data_p          data;
  unsigned char          alias;
  unsigned short         label;
  unsigned long long      clock;
  //  group data
  faster_buffer_reader_p group_reader;
  void*                  group_buffer;
  unsigned short         lsize;
  //faster_data_p          group_data;
  //  qdc tdc data  (from faster group)
  adc_data	adc2,adc4;
  
  
  
  /**********/
  /*  ROOT  */
  /**********/
  //  root leaves : channels 1
 // Int_t                  leaf_q1;
  Double_t        leaf_Si, leaf_TAC;
  //  root tree
  TTree*                 tree;
  char                   tree_title [256];
  //  root file
  TString                fName = ROOTFILENAME;
  TFile*                 root_file;


  //  print infos
  printf ("\n");
  printf ("  group2tree :\n");
  printf ("     - read faster file '%s'\n", DATAFILENAME);
  printf ("     - get grouped data (labels %d and %d)\n", 3,7767);
  printf ("     - output those data to root file '%s'\n", ROOTFILENAME);
  printf ("\n");




  //  open faster file reader
  reader = faster_file_reader_open (DATAFILENAME);
  if (reader == NULL) {
    printf ("error opening file %s\n", DATAFILENAME);
    return EXIT_FAILURE;
  }
  
  
  
  
 
  //  output root file
  root_file = new TFile (fName.Data (), "recreate");
  tree = new TTree ("DataTree", tree_title);
  tree->Branch ("Si", &leaf_Si, "Si/D");
  tree->Branch ("TAC", &leaf_TAC, "TAC/D");




  // main loop
  while ((data = faster_file_reader_next (reader)) != NULL) {                    //  read each data
    alias = faster_data_type_alias(data);
    lsize = faster_data_load_size(data);
    label = faster_data_label(data);
    clock = faster_data_clock_ns(data);
    
//    cout << alias<<"  "<<CRRC4_SPECTRO_TYPE_ALIAS << endl;
    
    
   // reading Si first 
   if (alias == CRRC4_SPECTRO_TYPE_ALIAS) {
        if (label == 2) {                                                   //  
          faster_data_load (data, &adc2);
          leaf_Si = adc2.measure;
          // next event should be TAC
          data = faster_file_reader_next (reader);      
          label = faster_data_label(data);
        	if (label == 4) {                                                   //  
          	faster_data_load (data, &adc4);
          	leaf_TAC = adc4.measure; 
          	    }    
          	else {
          	leaf_TAC = 5000.;
          	} 
		tree->Fill ();
        } 
      }
      
      

               

      
    }

  //  close files & quit
  faster_file_reader_close (reader);
  root_file->Write ();
  root_file->Close ();
  return EXIT_SUCCESS;
}


