
/*
 *  The file '/usr/share/fasterac/data/smallgrp.fast' is a faster data file
 *  containing four channels of qdc_tdc data.
 *	The channels 1 and 3 are systematically
 *  grouped in this file (coinc).
 *  This program reads the data file, put the grouped channels 1 & 3 in a root
 *tree and saves the given tree in 'group2tree_data.root'.
 *
 *
 *
 */

//  std includes
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//  root includes
#include "TFile.h"
#include "TTree.h"

//  fasterac includes
#include "fasterac/adc.h"
#include "fasterac/fast_data.h"
#include "fasterac/fasterac.h"
#include "fasterac/group.h"

using namespace std;

#define DATAFILENAME "data/calib_Si_5_0001.fast"
#define ROOTFILENAME "root/calib_Si_5.root"

int main(int argc, char **argv) {

  /************/
  /*  FASTER  */
  /************/
  //  file reader
  faster_file_reader_p reader;
  //  data
  faster_data_p data;
  unsigned char alias;
  unsigned short label;
  unsigned long long clock;
  //  group data
  faster_buffer_reader_p group_reader;
  void *group_buffer;
  unsigned short lsize;
  // faster_data_p          group_data;
  //   qdc tdc data  (from faster group)
  adc_data adc2, adc4;

  /**********/
  /*  ROOT  */
  /**********/
  //  root leaves : channels 1
  // Int_t                  leaf_q1;
  Double_t leaf_Si_1, leaf_TAC_1, leaf_Si_2, leaf_TAC_2, leaf_Si_B, leaf_TAC_B,
      leaf_Ge, leaf_faster_time_1, leaf_faster_time_2, leaf_faster_time_B,
      clock_Si_1, clock_tac_1, clock_Si_2, clock_tac_2, clock_Si_B, clock_tac_B,
      clock_Ge;
  //  root tree
  TTree *tree;
  char tree_title[256];
  //  root file
  TString fName = ROOTFILENAME;
  TFile *root_file;

  //  print infos
  printf("\n");
  printf("  group2tree :\n");
  printf("     - read faster file '%s'\n", DATAFILENAME);
  printf("     - get grouped data (labels %d and %d)\n", 3, 7767);
  printf("     - output those data to root file '%s'\n", ROOTFILENAME);
  printf("\n");

  //  open faster file reader
  reader = faster_file_reader_open(DATAFILENAME);
  if (reader == NULL) {
    printf("error opening file %s\n", DATAFILENAME);
    return EXIT_FAILURE;
  }

  //  output root file
  root_file = new TFile(fName.Data(), "recreate");
  tree = new TTree("DataTree", tree_title);
  tree->Branch("Si_1", &leaf_Si_1, "Si_1/D");
  tree->Branch("TAC_1", &leaf_TAC_1, "TAC_1/D");
  tree->Branch("faster_time_1", &leaf_faster_time_1, "faster_time_1/D");
  tree->Branch("Si_2", &leaf_Si_2, "Si_2/D");
  tree->Branch("TAC_2", &leaf_TAC_2, "TAC_2/D");
  tree->Branch("faster_time_2", &leaf_faster_time_2, "faster_time_2/D");
  tree->Branch("Si_B", &leaf_Si_B, "Si_B/D");
  tree->Branch("TAC_B", &leaf_TAC_B, "TAC_B/D");
  tree->Branch("faster_time_B", &leaf_faster_time_B, "faster_time_B/D");
  tree->Branch("Ge", &leaf_Ge, "Ge/D");

  // main loop
  while ((data = faster_file_reader_next(reader)) != NULL) { //  read each data
    alias = faster_data_type_alias(data);
    lsize = faster_data_load_size(data);
    label = faster_data_label(data);
    clock = faster_data_clock_ns(data);

    //    cout << alias<<"  "<<CRRC4_SPECTRO_TYPE_ALIAS << endl;

    // reading Si first
    if (alias == CRRC4_SPECTRO_TYPE_ALIAS) {
      leaf_Si_1 = 0;
      leaf_TAC_1 = 0;
      leaf_faster_time_1 = 0;
      leaf_Si_2 = 0;
      leaf_TAC_2 = 0;
      leaf_faster_time_2 = 0;
      leaf_Si_B = 0;
      leaf_TAC_B = 0;

      if (label == 3) { //
        faster_data_load(data, &adc2);
        leaf_Si_1 = adc2.measure;
        clock_Si_1 = faster_data_clock_ns(data);
        // next event should be TAC
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_1 = faster_data_clock_ns(data);
          // std::cout << "diff = " << clock_tac_1 - clock_Si_1 << "\n";
          if (label == 4 && clock_tac_1 - clock_Si_1 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_faster_time_1 = clock_tac_1 - clock_Si_1;
            leaf_TAC_1 = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_2 = faster_data_clock_ns(data);
          if (label == 4 && clock_tac_2 - clock_Si_2 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_faster_time_1 = clock_tac_1 - clock_Si_1;
            leaf_TAC_1 = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_2 = faster_data_clock_ns(data);
          if (label == 4 && clock_tac_2 - clock_Si_2 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_faster_time_1 = clock_tac_1 - clock_Si_1;
            leaf_TAC_1 = adc4.measure;
          }
        } else {
          leaf_TAC_1 = 5000.;
        }
        tree->Fill();
      }

      if (label == 5) { //
        faster_data_load(data, &adc2);
        leaf_Si_2 = adc2.measure;
        clock_Si_2 = faster_data_clock_ns(data);
        // next event should be TAC
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_2 = faster_data_clock_ns(data);
          if (label == 6 && clock_tac_2 - clock_Si_2 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_2 = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_2 = faster_data_clock_ns(data);
          if (label == 6 && clock_tac_2 - clock_Si_2 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_2 = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_2 = faster_data_clock_ns(data);
          if (label == 6 && clock_tac_2 - clock_Si_2 < 1000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_2 = adc4.measure;
          }
        } else {
          leaf_TAC_2 = 5000.;
        }
        tree->Fill();
      }

      if (label == 7) { //
        faster_data_load(data, &adc2);
        leaf_Si_B = adc2.measure;
        clock_Si_B = faster_data_clock_ns(data);
        // next event should be TAC
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_B = faster_data_clock_ns(data);
          // std::cout << "diff = " << clock_tac_B - clock_Si_B << "\n";
          if (label == 8 && clock_tac_B - clock_Si_B < 10000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_B = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_B = faster_data_clock_ns(data);
          if (label == 8 && clock_tac_B - clock_Si_B < 10000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_B = adc4.measure;
          }
        }
        if ((data = faster_file_reader_next(reader)) !=
            NULL) { //  read each data
          label = faster_data_label(data);
          clock_tac_B = faster_data_clock_ns(data);
          if (label == 8 && clock_tac_B - clock_Si_B < 10000) { //
            faster_data_load(data, &adc4);
            leaf_TAC_B = adc4.measure;
          }
        } else {
          leaf_TAC_B = 5000.;
        }
        tree->Fill();
      }

      if (label == 9) { //
        faster_data_load(data, &adc2);
        leaf_Ge = adc2.measure;
        clock_Ge = faster_data_clock_ns(data);
        tree->Fill();
      }

      // for time calibration//
      // if (label == 4) {
      //   faster_data_load(data, &adc4);
      //   leaf_TAC_1 = adc4.measure;
      //   //  next event should be TAC
      //   if ((data = faster_file_reader_next(reader)) !=
      //       NULL) { //  read each data
      //     label = faster_data_label(data);
      //     faster_data_load(data, &adc4);
      //     leaf_TAC_1 = adc4.measure;
      //     tree->Fill();
      //   } else
      //     break;
      // }

      // if (label == 6) {
      //   faster_data_load(data, &adc4);
      //   leaf_TAC_2 = adc4.measure;
      //   //  next event should be TAC
      //   if ((data = faster_file_reader_next(reader)) !=
      //       NULL) { //  read each data
      //     label = faster_data_label(data);
      //     faster_data_load(data, &adc4);
      //     leaf_TAC_2 = adc4.measure;
      //     tree->Fill();
      //   } else
      //     break;
      // }

      // if (label == 8) {
      //   faster_data_load(data, &adc4);
      //   leaf_TAC_B = adc4.measure;
      //   //  next event should be TAC
      //   if ((data = faster_file_reader_next(reader)) !=
      //       NULL) { //  read each data
      //     label = faster_data_label(data);
      //     faster_data_load(data, &adc4);
      //     leaf_TAC_B = adc4.measure;
      //     tree->Fill();
      //   } else
      //     break;
      // }
    }
  }

  //  close files & quit
  faster_file_reader_close(reader);
  root_file->Write();
  root_file->Close();
  return EXIT_SUCCESS;
}
