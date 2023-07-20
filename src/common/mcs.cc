 /**
 * @file mcs.cc
 * @brief Class implementation for mcs handling
 */

#include "utils.h"
#include "mcs.h"
#include "modulation.h"



struct Modulation_Tables {
    Table<complex_float> ul_tables[4];
    Table<complex_float> dl_tables[4];
}
 

Mcs::Mcs (size_t initial_modulation) {
  current_modulation_.frame_number = 0;
  current_modulation_.modulation_type = initial_modulation;
  Modulation_Tables modulation_tables;
  Create_Modulation_Tables(modulation_tables);

}

Mcs::~Mcs()=default;

Mcs::Create_Modulation_Tables (Modulation_Tables modulation_tables){
  for (int i = 0; i < modulation_tables.dl_tables.size()){
      InitModulationTable(modulation_tables.dl_tables[i], (i+1)*2);
      InitModulationTable(modulation_tables.ul_tables[i], (i+1)*2);
  }
}

Mcs::Update_Modulation_Scheme(size_t current_frame_number){
  if (current_frame_number < next_modulation_.frame_number) {
    current_modulation_.frame_number = current_frame_number;
    current_modulation_.modulation_type = next_modulation_.modulation_type;
    // UPDATE THE MCS
  }
}

Mcs::Set_Next_Modulation_Scheme(Modulation_Scheme next_modulation_scheme) {
  next_modulation_.frame_number = next_modulation_scheme.frame_number;
  next_modulation_.modulation_type = next_modulation_scheme.modulation_type;
}