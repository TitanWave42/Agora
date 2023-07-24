/**
 * @file mcs.h
 * @brief Class definition for mcs.
 */
#ifndef MCS_H_
#define MCS_H_

#include "config.h"
#include "utils.h"
#include "nlohmann/json.hpp"

const size_t kNumTables = 5;

struct MCS_Scheme {
  size_t frame_number;
  size_t mcs_index;
  size_t modulation;
  size_t code_rate;
};

struct Modulation_Tables {
  Table<complex_float> ul_tables[kNumTables];
  Table<complex_float> dl_tables[kNumTables];
};

//I'm adding this as a struct here
//so that instead of reading from the conf file
// every time we want to change the modulation
// we can read once at the begining and access these
//variables whenever needed.

struct Initial_Mcs_Properties {
  uint16_t base_graph;
  bool early_term;
  int16_t max_decoder_iter;
  size_t ofdm_data_num;
};

class Mcs {

  public:
    explicit Mcs(nlohmann::json ul_mcs_params, nlohmann::json dl_mcs_params, size_t ofdm_data_num);
    ~Mcs();


    // LDPCconfig get_Ul_Ldpc_Config();
    // LDPCconfig get_Dl_Ldpc_Config();

    void Create_Modulation_Tables();
    void Update_MCS_Schemes(size_t current_frame_number);

    void Set_Next_Ul_MCS_Scheme(MCS_Scheme next_mcs_scheme);
    void Set_Next_Dl_MCS_Scheme(MCS_Scheme next_mcs_scheme);

    void Update_Ul_Ldpc_Config();
    void Update_Dl_Ldpc_Config();
    
    
  private:
    nlohmann::json ul_mcs_;
    nlohmann::json dl_mcs_;
    MCS_Scheme current_ul_mcs_;
    MCS_Scheme current_dl_mcs_;

    MCS_Scheme next_ul_mcs_;
    MCS_Scheme next_dl_mcs_;

    Modulation_Tables modulation_tables_;
    LDPCconfig* ul_ldpc_config_; //Uplink LDPC parameters
    LDPCconfig* dl_ldpc_config_; //Downlink LDPC parameters

    Initial_Mcs_Properties initial_ul_mcs_properties_;
    Initial_Mcs_Properties initial_dl_mcs_properties_;
    
    void Initialize_Ul_Mcs(const  nlohmann::json ul_mcs);
    void Initialize_Dl_Mcs(const  nlohmann::json dl_mcs);

    void Update_Ul_MCS_Scheme(size_t current_frame_number);
    void Update_Dl_MCS_Scheme(size_t current_frame_number);


};

#endif  // MCS_H_
