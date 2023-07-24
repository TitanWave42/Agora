 /**
 * @file mcs.cc
 * @brief Class implementation for mcs handling
 */

#include "utils.h"
#include "mcs.h"
#include "modulation.h"
#include "config.h"
#include "ldpc_config.h"
#include "comms-lib.h"
#include "utils_ldpc.h"
#include "logger.h"
#include <stddef.h>

Mcs::Mcs(nlohmann::json ul_mcs_params, nlohmann::json dl_mcs_params, size_t ofdm_data_num) {
  
  ul_mcs_ = ul_mcs_params;
  dl_mcs_ = dl_mcs_params;

  Create_Modulation_Tables();

  initial_ul_mcs_properties_.base_graph = ul_mcs_params.value("base_graph", 1);
  initial_ul_mcs_properties_.early_term = ul_mcs_params.value("earlyTermination", true);
  initial_ul_mcs_properties_.max_decoder_iter = ul_mcs_params.value("decoderIter", 5);
  initial_ul_mcs_properties_.ofdm_data_num = ofdm_data_num;

  initial_dl_mcs_properties_.base_graph = dl_mcs_params.value("base_graph", 1);
  initial_dl_mcs_properties_.early_term = dl_mcs_params.value("earlyTermination", true);
  initial_dl_mcs_properties_.max_decoder_iter = dl_mcs_params.value("decoderIter", 5);
  initial_dl_mcs_properties_.ofdm_data_num = ofdm_data_num;


  //Initialize UL MCS
  Initialize_Ul_Mcs(ul_mcs_);

  //Initialize DL MCS
  Initialize_Dl_Mcs(dl_mcs_);

}

Mcs::~Mcs()=default;

void Mcs::Initialize_Ul_Mcs(const nlohmann::json ul_mcs) {
  current_ul_mcs_.frame_number = 0;
  if (ul_mcs.find("mcs_index") == ul_mcs.end()) {
  std::string ul_modulation_type = ul_mcs.value("modulation", "16QAM");
  current_ul_mcs_.modulation = kModulStringMap.at(ul_modulation_type);

  double ul_code_rate_usr = ul_mcs.value("code_rate", 0.333);
  size_t code_rate_int =
      static_cast<size_t>(std::round(ul_code_rate_usr * 1024.0));

  current_ul_mcs_.mcs_index = CommsLib::GetMcsIndex(current_ul_mcs_.modulation, code_rate_int);
  current_ul_mcs_.code_rate = GetCodeRate(current_ul_mcs_.mcs_index);
  if (current_ul_mcs_.code_rate / 1024.0 != ul_code_rate_usr) {
    AGORA_LOG_WARN(
        "Rounded the user-defined uplink code rate to the closest standard "
        "rate %zu/1024.\n",
        current_ul_mcs_.code_rate);
  }
  } else {
    current_ul_mcs_.mcs_index = ul_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_ul_mcs_.modulation = GetModOrderBits(current_ul_mcs_.mcs_index);
    current_ul_mcs_.code_rate = GetCodeRate(current_ul_mcs_.mcs_index);
  }
}

void Mcs::Initialize_Dl_Mcs(const nlohmann::json dl_mcs) {
  current_dl_mcs_.frame_number = 0;
  if (dl_mcs.find("mcs_index") == dl_mcs.end()) {
  std::string dl_modulation_type = dl_mcs.value("modulation", "16QAM");
  current_dl_mcs_.modulation = kModulStringMap.at(dl_modulation_type);

  double dl_code_rate_usr = dl_mcs.value("code_rate", 0.333);
  size_t code_rate_int =
      static_cast<size_t>(std::round(dl_code_rate_usr * 1024.0));
  current_dl_mcs_.mcs_index = CommsLib::GetMcsIndex(current_dl_mcs_.modulation, code_rate_int);
  current_dl_mcs_.code_rate = GetCodeRate(current_dl_mcs_.mcs_index);
  if (current_dl_mcs_.code_rate / 1024.0 != dl_code_rate_usr) {
    AGORA_LOG_WARN(
        "Rounded the user-defined downlink code rate to the closest standard "
        "rate %zu/1024.\n",
        current_dl_mcs_.code_rate);
  }
  } else {
    current_dl_mcs_.mcs_index = dl_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_dl_mcs_.modulation = GetModOrderBits(current_dl_mcs_.mcs_index);
    current_dl_mcs_.code_rate = GetCodeRate(current_dl_mcs_.mcs_index);
  }
}

void Mcs::Create_Modulation_Tables() {
  for (int i = 0; i < modulation_tables_.dl_tables.size(); i++){
      InitModulationTable(modulation_tables_.dl_tables[i], (i+1)*2);
      InitModulationTable(modulation_tables_.ul_tables[i], (i+1)*2);
  }
}

void Mcs::Update_MCS_Schemes(size_t current_frame_number) {
  Update_Ul_MCS_Scheme(size_t current_frame_number);
  Update_Dl_MCS_Scheme(size_t current_frame_number);
}

void Mcs::Update_Ul_MCS_Scheme(size_t current_frame_number) {
  if (current_frame_number >= next_ul_mcs_.frame_number) {
    current_ul_mcs_.frame_number = current_frame_number;
    current_ul_mcs_.mcs_index = next_ul_mcs_.mcs_index;
  }
  Update_Ul_Ldpc_Config();
}

void Mcs::Update_Dl_MCS_Scheme(size_t current_frame_number) {
  if (current_frame_number >= next_dl_mcs_.frame_number) {
    current_dl_mcs_.frame_number = current_frame_number;
    current_dl_mcs_.mcs_index = next_dl_mcs_.mcs_index;
  }
  Update_Dl_Ldpc_Config();
}

void Mcs::Set_Next_Ul_MCS_Scheme(MCS_Scheme next_mcs_scheme) {
  next_ul_mcs_.frame_number = next_mcs_scheme.frame_number;
  next_ul_mcs_.mcs_index = next_mcs_scheme.modulation;
}

void Mcs::Set_Next_Dl_MCS_Scheme(MCS_Scheme next_mcs_scheme) {
  next_dl_mcs_.frame_number = next_mcs_scheme.frame_number;
  next_dl_mcs_.mcs_index = next_mcs_scheme.mcs_index;
}


void Mcs::Update_Ul_Ldpc_Config() {

  size_t ul_mod_order_bits = GetModOrderBits(current_ul_mcs_.mcs_index);
  size_t ul_code_rate = GetCodeRate(current_ul_mcs_.mcs_index);

  Table<complex_float> ul_mod_table_ = modulation_tables_.ul_tables[ul_mod_order_bits / 2 - 1];

  size_t zc = SelectZc(initial_ul_mcs_properties_.base_graph, ul_code_rate, ul_mod_order_bits,
                       initial_dl_mcs_properties_.ofdm_data_num, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows =
      static_cast<size_t>(
          std::round(1024.0 * LdpcNumInputCols(base_graph) / ul_code_rate)) -
      (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  *ul_ldpc_config_ = LDPCconfig(base_graph, zc, max_decoder_iter, early_term,
                               num_cb_len, num_cb_codew_len, num_rows, 0);

  ul_ldpc_config_.NumBlocksInSymbol((initial_dl_mcs_properties_.ofdm_data_num* ul_mod_order_bits) /
                                    ul_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumULSyms() == 0) || (ul_ldpc_config_.NumBlocksInSymbol() > 0),
      "Uplink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");
}

void Mcs::Update_Dl_Ldpc_Config(){

  size_t dl_mod_order_bits = GetModOrderBits(current_dl_mcs_.mcs_index);
  size_t dl_code_rate = GetCodeRate(current_dl_mcs_.mcs_index);

  Table<complex_float> dl_mod_table_ = modulation_tables_.dl_tables[dl_mod_order_bits / 2 - 1];

  size_t zc = SelectZc(initial_ul_mcs_properties_.base_graph, dl_code_rate, dl_mod_order_bits,
                       initial_dl_mcs_properties_.ofdm_data_num, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows =
      static_cast<size_t>(
          std::round(1024.0 * LdpcNumInputCols(base_graph) / dl_code_rate)) -
      (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  *dl_ldpc_config_ = LDPCconfig(base_graph, zc, max_decoder_iter, early_term,
                               num_cb_len, num_cb_codew_len, num_rows, 0);

  dl_ldpc_config_.NumBlocksInSymbol((initial_dl_mcs_properties_.ofdm_data_num * dl_mod_order_bits) /
                                    ul_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumULSyms() == 0) || (ul_ldpc_config_.NumBlocksInSymbol() > 0),
      "Uplink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");

}

// LDPCconfig Mcs::get_Ul_Ldpc_Config() {
//   return ul_ldpc_config_;
// }

// LDPCconfig Mcs::get_Dl_Ldpc_Config() {
//   return dl_ldpc_config_;
// }

