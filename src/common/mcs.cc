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

Mcs::Mcs(nlohmann::json ul_mcs_params_, nlohmann::json dl_mcs_params_, size_t ofdm_data_num) {
  
  Create_Modulation_Tables();

  initial_ul_mcs_properties_.base_graph = ul_mcs_params_.value("base_graph", 1);
  initial_ul_mcs_properties_.early_term = ul_mcs_params_.value("earlyTermination", true);
  initial_ul_mcs_properties_.max_decoder_iter = ul_mcs_params_.value("decoderIter", 5);
  inidiaul_ul_mcs_properties_.ofdm_data_num = ofdm_data_num;

  initial_dl_mcs_properties_.base_graph = dl_mcs_params_.value("base_graph", 1);
  initial_dl_mcs_properties_.early_term = dl_mcs_params_.value("earlyTermination", true);
  initial_dl_mcs_properties_.max_decoder_iter = dl_mcs_params_.value("decoderIter", 5);
  inidiaul_dl_mcs_properties_.ofdm_data_num = ofdm_data_num;


  //Initialize UL MCS
  Initialize_Ul_Mcs(ul_mcs_params_);

  //Initialize DL MCS
  Initialize_Dl_Mcs(dl_mcs_params_);

}

Mcs::~Mcs()=default;

void Mcs::Initialize_Ul_Mcs(const nlohmann::json& ul_mcs) {
  current_ul_mcs_.frame_number = 0;
  if (ul_mcs.find("mcs_index") == ul_mcs.end()) {
  std::string ul_modulation_type = ul_mcs.value("modulation", "16QAM");
  current_mcs_.modulation = kModulStringMap.at(ul_modulation_type);

  double ul_code_rate_usr = ul_mcs.value("code_rate", 0.333);
  size_t code_rate_int =
      static_cast<size_t>(std::round(ul_code_rate_usr * 1024.0));

  ul_mcs_index_ = CommsLib::GetMcsIndex(ul_mod_order_bits_, code_rate_int);
  ul_code_rate_ = GetCodeRate(ul_mcs_index_);
  if (ul_code_rate_ / 1024.0 != ul_code_rate_usr) {
    AGORA_LOG_WARN(
        "Rounded the user-defined uplink code rate to the closest standard "
        "rate %zu/1024.\n",
        ul_code_rate_);
  }
  } else {
    current_mcs_.mcs_index = ul_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_mcs_.modulation = GetModOrderBits(current_mcs_.mcs_index);
    current_mcs.code_rate = GetCodeRate(current_mcs_.mcs_index);
  }
}

void Mcs::Initialize_Dl_Mcs(const nlohmann::json& dl_mcs) {
  current_dl_mcs_.frame_number = 0;
  if (dl_mcs.find("mcs_index") == dl_mcs.end()) {
  dl_modulation_ = dl_mcs.value("modulation", "16QAM");
  dl_mod_order_bits_ = kModulStringMap.at(dl_modulation_);

  double dl_code_rate_usr = dl_mcs.value("code_rate", 0.333);
  size_t code_rate_int =
      static_cast<size_t>(std::round(dl_code_rate_usr * 1024.0));
  dl_mcs_index_ = CommsLib::GetMcsIndex(dl_mod_order_bits_, code_rate_int);
  dl_code_rate_ = GetCodeRate(dl_mcs_index_);
  if (dl_code_rate_ / 1024.0 != dl_code_rate_usr) {
    AGORA_LOG_WARN(
        "Rounded the user-defined downlink code rate to the closest standard "
        "rate %zu/1024.\n",
        dl_code_rate_);
  }
  } else {
    dl_mcs_index_ = dl_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    dl_mod_order_bits_ = GetModOrderBits(dl_mcs_index_);
    dl_modulation_ = MapModToStr(dl_mod_order_bits_);
    dl_code_rate_ = GetCodeRate(dl_mcs_index_);
    dl_modulation_ = MapModToStr(dl_mod_order_bits_);
  }
}

void Mcs::Create_Modulation_Tables() {
  for (int i = 0; i < modulation_tables_.dl_tables.size(); i++){
      InitModulationTable(modulation_tables_.dl_tables[i], (i+1)*2);
      InitModulationTable(modulation_tables_.ul_tables[i], (i+1)*2);
  }
}

void Mcs::Update_MCS_Schemes(size_t current_frame_number) {
  Update_MCS_ul_Scheme(size_t current_frame_number);
  Update_MCS_dl_Scheme(size_t current_frame_number);
}

void Mcs::Update_Ul_MCS_Scheme(size_t current_frame_number) {
  if (current_frame_number >= next_ul_mcs_.frame_number) {
    current_ul_mcs_.frame_number = current_frame_number;
    current_ul_mcs_.mcs_index = next_ul_mcs_.mcs_index;
  }
  ul_ldpc_config_ = Update_Ul_Ldpc_Config();
}

void Mcs::Update_Dl_MCS_Scheme(size_t current_frame_number) {
  if (current_frame_number >= next_dl_mcs_.frame_number) {
    current_dl_mcs_.frame_number = current_frame_number;
    current_dl_mcs_.mcs_index = next_dl_mcs_.mcs_index;
  }
  dl_ldpc_config_ = Update_Dl_Ldpc_Config();
}

void Mcs::Set_Next_Ul_MCS_Scheme(MCS_Scheme next_mcs_scheme) {
  next_ul_mcs_.frame_number = next_mcs_scheme.frame_number;
  next_ul_mcs_.mcs_index = next_mcs_Scheme.modulation_type;
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
                       inidiaul_dl_mcs_properties_.ofdm_data_num, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows =
      static_cast<size_t>(
          std::round(1024.0 * LdpcNumInputCols(base_graph) / ul_code_rate)) -
      (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  ul_ldpc_config_ = LDPCconfig(base_graph, zc, max_decoder_iter, early_term,
                               num_cb_len, num_cb_codew_len, num_rows, 0);

  ul_ldpc_config_.NumBlocksInSymbol((inidiaul_dl_mcs_properties_.ofdm_data_num* ul_mod_order_bits) /
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
                       inidiaul_dl_mcs_properties_.ofdm_data_num, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows =
      static_cast<size_t>(
          std::round(1024.0 * LdpcNumInputCols(base_graph) / dl_code_rate)) -
      (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  dl_ldpc_config_ = LDPCconfig(base_graph, zc, max_decoder_iter, early_term,
                               num_cb_len, num_cb_codew_len, num_rows, 0);

  dl_ldpc_config_.NumBlocksInSymbol((inidiaul_dl_mcs_properties_.ofdm_data_num * dl_mod_order_bits) /
                                    ul_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumULSyms() == 0) || (ul_ldpc_config_.NumBlocksInSymbol() > 0),
      "Uplink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");

}

LDPCconfig Mcs::Ul_Ldpc_Config() {
  return ul_ldpc_config_;
}

LDPCconfig Mcs::Dl_Ldpc_Config() {
  return dl_ldpc_config_;
}

