/**
 * @file ldpc_update.h
 * @brief Class definition for ldpc_update.
 */

#ifndef LDPC_UPDATE_H_
#define LDPC_UPDATE_H_

#include "ldpc_config.h"

static constexpr size_t NumMcsIndices = 32;


struct UlMcsParams {
  LDPCconfig ul_ldpc_config_;
  bool ul_early_term_;
  std::string ul_modulation_ = NULL;
  size_t ul_mod_order_bits_;
  size_t ul_base_graph_;
  size_t ul_max_decoder_iter_;
  size_t ul_mcs_index;
  size_t ul_num_bytes_per_cb_;
  size_t ul_num_padding_bytes_per_cb_;
  size_t ul_data_bytes_num_persymbol_;
  size_t ul_mac_packet_length_;
  size_t ul_mac_packets_perframe_;
  size_t ul_mac_data_bytes_num_perframe_;
  size_t ul_mac_data_length_max_;
  double ul_code_rate_usr_;
};

struct DlMcsParams {
  LDPCconfig dl_ldpc_config_;
  bool dl_early_term_;
  std::string dl_modulation_ = NULL;
  size_t dl_mod_order_bits_;
  size_t dl_base_graph_;
  size_t dl_max_decoder_iter_;
  size_t dl_mcs_index;
  size_t dl_num_bytes_per_cb_;
  size_t dl_num_padding_bytes_per_cb_;
  size_t dl_data_bytes_num_persymbol_;
  size_t dl_mac_packet_length_;
  size_t dl_mac_packet_length_;
  size_t dl_mac_packets_perframe_;
  size_t dl_mac_data_bytes_num_perframe_;
  size_t dl_mac_data_length_max_;
  double dl_code_rate_usr_;
};

// struct InitialMcsProperties {
//   uint16_t base_graph_;
//   bool early_term_;
//   int16_t max_decoder_iter_;
// };

//Not sure where else to put this.
struct OfdmConfig {
  size_t ofdm_data_num_ul_;
  size_t ofdm_data_num_dl_;
  size_t ofdm_ctrl_data_;
}

class LdpcUpdater {
 public:
  LdpcUpdater(const OfdmConfig ofdm_data, uint16_t ul_base_graph, bool ul_early_term,
                         int16_t ul_max_decoder_iter, uint16_t dl_base_graph,
                         bool dl_early_term, int16_t dl_max_decoder_iter);
  ~LdpcUpdater();

 private:
  const size_t ofdm_data_num_ul_;
  const size_t ofdm_data_num_dl_;
  const size_t ofdm_ctrl_data_;
  const uint16_t ul_base_graph_;
  const bool ul_early_term_;
  const int16_t ul_max_decoder_iter_; 
  const uint16_t dl_base_graph_; 
  const bool dl_early_term_; 
  const int16_t dl_max_decoder_iter_;
  McsScheme next_ul_mcs_;
  McsScheme next_dl_mcs_;
  McsScheme current_ul_mcs_;
  McsScheme current_dl_mcs_;
  

  std::vector<pair<UlMcsParams*, DlMcsParams*> ul_dl_ldpc_table_;

  LDPCconfig UpdateUlLdpcConfig(size_t ul_mcs_index) 
  LDPCconfig UpdateDlLdpcConfig(size_t dl_mcs_index);
  LDPCconfig UpdateCtrlMCS(size_t control_mcs_index);

  void GenerateLdpcTable();
  void CalculateUlLdpcProperties(UlMcsParams* ul_mcs_params);
  void CalculateDlLdpcProperties(DlMcsParams* dl_mcs_params);
  
}

#endif /* LDPC_UPDATE_H_ */