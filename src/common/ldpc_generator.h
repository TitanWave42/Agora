/**
 * @file ldpc_update.h
 * @brief Class definition for ldpc_update.
 */

#ifndef LDPC_GENERATOR_H_
#define LDPC_GENERATOR_H_

#include "ldpc_config.h"

static constexpr size_t NumMcsIndices = 32;

struct McsParams {
  LDPCconfig ldpc_config_(0, 0, 0, false, 0, 0, 0, 0);
  bool early_term_;
  std::string modulation_ = NULL;
  size_t mod_order_bits_;
  size_t base_graph_;
  size_t max_decoder_iter_;
  size_t mcs_index;
  size_t num_bytes_per_cb_;
  size_t num_padding_bytes_per_cb_;
  size_t data_bytes_num_persymbol_;
  size_t mac_packet_length_;
  size_t mac_packets_perframe_;
  size_t mac_data_bytes_num_perframe_;
  size_t mac_data_length_max_;
  double code_rate_usr_;
};

//Not sure where else to put this.
struct OfdmConfig {
  size_t ofdm_data_num_ul_;
  size_t ofdm_data_num_dl_;
  size_t ofdm_ctrl_data_;
}

static constexpr size_t kCbPerSymbol = 1;

class LdpcGenerator {
 public:
  LdpcGenerator(const OfdmConfig ofdm_data, uint16_t ul_base_graph,
                bool ul_early_term, int16_t ul_max_decoder_iter,
                uint16_t dl_base_graph, bool dl_early_term,
                int16_t dl_max_decoder_iter);
  ~LdpcGenerator();

  inline std::vector < pair<McsParams*, McsParams*> LdpcTable() {
    return this->ul_dl_ldpc_table_;
  }

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

  std::vector < pair<McsParams*, McsParams*> ul_dl_ldpc_table_;

  LDPCconfig UpdateUlLdpcConfig(size_t ul_mcs_index) LDPCconfig
      UpdateDlLdpcConfig(size_t dl_mcs_index);
  LDPCconfig UpdateCtrlMCS(size_t control_mcs_index);

  void GenerateLdpcTable();
  void CalculateUlLdpcProperties(McsParams* ul_mcs_params);
  void CalculateDlLdpcProperties(McsParams* dl_mcs_params);
  inline size_t SelectZc(size_t base_graph, size_t code_rate,
                         size_t mod_order_bits, size_t num_sc,
                         size_t cb_per_sym, const std::string& dir);
}

#endif /* LDPC_GENERATOR */