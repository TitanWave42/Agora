/**
 * @file ldpc_update.h
 * @brief Class definition for ldpc_update.
 */

#ifndef LDPC_UPDATE_H_
#define LDPC_UPDATE_H_

#include "ldpc_config.h"

struct UlMcsParams {
  LDPCconfig ul_ldpc_config_;
  bool ul_early_term_;
  size_t ul_modulation_;
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
  size_t dl_modulation_;
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

class LdpcUpdater {
 public:
  LdpcUpdater(uint16_t ul_base_graph, bool ul_early_term,
              int16_t ul_max_decoder_iter, uint16_t dl_base_graph,
              bool dl_early_term, int16_t dl_max_decoder_iter);
  ~LdpcUpdater();

  UpdateLdpc(LDPCconfig ldpc_config, Direction dir, size_t mcs_index);

 private:
  const uint16_t ul_base_graph_ const bool ul_early_term_ const int16_t ul_max_decoder_iter_ const
      uint16_t dl_base_graph_ const bool dl_early_term_ const int16_t dl_max_decoder_iter_

          // static constexpr size_t kCbPerSymbol = 1;
          // const size_t ofdm_data_num_ul_;
          // const size_t ofdm_data_num_dl_;
          // const size_t ofdm_ctrl_data_;

          // size_t dl_bcast_mod_order_bits_;

          // // The total number of downlink MAC payload data bytes in each Frame
          // size_t dl_mac_data_bytes_num_perframe_;
          // size_t ul_num_padding_bytes_per_cb_;
          // size_t dl_num_padding_bytes_per_cb_;
          // size_t dl_num_bytes_per_cb_;
          // size_t dl_mac_data_length_max_;
          // size_t dl_mac_bytes_num_perframe_;
          // size_t dl_mac_packet_length_;
          // size_t dl_mac_packets_perframe_;
          // // The length (in bytes) of a uplink MAC packet payload (data)
          // size_t ul_mac_data_length_max_;
          // size_t dl_data_bytes_num_persymbol_;

          // std::string ul_modulation_;
          // std::string dl_modulation_;
          // LDPCconfig dl_ldpc_config_;
          // LDPCconfig ul_ldpc_config_;
          // LDPCconfig dl_bcast_ldpc_config_;
          // FrameStats frame_;
          // const OfdmConfig ofdm_data_;
          // // The total number of uplink MAC payload data bytes in each Frame
          // size_t ul_mac_data_bytes_num_perframe_;
          // // The total number of uplink MAC packet bytes in each Frame
          // size_t ul_mac_bytes_num_perframe_;
          // // The length (in bytes) of a uplink MAC packet including the header
          // size_t ul_mac_packet_length_;
          // // The length (in bytes) of a uplink MAC packet payload (data)
          // size_t ul_mac_data_length_max_;
          // // The total number of downlink MAC payload data bytes in each Frame
          // size_t dl_mac_data_bytes_num_perframe_;
          // // The total number of downlink MAC packet bytes in each Frame
          // // The length (in bytes) of a downlink MAC packet including the header
          // size_t dl_mac_packet_length_;
          // // The length (in bytes) of a downlink MAC packet payload (data)
          // // The total number of downlink mac packets sent/received in each frame
          // size_t dl_mac_packets_perframe_;
          // // The total number of uplink mac packets sent/received in each frame
          // size_t ul_mac_packets_perframe_;
          // // size_t dl_bcast_mod_order_bits_;
          // // // Number of bytes per code block
          // // size_t ul_num_bytes_per_cb_;
          // // size_t dl_num_bytes_per_cb_;
          // // size_t ul_mcs_index_;
          // // size_t dl_mcs_index_;
          // // // Number of padding bytes per code block
          // // size_t ul_num_padding_bytes_per_cb_;
          // // size_t dl_num_padding_bytes_per_cb_;
          // // // The total number of uncoded uplink data bytes in each OFDM symbol
          // size_t ul_data_bytes_num_persymbol_;
          // // The total number of uncoded downlink data bytes in each OFDM symbol
          // size_t dl_data_bytes_num_persymbol_;
          // Table<std::complex<int16_t>> dl_iq_t_;
          // Table<std::complex<int16_t>> ul_iq_t_;
          // Table<std::complex<int16_t>> ue_specific_pilot_t_;
          // std::vector<std::string> ul_tx_f_data_files_;
          // std::vector<uint32_t> beacon_;
          // std::vector<uint32_t> coeffs_;
          // std::vector<uint32_t> pilot_;
          // /// I/Q samples of common pilot
          // std::vector<std::complex<int16_t>> pilot_ci16_;
          // std::vector<std::complex<int16_t>> beacon_ci16_;
          // std::vector<std::complex<float>> pilot_cf32_;
          // std::vector<std::complex<float>> gold_cf32_;
          // std::vector<std::complex<float>> common_pilot_;

          // /// I/Q samples of pilots per UE antenna per pilot symbol
          // std::vector<std::vector<std::vector<std::complex<int16_t>>>> pilot_ue_ci16_;
          // // List of subcarriers used per UE to transmit pilot
          // std::vector<arma::uvec> pilot_ue_sc_;
          // LDPCconfig ul_ldpc_config_;        //Uplink LDPC parameters
          // LDPCconfig dl_ldpc_config_;        //Downlink LDPC parameters
          // LDPCconfig dl_bcast_ldpc_config_;  // Downlink Broadcast LDPC parameters
          // nlohmann::json ul_mcs_params_;     // Uplink Modulation and Coding (MCS)
          // nlohmann::json dl_mcs_params_;     // Downlink Modulation and Coding (MCS)
          McsScheme next_ul_mcs_;
  McsScheme next_dl_mcs_;
  McsScheme current_ul_mcs_;
  McsScheme current_dl_mcs_;
  // ModulationTables modulation_tables_;
  // FrameStats frame_;
  // const InitialMcsProperties initial_ul_mcs_properties_;
  // const InitialMcsProperties initial_dl_mcs_properties_;

  UpdateUlLdpcConfig(LDPCconfig ul_ldpc_config_, size_t ul_mcs_index);
  UpdateDlLdpcConfig(LDPCconfig dl_ldpc_config_, size_t dl_mcs_index);
  UpdateCtrlMCS(LDPCconfig dl_bcast_ldpc_config_, size_t control_mcs_index);
  CalculateUlLdpcProperties(UlLdpcProperties ul_ldpc_properties,
                            LDPCconfig ul_ldpc_config);
  CalculateDlLdpcProperties(DlLdpcProperties dl_ldpc_properties,
                            LDPCconfig dl_ldpc_config);

  UlMcsParams ul_mcs_params_;
  DlMcsParams dl_mcs_params_;
}

#endif /* LDPC_UPDATE_H_ */