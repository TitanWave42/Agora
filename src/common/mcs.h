/**
 * @file mcs.h
 * @brief Class definition for mcs.
 */
#ifndef MCS_H_
#define MCS_H_

#include "config.h"
#include "framestats.h"
#include "ldpc_config.h"
#include "memory_manage.h"
#include "nlohmann/json.hpp"
#include "utils.h"

//#include "comms-constants.inc"

constexpr size_t kNumTables = 4;

//For determining if we need to change mcs
constexpr float kSnrDeltaTolerance = 0.1;

static const std::map<size_t, float> kMcsIndexToSpectralEffeciency = {
    {0, 0.2344},  {1, 0.3066},  {2, 0.3770},  {3, 0.4902},  {4, 0.6016},
    {5, 0.7402},  {6, 0.8770},  {7, 1.0273},  {8, 1.1758},  {9, 1.3262},
    {10, 1.3281}, {11, 1.4766}, {12, 1.6953}, {13, 1.9141}, {14, 2.1602},
    {15, 2.4063}, {16, 2.5703}, {17, 2.5664}, {18, 2.7305}, {19, 3.0293},
    {20, 3.3223}, {21, 3.6094}, {22, 3.9023}, {23, 4.2129}, {24, 4.5234},
    {25, 4.8164}, {26, 5.1152}, {27, 5.3320}, {28, 5.5547}, {29, 5.8906},
    {30, 6.2266}, {31, 6.5703}};

struct McsScheme {
  size_t frame_number_;
  size_t mcs_index_;
  size_t mod_order_bits_;
  size_t code_rate_;
  Direction dir;
};

struct ModulationTables {
  Table<complex_float> ul_tables_[kNumTables];
  Table<complex_float> dl_tables_[kNumTables];
};


class Mcs {
 public:
  explicit Mcs();
  ~Mcs();

  void CheckUlMcs(float snr, size_t frame_id);

 private:
  static constexpr size_t kCbPerSymbol = 1;
  const size_t ofdm_data_num_ul_;
  const size_t ofdm_data_num_dl_;
  const size_t ofdm_ctrl_data_;

  //   size_t dl_bcast_mod_order_bits_;

  // // The total number of downlink MAC payload data bytes in each Frame
  //   size_t dl_mac_data_bytes_num_perframe_;
  //   size_t ul_num_padding_bytes_per_cb_;
  //   size_t dl_num_padding_bytes_per_cb_;
  //   size_t dl_num_bytes_per_cb_;
  //   size_t dl_mac_data_length_max_;
  //   size_t dl_mac_bytes_num_perframe_;
  //   size_t dl_mac_packet_length_;
  //   size_t dl_mac_packets_perframe_;
  //   // The length (in bytes) of a uplink MAC packet payload (data)
  //   size_t ul_mac_data_length_max_;
  //   size_t dl_data_bytes_num_persymbol_;

  std::string ul_modulation_;
  std::string dl_modulation_;
  LDPCconfig dl_ldpc_config_;
  LDPCconfig ul_ldpc_config_;
  LDPCconfig dl_bcast_ldpc_config_;
  FrameStats frame_;
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
  // size_t dl_bcast_mod_order_bits_;
  // // Number of bytes per code block
  // size_t ul_num_bytes_per_cb_;
  // size_t dl_num_bytes_per_cb_;
  // size_t ul_mcs_index_;
  // size_t dl_mcs_index_;
  // // Number of padding bytes per code block
  // size_t ul_num_padding_bytes_per_cb_;
  // size_t dl_num_padding_bytes_per_cb_;
  // // The total number of uncoded uplink data bytes in each OFDM symbol
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


  std::unique_ptr<LdpcUpdater> ldpc_updater;

  void UpdateUlMcs(size_t current_frame_number);
  void UpdateDlMcs(size_t current_frame_number);
  void CalculateLdpcProperties();
  void CreateModulationTables();

  void SetNextUlMcs(size_t frame_number, size_t mod_order_bits);
  void SetNextDlMcs(size_t frame_number, size_t mod_order_bits);
  void InitializeUlMcs(const nlohmann::json ul_mcs);
  void InitializeDlMcs(const nlohmann::json dl_mcs);
  void UpdateUlLdpcConfig();
  void UpdateDlLdpcConfig();

  void UpdateMcs(size_t current_frame_number);
  void DumpMcsInfo();
  void UpdateCtrlMCS();
};

#endif /* MCS_H_ */
