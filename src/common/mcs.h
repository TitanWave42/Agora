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
constexpr float snr_delta_tolerance = 0.1;

static const std::map<size_t, float> kmcs_index_to_spectral_effeciency = {
    {0, 0.2344},  {1, 0.3066},  {2, 0.3770},  {3, 0.4902},  {4, 0.6016},
    {5, 0.7402},  {6, 0.8770},  {7, 1.0273},  {8, 1.1758},  {9, 1.3262},
    {10, 1.3281}, {11, 1.4766}, {12, 1.6953}, {13, 1.9141}, {14, 2.1602},
    {15, 2.4063}, {16, 2.5703}, {17, 2.5664}, {18, 2.7305}, {19, 3.0293},
    {20, 3.3223}, {21, 3.6094}, {22, 3.9023}, {23, 4.2129}, {24, 4.5234},
    {25, 4.8164}, {26, 5.1152}, {27, 5.3320}, {28, 5.5547}, {29, 5.8906},
    {30, 6.2266}, {31, 6.5703}};

struct MCS_Scheme {
  size_t frame_number;
  size_t mcs_index;
  size_t mod_order_bits;
  size_t code_rate;
};

struct Modulation_Tables {
  Table<complex_float> ul_tables[kNumTables];
  Table<complex_float> dl_tables[kNumTables];
};

struct Initial_Mcs_Properties {
  uint16_t base_graph;
  bool early_term;
  int16_t max_decoder_iter;
  size_t ofdm_data_num;
};

class Mcs {
 public:
  explicit Mcs(Config* const cfg);
  ~Mcs();

  void CreateModulationTables();
  void InitializeUlMcs(const nlohmann::json ul_mcs);
  void InitializeDlMcs(const nlohmann::json dl_mcs);
  void UpdateMcs(size_t current_frame_number);
  void SetNextUlMcs(size_t frame_number, size_t mod_order_bits);
  void SetNextDlMcs(size_t frame_number, size_t mod_order_bits);
  void UpdateUlLdpcConfig();
  void UpdateDlLdpcConfig();

  inline void Running(bool value) { this->running_.store(value); }
  inline bool Running() const { return this->running_.load(); }
  inline MCS_Scheme CurrentDlMcs() { return this->current_dl_mcs_; }

  //spectral effeciency calculation is log_2(1+0.902)
  inline float SpectralEffeciency(float snr) { return log2(1 + snr); }

  inline const complex_float* Pilots(void) const { return this->pilots_; }
  inline const complex_float* PilotsSgn() const { return this->pilots_sgn_; }
  inline const std::vector<uint32_t>& Beacon() const { return this->beacon_; }
  inline const std::vector<uint32_t>& Coeffs() const { return this->coeffs_; }
  inline std::vector<std::complex<int16_t>>& PilotCi16() {
    return this->pilot_ci16_;
  }

  inline Table<complex_float>& UeSpecificPilot() {
    return this->ue_specific_pilot_;
  }
  inline Table<std::complex<int16_t>>& UeSpecificPilotT() {
    return this->ue_specific_pilot_t_;
  }

  inline const std::vector<std::complex<float>>& GoldCf32() const {
    return this->gold_cf32_;
  }

  inline size_t BeaconLen() const { return this->beacon_len_; }

  inline const arma::uvec& PilotUeSc(size_t ue_id) const {
    return this->pilot_ue_sc_.at(ue_id);
  }
  inline std::vector<std::complex<int16_t>>& BeaconCi16() {
    return this->beacon_ci16_;
  }
  inline std::vector<std::complex<int16_t>>& PilotUeCi16(size_t ue_id,
                                                         size_t pilot_idx) {
    return this->pilot_ue_ci16_.at(ue_id).at(pilot_idx);
  }
  inline const std::vector<uint32_t>& Pilot() const { return this->pilot_; };
  inline const std::vector<std::complex<float>>& PilotCf32() const {
    return this->pilot_cf32_;
  }
  inline float Scale() const { return this->scale_; }
  inline Table<std::complex<int16_t>>& DlIqT() { return this->dl_iq_t_; }
  inline const std::vector<std::string>& UlTxFreqDataFiles() const {
    return ul_tx_f_data_files_;
  }
  inline std::string Modulation(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_modulation_
                                     : this->dl_modulation_;
  }
  inline size_t McsUpdateFrame(Direction dir) {
    return dir == Direction::kUplink ? next_ul_mcs_.frame_number
                                     : next_dl_mcs_.frame_number;
  }
  inline size_t ModOrderBits(Direction dir) const {
    return dir == Direction::kUplink ? this->current_ul_mcs_.mod_order_bits
                                     : this->current_dl_mcs_.mod_order_bits;
  }
  inline size_t NumBytesPerCb(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_num_bytes_per_cb_
                                     : this->dl_num_bytes_per_cb_;
  }
  inline size_t NumPaddingBytesPerCb(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_num_padding_bytes_per_cb_
                                     : this->dl_num_padding_bytes_per_cb_;
  }
  inline size_t MacDataBytesNumPerframe(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mac_data_bytes_num_perframe_
                                     : this->dl_mac_data_bytes_num_perframe_;
  }
  inline size_t MacBytesNumPerframe(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mac_bytes_num_perframe_
                                     : this->dl_mac_bytes_num_perframe_;
  }
  inline size_t MacPacketLength(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mac_packet_length_
                                     : this->dl_mac_packet_length_;
  }
  inline size_t MacPayloadMaxLength(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mac_data_length_max_
                                     : this->dl_mac_data_length_max_;
  }
  inline size_t MacPacketsPerframe(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mac_packets_perframe_
                                     : this->dl_mac_packets_perframe_;
  }
  inline const LDPCconfig& LdpcConfig(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_ldpc_config_
                                     : this->dl_ldpc_config_;
  }
  inline const LDPCconfig& BcLdpcConfig() const {
    return dl_bcast_ldpc_config_;
  }
  inline const nlohmann::json& McsParams(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mcs_params_
                                     : this->dl_mcs_params_;
  }
  inline size_t SubcarrierPerCodeBlock(Direction dir) const {
    return this->LdpcConfig(dir).NumCbCodewLen() / this->ModOrderBits(dir);
  }
  inline size_t McsIndex(Direction dir) const {
    return dir == Direction::kUplink ? current_ul_mcs_.mcs_index
                                     : current_dl_mcs_.mcs_index;
  }
  inline Table<int8_t>& DlModBits() { return this->dl_mod_bits_; }
  inline Table<int8_t>& UlModBits() { return this->ul_mod_bits_; }
  inline Table<int8_t>& DlBits() { return this->dl_bits_; }
  inline Table<int8_t>& UlBits() { return this->ul_bits_; }
  inline Table<std::complex<int16_t>>& UlIqT() { return this->ul_iq_t_; }
  inline Table<complex_float>& UlIqF() { return this->ul_iq_f_; }
  inline Table<complex_float>& DlIqF() { return this->dl_iq_f_; }
  inline Table<complex_float>& ModTable(Direction dir) {
    // std::cout << "DL Mod Order Bits: " << current_dl_mcs_.mod_order_bits
    //           << std::endl
    //           << std::flush;
    // std::cout << "DL Mcs index: " << current_dl_mcs_.mcs_index << std::endl
    //           << std::flush;
    // std::cout
    //     << "Size DL: "
    //     << sizeof(modulation_tables_.dl_tables[current_dl_mcs_.mod_order_bits])
    //     << std::endl
    //     << std::flush;

    // std::cout << "UL Mod Order bits " << current_ul_mcs_.mod_order_bits
    //           << std::endl
    //           << std::flush;

    // std::cout
    //     << "Ul Mod Table: "
    //     << (modulation_tables_.ul_tables[current_ul_mcs_.mod_order_bits])[0][0]
    //            .re
    //     << " "
    //     << (modulation_tables_.ul_tables[current_ul_mcs_.mod_order_bits])[0][0]
    //            .im
    //     << std::endl
    //     << std::flush;

    // std::cout << "checking Dl Mod Table at position 0: "
    //           << (modulation_tables_.dl_tables[0])[0][0].re << " "
    //           << (modulation_tables_.dl_tables[0])[0][0].im << std::endl
    //           << std::flush;
    // std::cout
    //     << "Size UL: "
    //     << sizeof(modulation_tables_.ul_tables[current_ul_mcs_.mod_order_bits])
    //     << std::endl
    //     << std::flush;

    // std::cout << "current_ul_mcs_.mod_order_bits: "
    //           << current_ul_mcs_.mod_order_bits << std::endl
    //           << std::flush;
    // std::cout << " current_dl_mcs_.mod_order_bits: "
    //           << current_dl_mcs_.mod_order_bits << std::endl
    //           << std::flush;
    return dir == Direction::kUplink
               ? this->modulation_tables_
                     .ul_tables[current_ul_mcs_.mod_order_bits / 2 - 1]
               : this->modulation_tables_
                     .dl_tables[current_dl_mcs_.mod_order_bits / 2 - 1];
  }
  /// Get info bits for this symbol, user and code block ID
  inline int8_t* GetInfoBits(Table<int8_t>& info_bits, Direction dir,
                             size_t symbol_id, size_t ue_id,
                             size_t cb_id) const {
    size_t num_bytes_per_cb;
    size_t num_blocks_in_symbol;
    if (dir == Direction::kDownlink) {
      num_bytes_per_cb = this->dl_num_bytes_per_cb_;
      num_blocks_in_symbol = this->dl_ldpc_config_.NumBlocksInSymbol();
    } else {
      num_bytes_per_cb = this->ul_num_bytes_per_cb_;
      num_blocks_in_symbol = this->ul_ldpc_config_.NumBlocksInSymbol();
    }
    return &info_bits[symbol_id][Roundup<64>(num_bytes_per_cb) *
                                 (num_blocks_in_symbol * ue_id + cb_id)];
  }
  //EVENTUALLY UPDATE THIS FUNCTION TO TAKE num_byte_per_cb as an argument
  //because the num_bytes_per_cb is determined by the MCS and is effectively
  //How much data we can read from the date buffer.
  /// Get info bits for this symbol, user and code block ID
  /// Get mac bits for this frame, symbol, user and code block ID
  inline int8_t* GetMacBits(Table<int8_t>& info_bits, Direction dir,
                            size_t frame_id, size_t symbol_id, size_t ue_id,
                            size_t cb_id) const {
    size_t mac_bytes_perframe;
    size_t num_bytes_per_cb;
    size_t mac_packet_length;
    if (dir == Direction::kDownlink) {
      mac_bytes_perframe = this->dl_mac_bytes_num_perframe_;
      num_bytes_per_cb = this->dl_num_bytes_per_cb_;
      mac_packet_length = this->dl_mac_packet_length_;
    } else {
      mac_bytes_perframe = ul_mac_bytes_num_perframe_;
      num_bytes_per_cb = this->ul_num_bytes_per_cb_;
      mac_packet_length = this->ul_mac_packet_length_;
    }
    return &info_bits[ue_id][(frame_id % kFrameWnd) * mac_bytes_perframe +
                             symbol_id * mac_packet_length +
                             cb_id * num_bytes_per_cb];
  }
  /// Get encoded_buffer for this frame, symbol, user and code block ID
  inline int8_t* GetModBitsBuf(Table<int8_t>& mod_bits_buffer, Direction dir,
                               size_t frame_id, size_t symbol_id, size_t ue_id,
                               size_t sc_id) const {
    size_t total_data_symbol_id;
    size_t ofdm_data_num;
    if (dir == Direction::kDownlink) {
      total_data_symbol_id = cfg_->GetTotalDataSymbolIdxDl(frame_id, symbol_id);
      ofdm_data_num = cfg_->GetOFDMDataNum();
    } else {
      total_data_symbol_id = cfg_->GetTotalDataSymbolIdxUl(frame_id, symbol_id);
      ofdm_data_num = cfg_->OfdmDataNum();
    }

    return &mod_bits_buffer[total_data_symbol_id]
                           [Roundup<64>(ofdm_data_num) * ue_id + sc_id];
  }

  void GenData();
  void GenPilots();
  void DumpMcsInfo();
  void UpdateCtrlMCS();
  LDPCconfig UlLdpcConfig();
  LDPCconfig DlLdpcConfig();
  size_t DecodeBroadcastSlots(const int16_t* const bcast_iq_samps);
  void GenBroadcastSlots(std::vector<std::complex<int16_t>*>& bcast_iq_samps,
                         std::vector<size_t> ctrl_msg);
  void CheckUlMcs(float snr, size_t frame_id);

 private:
  static constexpr size_t kCbPerSymbol = 1;
  std::atomic<bool> running_;

  size_t beacon_len_;
  // The total number of uplink MAC payload data bytes in each Frame
  size_t ul_mac_data_bytes_num_perframe_;
  // The total number of uplink MAC packet bytes in each Frame
  size_t ul_mac_bytes_num_perframe_;
  // The length (in bytes) of a uplink MAC packet including the header
  size_t ul_mac_packet_length_;
  // The length (in bytes) of a uplink MAC packet payload (data)
  size_t ul_mac_data_length_max_;
  // The total number of downlink MAC payload data bytes in each Frame
  size_t dl_mac_data_bytes_num_perframe_;
  // The total number of downlink MAC packet bytes in each Frame
  size_t dl_mac_bytes_num_perframe_;
  // The length (in bytes) of a downlink MAC packet including the header
  size_t dl_mac_packet_length_;
  // The length (in bytes) of a downlink MAC packet payload (data)
  size_t dl_mac_data_length_max_;
  // The total number of downlink mac packets sent/received in each frame
  size_t dl_mac_packets_perframe_;
  // The total number of uplink mac packets sent/received in each frame
  size_t ul_mac_packets_perframe_;
  std::string ul_modulation_;  // Modulation order as a string, e.g., "16QAM"
  std::string dl_modulation_;
  size_t dl_bcast_mod_order_bits_;
  // Number of bytes per code block
  size_t ul_num_bytes_per_cb_;
  size_t dl_num_bytes_per_cb_;
  size_t ul_mcs_index_;
  size_t dl_mcs_index_;
  // Number of padding bytes per code block
  size_t ul_num_padding_bytes_per_cb_;
  size_t dl_num_padding_bytes_per_cb_;
  // The total number of uncoded uplink data bytes in each OFDM symbol
  size_t ul_data_bytes_num_persymbol_;
  // The total number of uncoded downlink data bytes in each OFDM symbol
  size_t dl_data_bytes_num_persymbol_;
  size_t ofdm_ca_num_;
  float scale_;  // Scaling factor for all transmit symbols
  complex_float* pilots_;
  complex_float* pilots_sgn_;
  complex_float* pilot_ifft_;
  Table<complex_float> dl_iq_f_;
  Table<complex_float> ul_iq_f_;
  Table<complex_float> ue_specific_pilot_;
  Table<complex_float> ue_pilot_ifft_;
  Table<int8_t> ul_bits_;
  Table<int8_t> dl_bits_;
  Table<int8_t> ul_mod_bits_;
  Table<int8_t> dl_mod_bits_;
  Table<std::complex<int16_t>> dl_iq_t_;
  Table<std::complex<int16_t>> ul_iq_t_;
  Table<std::complex<int16_t>> ue_specific_pilot_t_;
  std::vector<std::string> ul_tx_f_data_files_;
  std::vector<uint32_t> beacon_;
  std::vector<uint32_t> coeffs_;
  std::vector<uint32_t> pilot_;
  /// I/Q samples of common pilot
  std::vector<std::complex<int16_t>> pilot_ci16_;
  std::vector<std::complex<int16_t>> beacon_ci16_;
  std::vector<std::complex<float>> pilot_cf32_;
  std::vector<std::complex<float>> gold_cf32_;
  std::vector<std::complex<float>> common_pilot_;

  /// I/Q samples of pilots per UE antenna per pilot symbol
  std::vector<std::vector<std::vector<std::complex<int16_t>>>> pilot_ue_ci16_;
  // List of subcarriers used per UE to transmit pilot
  std::vector<arma::uvec> pilot_ue_sc_;
  LDPCconfig ul_ldpc_config_;        //Uplink LDPC parameters
  LDPCconfig dl_ldpc_config_;        //Downlink LDPC parameters
  LDPCconfig dl_bcast_ldpc_config_;  // Downlink Broadcast LDPC parameters
  nlohmann::json ul_mcs_params_;     // Uplink Modulation and Coding (MCS)
  nlohmann::json dl_mcs_params_;     // Downlink Modulation and Coding (MCS)
  MCS_Scheme next_ul_mcs_;
  MCS_Scheme next_dl_mcs_;
  MCS_Scheme current_ul_mcs_;
  MCS_Scheme current_dl_mcs_;
  Modulation_Tables modulation_tables_;
  Config* const cfg_;
  FrameStats frame_;
  Initial_Mcs_Properties initial_ul_mcs_properties_;
  Initial_Mcs_Properties initial_dl_mcs_properties_;

  void UpdateUlMcs(size_t current_frame_number);
  void UpdateDlMcs(size_t current_frame_number);
  void CalculateLdpcProperties();
};

#endif  // MCS_H_
