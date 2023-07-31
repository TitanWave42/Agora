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

const size_t kNumTables = 5;

struct MCS_Scheme {
  size_t frame_number;
  size_t mcs_index;
  size_t modulation_table_index;
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
  explicit Mcs(Config* const cfg);
  ~Mcs();

  void CreateModulationTables();
  void UpdateMcsSchemes(size_t current_frame_number);

  void Set_Next_Ul_MCS_Scheme(MCS_Scheme next_mcs_scheme);
  void Set_Next_Dl_MCS_Scheme(MCS_Scheme next_mcs_scheme);

  void Update_Ul_Ldpc_Config();
  void Update_Dl_Ldpc_Config();

  void InitializeUlMcs(const nlohmann::json ul_mcs);
  void InitializeDlMcs(const nlohmann::json dl_mcs);

  inline std::string Modulation(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_modulation_
                                     : this->dl_modulation_;
  }
  inline size_t ModOrderBits(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mod_order_bits_
                                     : this->dl_mod_order_bits_;
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
  inline Table<complex_float>& ModTable(Direction dir) {
    return dir == Direction::kUplink
               ? this->modulation_tables_
                     .ul_tables[current_ul_mcs_.modulation_table_index]
               : this->modulation_tables_
                     .dl_tables[current_dl_mcs_.modulation_table_index];
  }
  inline const nlohmann::json& MCSParams(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mcs_params_
                                     : this->dl_mcs_params_;
  }
  inline size_t SubcarrierPerCodeBlock(Direction dir) const {
    return this->LdpcConfig(dir).NumCbCodewLen() / this->ModOrderBits(dir);
  }
  inline size_t McsIndex(Direction dir) const {
    return dir == Direction::kUplink ? this->ul_mcs_index_
                                     : this->dl_mcs_index_;
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
  };

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
  };

  void GenData();
  void GenPilots();

  void DumpMcsInfo();
  void UpdateCtrlMCS();

  inline const arma::uvec& PilotUeSc(size_t ue_id) const {
    return this->pilot_ue_sc_.at(ue_id);
  }

  inline std::vector<std::complex<int16_t>>& BeaconCi16() {
    return this->beacon_ci16_;
  };

  LDPCconfig Ul_Ldpc_Config();
  LDPCconfig Dl_Ldpc_Config();

  inline std::vector<std::complex<int16_t>>& PilotUeCi16(size_t ue_id,
                                                         size_t pilot_idx) {
    return this->pilot_ue_ci16_.at(ue_id).at(pilot_idx);
  };

  inline const std::vector<uint32_t>& Pilot() const { return this->pilot_; };

  inline const std::vector<std::complex<float>>& PilotCf32() const {
    return this->pilot_cf32_;
  };

  inline float Scale() const { return this->scale_; }

  inline Table<std::complex<int16_t>>& DlIqT() { return this->dl_iq_t_; }

  inline const std::vector<std::string>& UlTxFreqDataFiles() const {
    return ul_tx_f_data_files_;
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

  inline Table<int8_t>& DlModBits() { return this->dl_mod_bits_; }
  inline Table<int8_t>& UlModBits() { return this->ul_mod_bits_; }
  inline Table<int8_t>& DlBits() { return this->dl_bits_; }
  inline Table<int8_t>& UlBits() { return this->ul_bits_; }
  inline Table<std::complex<int16_t>>& UlIqT() { return this->ul_iq_t_; }

  inline Table<complex_float>& UlIqF() { return this->ul_iq_f_; }
  inline Table<complex_float>& DlIqF() { return this->dl_iq_f_; }

  inline const complex_float* Pilots(void) const { return this->pilots_; };
  inline const complex_float* PilotsSgn() const { return this->pilots_sgn_; };

  size_t DecodeBroadcastSlots(const int16_t* const bcast_iq_samps);
  void GenBroadcastSlots(std::vector<std::complex<int16_t>*>& bcast_iq_samps,
                         std::vector<size_t> ctrl_msg);

 private:
  complex_float* pilots_;
  complex_float* pilots_sgn_;

  Table<complex_float> dl_iq_f_;
  Table<complex_float> ul_iq_f_;

  Table<int8_t> ul_bits_;

  Table<int8_t> dl_bits_;

  Table<int8_t> ul_mod_bits_;

  std::vector<std::string> ul_tx_f_data_files_;

  float scale_;  // Scaling factor for all transmit symbols

  Table<std::complex<int16_t>> dl_iq_t_;
  Table<std::complex<int16_t>> ul_iq_t_;

  Table<int8_t> dl_mod_bits_;

  std::vector<uint32_t> pilot_;

  std::vector<std::complex<float>> pilot_cf32_;

  /// I/Q samples of pilots per UE antenna per pilot symbol
  std::vector<std::vector<std::vector<std::complex<int16_t>>>> pilot_ue_ci16_;

  complex_float* pilot_ifft_;

  std::vector<std::complex<float>> gold_cf32_;

  static constexpr size_t kCbPerSymbol = 1;

  size_t beacon_len_;

  std::vector<uint32_t> beacon_;
  std::vector<uint32_t> coeffs_;

  std::vector<std::complex<float>> common_pilot_;

  std::vector<std::complex<int16_t>> beacon_ci16_;

  LDPCconfig ul_ldpc_config_;  //Uplink LDPC parameters
  LDPCconfig dl_ldpc_config_;  //Downlink LDPC parameters

  /// List of subcarriers used per UE to transmit pilot
  std::vector<arma::uvec> pilot_ue_sc_;

  Table<complex_float> ue_specific_pilot_;
  Table<std::complex<int16_t>> ue_specific_pilot_t_;

  LDPCconfig dl_bcast_ldpc_config_;  // Downlink Broadcast LDPC parameters

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
  size_t
      ul_mod_order_bits_;  // Number of binary bits used for a modulation order
  std::string dl_modulation_;
  size_t dl_mod_order_bits_;
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

  MCS_Scheme current_ul_mcs_;
  MCS_Scheme current_dl_mcs_;

  nlohmann::json ul_mcs_params_;  // Uplink Modulation and Coding (MCS)
  nlohmann::json dl_mcs_params_;  // Downlink Modulation and Coding (MCS)

  MCS_Scheme next_ul_mcs_;
  MCS_Scheme next_dl_mcs_;

  Modulation_Tables modulation_tables_;

  Table<complex_float> ue_pilot_ifft_;

  Config* const cfg_;
  FrameStats frame_;
  size_t ofdm_ca_num_;

  Initial_Mcs_Properties initial_ul_mcs_properties_;
  Initial_Mcs_Properties initial_dl_mcs_properties_;

  void UpdateUlMcsScheme(size_t current_frame_number);
  void UpdateDlMcsScheme(size_t current_frame_number);
  void Update_Ldpc_Properties();
  void CalculateLdpcProperties();
};

#endif  // MCS_H_
