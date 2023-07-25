/**
 * @file mcs.cc
 * @brief Class implementation for mcs handling
 */

#include "mcs.h"

#include <stddef.h>

#include "comms-lib.h"
#include "config.h"
#include "ldpc_config.h"
#include "logger.h"
#include "modulation.h"
#include "utils.h"
#include "utils_ldpc.h"

Mcs::Mcs(Config* const cfg)
    : ul_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      dl_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      cfg_(cfg) {
  nlohmann::json ul_mcs_params = cfg_->UlMcsParams();
  nlohmann::json dl_mcs_params = cfg_->DlMcsParams();
  size_t ofdm_data_num = cfg_->OfdmCaNum;

  ul_mcs_ = cfg_->UlMcsParams();
  dl_mcs_ = cfg_->DlMcsParams();

  Create_Modulation_Tables();

  this->frame_ = cfg_->Frame();
  ofdm_ca_num_ = cfg_->ofdm_ca_num_;

  initial_ul_mcs_properties_.base_graph = ul_mcs_params.value("base_graph", 1);
  initial_ul_mcs_properties_.early_term =
      ul_mcs_params.value("earlyTermination", true);
  initial_ul_mcs_properties_.max_decoder_iter =
      ul_mcs_params.value("decoderIter", 5);
  initial_ul_mcs_properties_.ofdm_data_num = ofdm_data_num;

  initial_dl_mcs_properties_.base_graph = dl_mcs_params.value("base_graph", 1);
  initial_dl_mcs_properties_.early_term =
      dl_mcs_params.value("earlyTermination", true);
  initial_dl_mcs_properties_.max_decoder_iter =
      dl_mcs_params.value("decoderIter", 5);

  //Initialize UL MCS
  Initialize_Ul_Mcs(ul_mcs_);

  //Initialize DL MCS
  Initialize_Dl_Mcs(dl_mcs_);

  //Update the LDPC configs
  Update_Ul_Ldpc_Config();
  Update_Dl_Ldpc_Config();

  CalculateLdpcProperties();
}

Mcs::~Mcs() = default;

void Mcs::Initialize_Ul_Mcs(const nlohmann::json ul_mcs) {
  current_ul_mcs_.frame_number = 0;
  if (ul_mcs.find("mcs_index") == ul_mcs.end()) {
    std::string ul_modulation_type = ul_mcs.value("modulation", "16QAM");
    current_ul_mcs_.modulation = kModulStringMap.at(ul_modulation_type);

    double ul_code_rate_usr = ul_mcs.value("code_rate", 0.333);
    size_t code_rate_int =
        static_cast<size_t>(std::round(ul_code_rate_usr * 1024.0));

    current_ul_mcs_.mcs_index =
        CommsLib::GetMcsIndex(current_ul_mcs_.modulation, code_rate_int);
    current_ul_mcs_.code_rate = GetCodeRate(current_ul_mcs_.mcs_index);
    if (current_ul_mcs_.code_rate / 1024.0 != ul_code_rate_usr) {
      AGORA_LOG_WARN(
          "Rounded the user-defined uplink code rate to the closest standard "
          "rate %zu/1024.\n",
          current_ul_mcs_.code_rate);
    }
  } else {
    current_ul_mcs_.mcs_index =
        ul_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
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
    current_dl_mcs_.mcs_index =
        CommsLib::GetMcsIndex(current_dl_mcs_.modulation, code_rate_int);
    current_dl_mcs_.code_rate = GetCodeRate(current_dl_mcs_.mcs_index);
    if (current_dl_mcs_.code_rate / 1024.0 != dl_code_rate_usr) {
      AGORA_LOG_WARN(
          "Rounded the user-defined downlink code rate to the closest standard "
          "rate %zu/1024.\n",
          current_dl_mcs_.code_rate);
    }
  } else {
    current_dl_mcs_.mcs_index =
        dl_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_dl_mcs_.modulation = GetModOrderBits(current_dl_mcs_.mcs_index);
    current_dl_mcs_.code_rate = GetCodeRate(current_dl_mcs_.mcs_index);
  }
}

void Mcs::Create_Modulation_Tables() {
  for (int i = 0; i < modulation_tables_.dl_tables.size(); i++) {
    InitModulationTable(modulation_tables_.dl_tables[i], (i + 1) * 2);
    InitModulationTable(modulation_tables_.ul_tables[i], (i + 1) * 2);
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
  uint16_t base_graph = initial_ul_mcs_properties_.base_graph;

  size_t ul_mod_order_bits = GetModOrderBits(current_ul_mcs_.mcs_index);
  size_t ul_code_rate = GetCodeRate(current_ul_mcs_.mcs_index);

  Table<complex_float> ul_mod_table_ =
      modulation_tables_.ul_tables[ul_mod_order_bits / 2 - 1];

  size_t zc = SelectZc(base_graph, ul_code_rate, ul_mod_order_bits,
                       cfg_->OfdmDataNum, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows = static_cast<size_t>(std::round(
                        1024.0 * LdpcNumInputCols(base_graph) / ul_code_rate)) -
                    (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  ul_ldpc_config_ =
      LDPCconfig(base_graph, zc, initial_ul_mcs_properties_.max_decoder_iter,
                 initial_ul_mcs_properties_.early_term, num_cb_len,
                 num_cb_codew_len, num_rows, 0);

  ul_ldpc_config_.NumBlocksInSymbol((cfg_->OfdmDataNum() * ul_mod_order_bits) /
                                    ul_ldpc_config_.NumCbCodewLen());
  RtAssert((frame_.NumULSyms() == 0) ||
               (ul_ldpc_config_.NumBlocksInSymbol() > 0),
           "Uplink LDPC expansion factor is too large for number of OFDM data "
           "subcarriers.");
}

void Mcs::Update_Dl_Ldpc_Config() {
  uint16_t base_graph = initial_dl_mcs_properties_.base_graph;

  size_t dl_mod_order_bits = GetModOrderBits(current_dl_mcs_.mcs_index);
  size_t dl_code_rate = GetCodeRate(current_dl_mcs_.mcs_index);

  Table<complex_float> dl_mod_table_ =
      modulation_tables_.dl_tables[dl_mod_order_bits / 2 - 1];

  size_t zc = SelectZc(base_graph, dl_code_rate, dl_mod_order_bits,
                       cfg_->OfdmDataNum, Config::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows = static_cast<size_t>(std::round(
                        1024.0 * LdpcNumInputCols(base_graph) / dl_code_rate)) -
                    (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  dl_ldpc_config_ =
      LDPCconfig(base_graph, zc, initial_ul_mcs_properties_.max_decoder_iter,
                 initial_ul_mcs_properties_.early_term, num_cb_len,
                 num_cb_codew_len, num_rows, 0);

  dl_ldpc_config_.NumBlocksInSymbol(
      (initial_dl_mcs_properties_.ofdm_data_num * dl_mod_order_bits) /
      ul_ldpc_config_.NumCbCodewLen());
  RtAssert((frame_.NumULSyms() == 0) ||
               (ul_ldpc_config_.NumBlocksInSymbol() > 0),
           "Uplink LDPC expansion factor is too large for number of OFDM data "
           "subcarriers.");
}

void Mcs::CalculateLdpcProperties() {
  this->ul_num_bytes_per_cb_ = ul_ldpc_config_.NumCbLen() / 8;
  this->ul_num_padding_bytes_per_cb_ =
      Roundup<64>(ul_num_bytes_per_cb_) - ul_num_bytes_per_cb_;
  this->ul_data_bytes_num_persymbol_ =
      ul_num_bytes_per_cb_ * ul_ldpc_config_.NumBlocksInSymbol();
  size_t ul_mac_packet_length_ = ul_data_bytes_num_persymbol_;

  //((cb_len_bits / zc_size) - 1) * (zc_size / 8) + kProcBytes(32)
  const size_t ul_ldpc_input_min =
      (((ul_ldpc_config_.NumCbLen() / ul_ldpc_config_.ExpansionFactor()) - 1) *
           (ul_ldpc_config_.ExpansionFactor() / 8) +
       32);
  const size_t ul_ldpc_sugg_input = LdpcEncodingInputBufSize(
      ul_ldpc_config_.BaseGraph(), ul_ldpc_config_.ExpansionFactor());

  if (ul_ldpc_input_min >
      (ul_num_bytes_per_cb_ + ul_num_padding_bytes_per_cb_)) {
    // Can cause a lot of wasted space, specifically the second argument of the max
    const size_t increased_padding =
        Roundup<64>(ul_ldpc_sugg_input) - ul_num_bytes_per_cb_;

    AGORA_LOG_WARN(
        "LDPC required Input Buffer size exceeds uplink code block size!, "
        "Increased cb padding from %zu to %zu uplink CB Bytes %zu, LDPC "
        "Input Min for zc 64:256: %zu\n",
        ul_num_padding_bytes_per_cb_, increased_padding,
        ul_num_bytes_per_cb_, ul_ldpc_input_min);
    ul_num_padding_bytes_per_cb_ = increased_padding;
  }

  // Smallest over the air packet structure
  RtAssert(frame_.NumULSyms() == 0 ||
               ul_mac_packet_length_ > sizeof(MacPacketHeaderPacked),
           "Uplink MAC Packet size must be larger than MAC header size");
  cfg_->ul_mac_data_length_max_ =
      ul_mac_packet_length_ - sizeof(MacPacketHeaderPacked);

  cfg_->ul_mac_packets_perframe_ = frame_.NumUlDataSyms();
  cfg_->ul_mac_data_bytes_num_perframe_ =
      cfg_->ul_mac_data_length_max_ * cfg_->ul_mac_packets_perframe_;
  cfg_->ul_mac_bytes_num_perframe_ =
      ul_mac_packet_length_ * cfg_->ul_mac_packets_perframe_;

  cfg_->dl_num_bytes_per_cb_ = dl_ldpc_config_.NumCbLen() / 8;
  cfg_->dl_num_padding_bytes_per_cb_ =
      Roundup<64>(cfg_->dl_num_bytes_per_cb_) - cfg_->dl_num_bytes_per_cb_;
  cfg_->dl_data_bytes_num_persymbol_ =
      cfg_->dl_num_bytes_per_cb_ * dl_ldpc_config_.NumBlocksInSymbol();
  cfg_->dl_mac_packet_length_ = cfg_->dl_data_bytes_num_persymbol_;
  // Smallest over the air packet structure
  RtAssert(frame_.NumDLSyms() == 0 ||
               cfg_->dl_mac_packet_length_ > sizeof(MacPacketHeaderPacked),
           "Downlink MAC Packet size must be larger than MAC header size");
  cfg_->dl_mac_data_length_max_ =
      cfg_->dl_mac_packet_length_ - sizeof(MacPacketHeaderPacked);

  cfg_->dl_mac_packets_perframe_ = frame_.NumDlDataSyms();
  cfg_->dl_mac_data_bytes_num_perframe_ =
      cfg_->dl_mac_data_length_max_ * cfg_->dl_mac_packets_perframe_;
  cfg_->dl_mac_bytes_num_perframe_ =
      cfg_->dl_mac_packet_length_ * cfg_->dl_mac_packets_perframe_;

  //((cb_len_bits / zc_size) - 1) * (zc_size / 8) + kProcBytes(32)
  const size_t dl_ldpc_input_min =
      (((dl_ldpc_config_.NumCbLen() / dl_ldpc_config_.ExpansionFactor()) - 1) *
           (dl_ldpc_config_.ExpansionFactor() / 8) +
       32);
  const size_t dl_ldpc_sugg_input = LdpcEncodingInputBufSize(
      dl_ldpc_config_.BaseGraph(), dl_ldpc_config_.ExpansionFactor());

  if (dl_ldpc_input_min >
      (cfg_->dl_num_bytes_per_cb_ + cfg_->dl_num_padding_bytes_per_cb_)) {
    // Can cause a lot of wasted space, specifically the second argument of the max
    const size_t increased_padding =
        Roundup<64>(dl_ldpc_sugg_input) - cfg_->dl_num_bytes_per_cb_;

    AGORA_LOG_WARN(
        "LDPC required Input Buffer size exceeds downlink code block size!, "
        "Increased cb padding from %zu to %zu Downlink CB Bytes %zu, LDPC "
        "Input Min for zc 64:256: %zu\n",
        cfg_->dl_num_padding_bytes_per_cb_, increased_padding,
        cfg_->dl_num_bytes_per_cb_, dl_ldpc_input_min);
    cfg_->dl_num_padding_bytes_per_cb_ = increased_padding;
  }
}

void Config::GenData() {
  this->GenPilots();
  // Get uplink and downlink raw bits either from file or random numbers
  const size_t dl_num_bytes_per_ue_pad =
      Roundup<64>(this->dl_num_bytes_per_cb_) *
      this->dl_ldpc_config_.NumBlocksInSymbol();
  dl_bits_.Calloc(frame_.NumDLSyms(),
                  dl_num_bytes_per_ue_pad * this->ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);
  dl_iq_f_.Calloc(frame_.NumDLSyms(), ofdm_data_num_ * ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);
  dl_iq_t_.Calloc(frame_.NumDLSyms(),
                  this->samps_per_symbol_ * this->ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);

  const size_t ul_num_bytes_per_ue_pad =
      Roundup<64>(this->ul_num_bytes_per_cb_) *
      this->ul_ldpc_config_.NumBlocksInSymbol();
  ul_bits_.Calloc(frame_.NumULSyms(),
                  ul_num_bytes_per_ue_pad * this->ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);
  ul_iq_f_.Calloc(frame_.NumULSyms(),
                  this->ofdm_data_num_ * this->ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);
  ul_iq_t_.Calloc(frame_.NumULSyms(),
                  this->samps_per_symbol_ * this->ue_ant_num_,
                  Agora_memory::Alignment_t::kAlign64);

#ifdef GENERATE_DATA
  for (size_t ue_id = 0; ue_id < this->ue_ant_num_; ue_id++) {
    for (size_t j = 0; j < num_bytes_per_ue_pad; j++) {
      int cur_offset = j * ue_ant_num_ + ue_id;
      for (size_t i = 0; i < frame_.NumULSyms(); i++) {
        this->ul_bits_[i][cur_offset] = rand() % mod_order;
      }
      for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
        this->dl_bits_[i][cur_offset] = rand() % mod_order;
      }
    }
  }
#else
  if (frame_.NumUlDataSyms() > 0) {
    const std::string ul_data_file =
        kUlDataFilePrefix + std::to_string(this->ofdm_ca_num_) + "_ant" +
        std::to_string(this->ue_ant_total_) + ".bin";
    AGORA_LOG_SYMBOL("Config: Reading raw ul data from %s\n",
                     ul_data_file.c_str());
    FILE* fd = std::fopen(ul_data_file.c_str(), "rb");
    if (fd == nullptr) {
      AGORA_LOG_ERROR("Failed to open antenna file %s. Error %s.\n",
                      ul_data_file.c_str(), strerror(errno));
      throw std::runtime_error("Config: Failed to open antenna file");
    }

    for (size_t i = frame_.ClientUlPilotSymbols();
         i < frame_.NumULSyms(); i++) {
      if (std::fseek(fd, (ul_data_bytes_num_persymbol_ * this->ue_ant_offset_),
                     SEEK_CUR) != 0) {
        AGORA_LOG_ERROR(
            " *** Error: failed to seek propertly (pre) into %s file\n",
            ul_data_file.c_str());
        RtAssert(false,
                 "Failed to seek propertly into " + ul_data_file + "file\n");
      }
      for (size_t j = 0; j < this->ue_ant_num_; j++) {
        size_t r = std::fread(this->ul_bits_[i] + (j * ul_num_bytes_per_ue_pad),
                              sizeof(int8_t), ul_data_bytes_num_persymbol_, fd);
        if (r < ul_data_bytes_num_persymbol_) {
          AGORA_LOG_ERROR(
              " *** Error: Uplink bad read from file %s (batch %zu : %zu) "
              "%zu : %zu\n",
              ul_data_file.c_str(), i, j, r, ul_data_bytes_num_persymbol_);
        }
      }
      if (std::fseek(fd,
                     ul_data_bytes_num_persymbol_ *
                         (this->ue_ant_total_ - this->ue_ant_offset_ -
                          this->ue_ant_num_),
                     SEEK_CUR) != 0) {
        AGORA_LOG_ERROR(
            " *** Error: failed to seek propertly (post) into %s file\n",
            ul_data_file.c_str());
        RtAssert(false,
                 "Failed to seek propertly into " + ul_data_file + "file\n");
      }
    }
    std::fclose(fd);
  }

  if (frame_.NumDlDataSyms() > 0) {
    const std::string dl_data_file =
        kDlDataFilePrefix + std::to_string(this->ofdm_ca_num_) + "_ant" +
        std::to_string(this->ue_ant_total_) + ".bin";

    AGORA_LOG_SYMBOL("Config: Reading raw dl data from %s\n",
                     dl_data_file.c_str());
    FILE* fd = std::fopen(dl_data_file.c_str(), "rb");
    if (fd == nullptr) {
      AGORA_LOG_ERROR("Failed to open antenna file %s. Error %s.\n",
                      dl_data_file.c_str(), strerror(errno));
      throw std::runtime_error("Config: Failed to open dl antenna file");
    }

    for (size_t i = this->frame_.ClientDlPilotSymbols();
         i < this->frame_.NumDLSyms(); i++) {
      for (size_t j = 0; j < this->ue_ant_num_; j++) {
        size_t r = std::fread(this->dl_bits_[i] + j * dl_num_bytes_per_ue_pad,
                              sizeof(int8_t), dl_data_bytes_num_persymbol_, fd);
        if (r < dl_data_bytes_num_persymbol_) {
          AGORA_LOG_ERROR(
              "***Error: Downlink bad read from file %s (batch %zu : %zu) "
              "\n",
              dl_data_file.c_str(), i, j);
        }
      }
    }
    std::fclose(fd);
  }
#endif

  auto scrambler = std::make_unique<AgoraScrambler::Scrambler>();

  const size_t ul_encoded_bytes_per_block =
      BitsToBytes(this->ul_ldpc_config_.NumCbCodewLen());
  const size_t ul_num_blocks_per_symbol =
      this->ul_ldpc_config_.NumBlocksInSymbol() * this->ue_ant_num_;

  SimdAlignByteVector ul_scramble_buffer(
      this->ul_num_bytes_per_cb_ + ul_num_padding_bytes_per_cb_, std::byte(0));

  int8_t* ldpc_input = nullptr;
  // Encode uplink bits
  Table<int8_t> ul_encoded_bits;
  ul_encoded_bits.Malloc(frame_.NumULSyms() * ul_num_blocks_per_symbol,
                         ul_encoded_bytes_per_block,
                         Agora_memory::Alignment_t::kAlign64);
  ul_mod_bits_.Calloc(frame_.NumULSyms(),
                      Roundup<64>(this->ofdm_data_num_) * this->ue_ant_num_,
                      Agora_memory::Alignment_t::kAlign32);
  auto* ul_temp_parity_buffer = new int8_t[LdpcEncodingParityBufSize(
      this->ul_ldpc_config_.BaseGraph(),
      this->ul_ldpc_config_.ExpansionFactor())];

  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t j = 0; j < ue_ant_num_; j++) {
      for (size_t k = 0; k < ul_ldpc_config_.NumBlocksInSymbol(); k++) {
        int8_t* coded_bits_ptr =
            ul_encoded_bits[i * ul_num_blocks_per_symbol +
                            j * ul_ldpc_config_.NumBlocksInSymbol() + k];

        if (scramble_enabled_) {
          scrambler->Scramble(
              ul_scramble_buffer.data(),
              GetInfoBits(ul_bits_, Direction::kUplink, i, j, k),
              this->ul_num_bytes_per_cb_);
          ldpc_input = reinterpret_cast<int8_t*>(ul_scramble_buffer.data());
        } else {
          ldpc_input = GetInfoBits(ul_bits_, Direction::kUplink, i, j, k);
        }
        //Clean padding
        if (this->ul_num_bytes_per_cb_ > 0) {
          std::memset(&ldpc_input[this->ul_num_bytes_per_cb_], 0u,
                      ul_num_padding_bytes_per_cb_);
        }
        LdpcEncodeHelper(ul_ldpc_config_.BaseGraph(),
                         ul_ldpc_config_.ExpansionFactor(),
                         ul_ldpc_config_.NumRows(), coded_bits_ptr,
                         ul_temp_parity_buffer, ldpc_input);
        int8_t* mod_input_ptr =
            GetModBitsBuf(ul_mod_bits_, Direction::kUplink, 0, i, j, k);
        AdaptBitsForMod(reinterpret_cast<uint8_t*>(coded_bits_ptr),
                        reinterpret_cast<uint8_t*>(mod_input_ptr),
                        ul_encoded_bytes_per_block, ul_mod_order_bits_);
      }
    }
  }

  // Generate freq-domain uplink symbols
  Table<complex_float> ul_iq_ifft;
  ul_iq_ifft.Calloc(frame_.NumULSyms(),
                    this->ofdm_ca_num_ * this->ue_ant_num_,
                    Agora_memory::Alignment_t::kAlign64);
  std::vector<FILE*> vec_fp_tx;
  if (kOutputUlScData) {
    for (size_t i = 0; i < this->ue_num_; i++) {
      const std::string filename_ul_data_f =
          kUlDataFreqPrefix + ul_modulation_ + "_" +
          std::to_string(ofdm_data_num_) + "_" + std::to_string(ofdm_ca_num_) +
          "_" + std::to_string(kOfdmSymbolPerSlot) + "_" +
          std::to_string(frame_.NumULSyms()) + "_" +
          std::to_string(kOutputFrameNum) + "_" + ue_channel_ + "_" +
          std::to_string(i) + ".bin";
      ul_tx_f_data_files_.push_back(filename_ul_data_f.substr(
          filename_ul_data_f.find_last_of("/\\") + 1));
      FILE* fp_tx_f = std::fopen(filename_ul_data_f.c_str(), "wb");
      if (fp_tx_f == nullptr) {
        AGORA_LOG_ERROR("Failed to create ul sc data file %s. Error %s.\n",
                        filename_ul_data_f.c_str(), strerror(errno));
        throw std::runtime_error("Config: Failed to create ul sc data file");
      }
      vec_fp_tx.push_back(fp_tx_f);
    }
  }
  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t u = 0; u < this->ue_ant_num_; u++) {
      const size_t q = u * ofdm_data_num_;

      for (size_t j = 0; j < ofdm_data_num_; j++) {
        const size_t sc = j + ofdm_data_start_;
        if (i >= frame_.ClientUlPilotSymbols()) {
          int8_t* mod_input_ptr =
              GetModBitsBuf(ul_mod_bits_, Direction::kUplink, 0, i, u, j);
          ul_iq_f_[i][q + j] = ModSingleUint8(*mod_input_ptr, ul_mod_table_);
        } else {
          ul_iq_f_[i][q + j] = ue_specific_pilot_[u][j];
        }
        // FFT Shift
        const size_t k = sc >= ofdm_ca_num_ / 2 ? sc - ofdm_ca_num_ / 2
                                                : sc + ofdm_ca_num_ / 2;
        ul_iq_ifft[i][u * ofdm_ca_num_ + k] = ul_iq_f_[i][q + j];
      }
      if (kOutputUlScData) {
        const auto write_status =
            std::fwrite(&ul_iq_ifft[i][u * ofdm_ca_num_], sizeof(complex_float),
                        ofdm_ca_num_, vec_fp_tx.at(u / num_ue_channels_));
        if (write_status != ofdm_ca_num_) {
          AGORA_LOG_ERROR("Config: Failed to write ul sc data file\n");
        }
      }
      CommsLib::IFFT(&ul_iq_ifft[i][u * ofdm_ca_num_], ofdm_ca_num_, false);
    }
  }
  if (kOutputUlScData) {
    for (size_t i = 0; i < vec_fp_tx.size(); i++) {
      const auto close_status = std::fclose(vec_fp_tx.at(i));
      if (close_status != 0) {
        AGORA_LOG_ERROR("Config: Failed to close ul sc data file %zu\n", i);
      }
    }
  }

  // Encode downlink bits
  const size_t dl_encoded_bytes_per_block =
      BitsToBytes(this->dl_ldpc_config_.NumCbCodewLen());
  const size_t dl_num_blocks_per_symbol =
      this->dl_ldpc_config_.NumBlocksInSymbol() * this->ue_ant_num_;

  SimdAlignByteVector dl_scramble_buffer(
      this->dl_num_bytes_per_cb_ + dl_num_padding_bytes_per_cb_, std::byte(0));

  Table<int8_t> dl_encoded_bits;
  dl_encoded_bits.Malloc(frame_.NumDLSyms() * dl_num_blocks_per_symbol,
                         dl_encoded_bytes_per_block,
                         Agora_memory::Alignment_t::kAlign64);
  dl_mod_bits_.Calloc(frame_.NumDLSyms(),
                      Roundup<64>(GetOFDMDataNum()) * ue_ant_num_,
                      Agora_memory::Alignment_t::kAlign32);
  auto* dl_temp_parity_buffer = new int8_t[LdpcEncodingParityBufSize(
      this->dl_ldpc_config_.BaseGraph(),
      this->dl_ldpc_config_.ExpansionFactor())];

  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t j = 0; j < this->ue_ant_num_; j++) {
      for (size_t k = 0; k < dl_ldpc_config_.NumBlocksInSymbol(); k++) {
        int8_t* coded_bits_ptr =
            dl_encoded_bits[i * dl_num_blocks_per_symbol +
                            j * dl_ldpc_config_.NumBlocksInSymbol() + k];

        if (scramble_enabled_) {
          scrambler->Scramble(
              dl_scramble_buffer.data(),
              GetInfoBits(dl_bits_, Direction::kDownlink, i, j, k),
              this->dl_num_bytes_per_cb_);
          ldpc_input = reinterpret_cast<int8_t*>(dl_scramble_buffer.data());
        } else {
          ldpc_input = GetInfoBits(dl_bits_, Direction::kDownlink, i, j, k);
        }
        if (dl_num_padding_bytes_per_cb_ > 0) {
          std::memset(&ldpc_input[this->dl_num_bytes_per_cb_], 0u,
                      dl_num_padding_bytes_per_cb_);
        }

        LdpcEncodeHelper(dl_ldpc_config_.BaseGraph(),
                         dl_ldpc_config_.ExpansionFactor(),
                         dl_ldpc_config_.NumRows(), coded_bits_ptr,
                         dl_temp_parity_buffer, ldpc_input);
        int8_t* mod_input_ptr =
            GetModBitsBuf(dl_mod_bits_, Direction::kDownlink, 0, i, j, k);
        AdaptBitsForMod(reinterpret_cast<uint8_t*>(coded_bits_ptr),
                        reinterpret_cast<uint8_t*>(mod_input_ptr),
                        dl_encoded_bytes_per_block, dl_mod_order_bits_);
      }
    }
  }

  // Generate freq-domain downlink symbols
  Table<complex_float> dl_iq_ifft;
  dl_iq_ifft.Calloc(frame_.NumDLSyms(), ofdm_ca_num_ * ue_ant_num_,
                    Agora_memory::Alignment_t::kAlign64);
  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t u = 0; u < ue_ant_num_; u++) {
      size_t q = u * ofdm_data_num_;

      for (size_t j = 0; j < ofdm_data_num_; j++) {
        size_t sc = j + ofdm_data_start_;
        if (IsDataSubcarrier(j) == true) {
          int8_t* mod_input_ptr =
              GetModBitsBuf(dl_mod_bits_, Direction::kDownlink, 0, i, u,
                            this->GetOFDMDataIndex(j));
          dl_iq_f_[i][q + j] = ModSingleUint8(*mod_input_ptr, dl_mod_table_);
        } else {
          dl_iq_f_[i][q + j] = ue_specific_pilot_[u][j];
        }
        // FFT Shift
        const size_t k = sc >= ofdm_ca_num_ / 2 ? sc - ofdm_ca_num_ / 2
                                                : sc + ofdm_ca_num_ / 2;
        dl_iq_ifft[i][u * ofdm_ca_num_ + k] = dl_iq_f_[i][q + j];
      }
      CommsLib::IFFT(&dl_iq_ifft[i][u * ofdm_ca_num_], ofdm_ca_num_, false);
    }
  }

  // Find normalization factor through searching for max value in IFFT results
  float ul_max_mag =
      CommsLib::FindMaxAbs(ul_iq_ifft, frame_.NumULSyms(),
                           this->ue_ant_num_ * this->ofdm_ca_num_);
  float dl_max_mag =
      CommsLib::FindMaxAbs(dl_iq_ifft, frame_.NumDLSyms(),
                           this->ue_ant_num_ * this->ofdm_ca_num_);
  float ue_pilot_max_mag = CommsLib::FindMaxAbs(
      ue_pilot_ifft_, this->ue_ant_num_, this->ofdm_ca_num_);
  float pilot_max_mag = CommsLib::FindMaxAbs(pilot_ifft_, this->ofdm_ca_num_);
  // additional 2^2 (6dB) power backoff
  this->scale_ =
      2 * std::max({ul_max_mag, dl_max_mag, ue_pilot_max_mag, pilot_max_mag});

  float dl_papr = dl_max_mag /
                  CommsLib::FindMeanAbs(dl_iq_ifft, frame_.NumDLSyms(),
                                        this->ue_ant_num_ * this->ofdm_ca_num_);
  float ul_papr = ul_max_mag /
                  CommsLib::FindMeanAbs(ul_iq_ifft, frame_.NumULSyms(),
                                        this->ue_ant_num_ * this->ofdm_ca_num_);
  std::printf("Uplink PAPR %2.2f dB, Downlink PAPR %2.2f dB\n",
              10 * std::log10(ul_papr), 10 * std::log10(dl_papr));

  // Generate time domain symbols for downlink
  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t u = 0; u < this->ue_ant_num_; u++) {
      size_t q = u * this->ofdm_ca_num_;
      size_t r = u * this->samps_per_symbol_;
      CommsLib::Ifft2tx(&dl_iq_ifft[i][q], &this->dl_iq_t_[i][r],
                        this->ofdm_ca_num_, this->ofdm_tx_zero_prefix_,
                        this->cp_len_, kDebugDownlink ? 1 : this->scale_);
    }
  }

  // Generate time domain uplink symbols
  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t u = 0; u < this->ue_ant_num_; u++) {
      size_t q = u * this->ofdm_ca_num_;
      size_t r = u * this->samps_per_symbol_;
      CommsLib::Ifft2tx(&ul_iq_ifft[i][q], &ul_iq_t_[i][r], this->ofdm_ca_num_,
                        this->ofdm_tx_zero_prefix_, this->cp_len_,
                        this->scale_);
    }
  }

  // Generate time domain ue-specific pilot symbols
  for (size_t i = 0; i < this->ue_ant_num_; i++) {
    CommsLib::Ifft2tx(ue_pilot_ifft_[i], this->ue_specific_pilot_t_[i],
                      this->ofdm_ca_num_, this->ofdm_tx_zero_prefix_,
                      this->cp_len_, kDebugDownlink ? 1 : this->scale_);
  }

  this->pilot_ci16_.resize(samps_per_symbol_, 0);
  CommsLib::Ifft2tx(pilot_ifft_, this->pilot_ci16_.data(), ofdm_ca_num_,
                    ofdm_tx_zero_prefix_, cp_len_, scale_);

  for (size_t i = 0; i < ofdm_ca_num_; i++) {
    this->pilot_cf32_.emplace_back(pilot_ifft_[i].re / scale_,
                                   pilot_ifft_[i].im / scale_);
  }
  this->pilot_cf32_.insert(this->pilot_cf32_.begin(),
                           this->pilot_cf32_.end() - this->cp_len_,
                           this->pilot_cf32_.end());  // add CP

  // generate a UINT32 version to write to FPGA buffers
  this->pilot_ = Utils::Cfloat32ToUint32(this->pilot_cf32_, false, "QI");

  std::vector<uint32_t> pre_uint32(this->ofdm_tx_zero_prefix_, 0);
  this->pilot_.insert(this->pilot_.begin(), pre_uint32.begin(),
                      pre_uint32.end());
  this->pilot_.resize(this->samps_per_symbol_);

  this->pilot_ue_sc_.resize(ue_ant_num_);
  this->pilot_ue_ci16_.resize(ue_ant_num_);
  for (size_t ue_id = 0; ue_id < this->ue_ant_num_; ue_id++) {
    this->pilot_ue_ci16_.at(ue_id).resize(frame_.NumPilotSyms());
    for (size_t pilot_idx = 0; pilot_idx < frame_.NumPilotSyms();
         pilot_idx++) {
      this->pilot_ue_ci16_.at(ue_id).at(pilot_idx).resize(samps_per_symbol_, 0);
      if (this->freq_orthogonal_pilot_ || ue_id == pilot_idx) {
        std::vector<arma::uword> pilot_sc_list;
        for (size_t sc_id = 0; sc_id < ofdm_data_num_; sc_id++) {
          const size_t org_sc = sc_id + ofdm_data_start_;
          const size_t center_sc = ofdm_ca_num_ / 2;
          // FFT Shift
          const size_t shifted_sc = (org_sc >= center_sc)
                                        ? (org_sc - center_sc)
                                        : (org_sc + center_sc);
          if (this->freq_orthogonal_pilot_ == false ||
              sc_id % this->pilot_sc_group_size_ == ue_id) {
            pilot_ifft_[shifted_sc] = this->pilots_[sc_id];
            pilot_sc_list.push_back(org_sc);
          } else {
            pilot_ifft_[shifted_sc].re = 0.0f;
            pilot_ifft_[shifted_sc].im = 0.0f;
          }
        }
        pilot_ue_sc_.at(ue_id) = arma::uvec(pilot_sc_list);
        CommsLib::IFFT(pilot_ifft_, this->ofdm_ca_num_, false);
        CommsLib::Ifft2tx(pilot_ifft_,
                          this->pilot_ue_ci16_.at(ue_id).at(pilot_idx).data(),
                          ofdm_ca_num_, ofdm_tx_zero_prefix_, cp_len_, scale_);
      }
    }
  }

  if (kDebugPrintPilot) {
    std::cout << "Pilot data = [" << std::endl;
    for (size_t sc_id = 0; sc_id < ofdm_data_num_; sc_id++) {
      std::cout << pilots_[sc_id].re << "+1i*" << pilots_[sc_id].im << " ";
    }
    std::cout << std::endl << "];" << std::endl;
    for (size_t ue_id = 0; ue_id < ue_ant_num_; ue_id++) {
      std::cout << "pilot_ue_sc_[" << ue_id << "] = [" << std::endl
                << pilot_ue_sc_.at(ue_id).as_row() << "];" << std::endl;
      std::cout << "ue_specific_pilot_[" << ue_id << "] = [" << std::endl;
      for (size_t sc_id = 0; sc_id < ofdm_data_num_; sc_id++) {
        std::cout << ue_specific_pilot_[ue_id][sc_id].re << "+1i*"
                  << ue_specific_pilot_[ue_id][sc_id].im << " ";
      }
      std::cout << std::endl << "];" << std::endl;
      std::cout << "ue_pilot_ifft_[" << ue_id << "] = [" << std::endl;
      for (size_t ifft_idx = 0; ifft_idx < ofdm_ca_num_; ifft_idx++) {
        std::cout << ue_pilot_ifft_[ue_id][ifft_idx].re << "+1i*"
                  << ue_pilot_ifft_[ue_id][ifft_idx].im << " ";
      }
      std::cout << std::endl << "];" << std::endl;
    }
  }

  if (pilot_ifft_ != nullptr) {
    FreeBuffer1d(&pilot_ifft_);
  }
  delete[](ul_temp_parity_buffer);
  delete[](dl_temp_parity_buffer);
  ul_iq_ifft.Free();
  dl_iq_ifft.Free();
  dl_encoded_bits.Free();
  ul_encoded_bits.Free();
}

inline size_t SelectZc(size_t base_graph, size_t code_rate,
                       size_t mod_order_bits, size_t num_sc, size_t cb_per_sym,
                       const std::string& dir) {
  size_t n_zc = sizeof(kZc) / sizeof(size_t);
  std::vector<size_t> zc_vec(kZc, kZc + n_zc);
  std::sort(zc_vec.begin(), zc_vec.end());
  // According to cyclic_shift.cc cyclic shifter for zc
  // larger than 256 has not been implemented, so we skip them here.
  size_t max_zc_index =
      (std::find(zc_vec.begin(), zc_vec.end(), kMaxSupportedZc) -
       zc_vec.begin());
  size_t max_uncoded_bits =
      static_cast<size_t>(num_sc * code_rate * mod_order_bits / 1024.0);
  size_t zc = SIZE_MAX;
  size_t i = 0;
  for (; i < max_zc_index; i++) {
    if ((zc_vec.at(i) * LdpcNumInputCols(base_graph) * cb_per_sym <
         max_uncoded_bits) &&
        (zc_vec.at(i + 1) * LdpcNumInputCols(base_graph) * cb_per_sym >
         max_uncoded_bits)) {
      zc = zc_vec.at(i);
      break;
    }
  }
  if (zc == SIZE_MAX) {
    AGORA_LOG_WARN(
        "Exceeded possible range of LDPC lifting Zc for " + dir +
            "! Setting lifting size to max possible value(%zu).\nThis may lead "
            "to too many unused subcarriers. For better use of the PHY "
            "resources, you may reduce your coding or modulation rate.\n",
        kMaxSupportedZc);
    zc = kMaxSupportedZc;
  }
  return zc;
}

LDPCconfig Mcs::Ul_Ldpc_Config() { return ul_ldpc_config_; }

LDPCconfig Mcs::Dl_Ldpc_Config() { return dl_ldpc_config_; }
