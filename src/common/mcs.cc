/**
 * @file mcs.cc
 * @brief Class implementation for mcs handling
 */

#include "mcs.h"

#include <stddef.h>

#include "comms-constants.inc"
#include "comms-lib.h"
#include "data_generator.h"
#include "datatype_conversion.h"
#include "logger.h"
#include "modulation.h"
#include "scrambler.h"
#include "simd_types.h"
#include "symbols.h"
#include "utils_ldpc.h"

static constexpr size_t kMaxSupportedZc = 256;

/// Print the I/Q samples in the pilots
static constexpr bool kDebugPrintPilot = false;

static const std::string kLogFilepath =
    TOSTRING(PROJECT_DIRECTORY) "/files/log/";
static const std::string kExperimentFilepath =
    TOSTRING(PROJECT_DIRECTORY) "/files/experiment/";
static const std::string kUlDataFilePrefix =
    kExperimentFilepath + "LDPC_orig_ul_data_";
static const std::string kDlDataFilePrefix =
    kExperimentFilepath + "LDPC_orig_dl_data_";
static const std::string kUlDataFreqPrefix = kExperimentFilepath + "ul_data_f_";
static constexpr size_t kControlMCS = 5;  // QPSK, 379/1024

Mcs::Mcs(Config* const cfg)
    : pilot_ifft_(nullptr),
      ul_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      dl_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      dl_bcast_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      cfg_(cfg),
      frame_("") {
  pilots_ = nullptr;
  pilots_sgn_ = nullptr;

  this->running_.store(true);

  size_t ofdm_data_num = cfg_->OfdmCaNum();
  ul_mcs_params_ = cfg_->UlMcsParams();
  dl_mcs_params_ = cfg_->DlMcsParams();
  CreateModulationTables();

  this->frame_ = cfg_->Frame();
  ofdm_ca_num_ = cfg_->OfdmCaNum();

  initial_ul_mcs_properties_.base_graph = ul_mcs_params_.value("base_graph", 1);
  initial_ul_mcs_properties_.early_term =
      ul_mcs_params_.value("earlyTermination", true);
  initial_ul_mcs_properties_.max_decoder_iter =
      ul_mcs_params_.value("decoderIter", 5);
  initial_ul_mcs_properties_.ofdm_data_num = ofdm_data_num;

  initial_dl_mcs_properties_.base_graph = dl_mcs_params_.value("base_graph", 1);
  initial_dl_mcs_properties_.early_term =
      dl_mcs_params_.value("earlyTermination", true);
  initial_dl_mcs_properties_.max_decoder_iter =
      dl_mcs_params_.value("decoderIter", 5);

  //Initialize UL MCS
  InitializeUlMcs(ul_mcs_params_);

  //Initialize DL MCS
  InitializeDlMcs(dl_mcs_params_);

  //Update the LDPC
  UpdateUlLdpcConfig();
  UpdateDlLdpcConfig();
  UpdateCtrlMCS();

  CalculateLdpcProperties();

  if (std::filesystem::is_directory(kLogFilepath) == false) {
    std::filesystem::create_directory(kLogFilepath);
  }

  this->DumpMcsInfo();
  this->UpdateCtrlMCS();
}

Mcs::~Mcs() {
  if (pilots_ != nullptr) {
    std::free(pilots_);
    pilots_ = nullptr;
  }
  if (pilots_sgn_ != nullptr) {
    std::free(pilots_sgn_);
    pilots_sgn_ = nullptr;
  }
  ue_specific_pilot_t_.Free();
  ue_specific_pilot_.Free();
  dl_bits_.Free();
  ul_bits_.Free();
  ul_mod_bits_.Free();
  dl_mod_bits_.Free();
  dl_iq_f_.Free();
  dl_iq_t_.Free();
  ul_iq_f_.Free();
  ul_iq_t_.Free();
}

void Mcs::CheckUlMcs(float snr, size_t frame_id) {
  float spectral_efficiency = SpectralEffeciency(snr);
  std::cout << "checking the dl mcs" << std::endl << std::flush;
  if (abs(kmcs_index_to_spectral_effeciency.at(current_dl_mcs_.mcs_index) -
          spectral_efficiency) > 0.1) {
    //Find the closest spectral effeciency to the measured average special
    //effeciency and update the dl mcs to the corresponding dl mcs.

    float min_spectral_effeciency_delta = MAXFLOAT;
    size_t optimal_mcs_index = 0;

    for (auto map_iter = kmcs_index_to_spectral_effeciency.begin();
         map_iter != kmcs_index_to_spectral_effeciency.end(); ++map_iter) {
      if (map_iter->second - spectral_efficiency < 0 &&
          abs(map_iter->second - spectral_efficiency) <
              min_spectral_effeciency_delta) {
        min_spectral_effeciency_delta =
            abs(map_iter->second - spectral_efficiency);
        optimal_mcs_index = map_iter->first;
      }
    }
    SetNextDlMcs(frame_id, GetModOrderBits(optimal_mcs_index));
  }
}

void Mcs::InitializeUlMcs(const nlohmann::json ul_mcs) {
  current_ul_mcs_.frame_number = 0;
  if (ul_mcs.find("mcs_index") == ul_mcs.end()) {
    ul_modulation_ = ul_mcs.value("modulation", "16QAM");
    current_ul_mcs_.mod_order_bits = kModulStringMap.at(ul_modulation_);

    double ul_code_rate_usr = ul_mcs.value("code_rate", 0.333);
    size_t code_rate_int =
        static_cast<size_t>(std::round(ul_code_rate_usr * 1024.0));

    current_ul_mcs_.mcs_index =
        CommsLib::GetMcsIndex(current_ul_mcs_.mod_order_bits, code_rate_int);
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
    current_ul_mcs_.mod_order_bits = GetModOrderBits(current_ul_mcs_.mcs_index);
    current_ul_mcs_.code_rate = GetCodeRate(current_ul_mcs_.mcs_index);
  }
}

void Mcs::InitializeDlMcs(const nlohmann::json dl_mcs) {
  current_dl_mcs_.frame_number = 0;
  if (dl_mcs.find("mcs_index") == dl_mcs.end()) {
    dl_modulation_ = dl_mcs.value("modulation", "16QAM");
    current_dl_mcs_.mod_order_bits = kModulStringMap.at(dl_modulation_);

    double dl_code_rate_usr = dl_mcs.value("code_rate", 0.333);
    size_t code_rate_int =
        static_cast<size_t>(std::round(dl_code_rate_usr * 1024.0));
    current_dl_mcs_.mcs_index =
        CommsLib::GetMcsIndex(current_dl_mcs_.mod_order_bits, code_rate_int);
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
    current_dl_mcs_.mod_order_bits = GetModOrderBits(current_dl_mcs_.mcs_index);
    current_dl_mcs_.code_rate = GetCodeRate(current_dl_mcs_.mcs_index);
  }
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

void Mcs::CreateModulationTables() {
  std::cout << "CREATING MODULATION TABLES" << std::endl << std::flush;

  for (size_t i = 0; i < kNumTables; i++) {
    std::cout << "creating ul table: " << std::endl << std::flush;
    InitModulationTable(modulation_tables_.ul_tables[i], (i + 1) * 2);
    std::cout << modulation_tables_.ul_tables[i][0][0].re << " "
              << modulation_tables_.ul_tables[i][0][0].im << std::endl
              << std::flush;

    std::cout << "creating dl table" << std::endl << std::flush;
    InitModulationTable(modulation_tables_.dl_tables[i], (i + 1) * 2);
    std::cout << "DL table: at 0, 0: "
              << modulation_tables_.dl_tables[i][0][0].re << " "
              << modulation_tables_.dl_tables[i][0][0].im << std::endl
              << std::flush;
  }
}

void Mcs::UpdateMcs(size_t current_frame_number) {
  UpdateUlMcs(current_frame_number);
  UpdateDlMcs(current_frame_number);

  CalculateLdpcProperties();
}

void Mcs::UpdateUlMcs(size_t current_frame_number) {
  if (current_frame_number >= next_ul_mcs_.frame_number) {
    current_ul_mcs_.frame_number = current_frame_number;
    current_ul_mcs_.mcs_index = next_ul_mcs_.mcs_index;
  }
  UpdateUlLdpcConfig();
}

void Mcs::UpdateDlMcs(size_t current_frame_number) {
  if (current_frame_number >= next_dl_mcs_.frame_number) {
    current_dl_mcs_.frame_number = current_frame_number;
    current_dl_mcs_.mcs_index = next_dl_mcs_.mcs_index;
  }
  UpdateDlLdpcConfig();
}

void Mcs::SetNextUlMcs(size_t frame_number, size_t mod_order_bits) {
  next_ul_mcs_.frame_number = frame_number;
  next_ul_mcs_.mcs_index = mod_order_bits;
}

void Mcs::SetNextDlMcs(size_t frame_number, size_t mod_order_bits) {
  next_dl_mcs_.frame_number = frame_number;
  next_dl_mcs_.mcs_index = mod_order_bits;
}

void Mcs::UpdateUlLdpcConfig() {
  uint16_t base_graph = initial_ul_mcs_properties_.base_graph;

  size_t ul_mod_order_bits = GetModOrderBits(current_ul_mcs_.mcs_index);
  size_t ul_code_rate = GetCodeRate(current_ul_mcs_.mcs_index);

  size_t zc = SelectZc(base_graph, ul_code_rate, ul_mod_order_bits,
                       cfg_->OfdmDataNum(), this->kCbPerSymbol, "uplink");

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
  RtAssert(
      (frame_.NumULSyms() == 0) || (ul_ldpc_config_.NumBlocksInSymbol() > 0),
      "Uplink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");
}

void Mcs::UpdateDlLdpcConfig() {
  uint16_t base_graph = initial_dl_mcs_properties_.base_graph;

  size_t dl_mod_order_bits = GetModOrderBits(current_dl_mcs_.mcs_index);
  size_t dl_code_rate = GetCodeRate(current_dl_mcs_.mcs_index);

  size_t zc = SelectZc(base_graph, dl_code_rate, dl_mod_order_bits,
                       cfg_->GetOFDMDataNum(), this->kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows = static_cast<size_t>(std::round(
                        1024.0 * LdpcNumInputCols(base_graph) / dl_code_rate)) -
                    (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  dl_ldpc_config_ =
      LDPCconfig(base_graph, zc, initial_dl_mcs_properties_.max_decoder_iter,
                 initial_dl_mcs_properties_.early_term, num_cb_len,
                 num_cb_codew_len, num_rows, 0);

  dl_ldpc_config_.NumBlocksInSymbol((cfg_->GetOFDMDataNum() * dl_mod_order_bits) /
                                    dl_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumDLSyms() == 0) || (dl_ldpc_config_.NumBlocksInSymbol() > 0),
      "Downlink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");
}

void Mcs::CalculateLdpcProperties() {
  this->ul_num_bytes_per_cb_ = ul_ldpc_config_.NumCbLen() / 8;
  this->ul_num_padding_bytes_per_cb_ =
      Roundup<64>(ul_num_bytes_per_cb_) - ul_num_bytes_per_cb_;
  this->ul_data_bytes_num_persymbol_ =
      ul_num_bytes_per_cb_ * ul_ldpc_config_.NumBlocksInSymbol();
  ul_mac_packet_length_ = ul_data_bytes_num_persymbol_;

  //((cb_len_bits / zc_size) - 1) * (zc_size / 8) + kProcBytes(32)
  // std::cout << "ul_ldpc_config_.NumCbLen: "
  //           << std::to_string(ul_ldpc_config_.NumCbLen()) << std::endl
  //           << std::flush;
  // std::cout << "ul_ldpc_config_.ExpansionFactor: "
  //           << std::to_string(ul_ldpc_config_.ExpansionFactor()) << std::endl
  //           << std::flush;

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
        ul_num_padding_bytes_per_cb_, increased_padding, ul_num_bytes_per_cb_,
        ul_ldpc_input_min);
    ul_num_padding_bytes_per_cb_ = increased_padding;
  }

  // Smallest over the air packet structure
  RtAssert(frame_.NumULSyms() == 0 ||
               ul_mac_packet_length_ > sizeof(MacPacketHeaderPacked),
           "Uplink MAC Packet size must be larger than MAC header size");
  this->ul_mac_data_length_max_ =
      ul_mac_packet_length_ - sizeof(MacPacketHeaderPacked);

  this->ul_mac_packets_perframe_ = frame_.NumUlDataSyms();
  this->ul_mac_data_bytes_num_perframe_ =
      this->ul_mac_data_length_max_ * this->ul_mac_packets_perframe_;
  this->ul_mac_bytes_num_perframe_ =
      ul_mac_packet_length_ * this->ul_mac_packets_perframe_;

  this->dl_num_bytes_per_cb_ = dl_ldpc_config_.NumCbLen() / 8;
  this->dl_num_padding_bytes_per_cb_ =
      Roundup<64>(this->dl_num_bytes_per_cb_) - this->dl_num_bytes_per_cb_;
  this->dl_data_bytes_num_persymbol_ =
      this->dl_num_bytes_per_cb_ * dl_ldpc_config_.NumBlocksInSymbol();
  this->dl_mac_packet_length_ = this->dl_data_bytes_num_persymbol_;
  // Smallest over the air packet structure
  RtAssert(frame_.NumDLSyms() == 0 ||
               this->dl_mac_packet_length_ > sizeof(MacPacketHeaderPacked),
           "Downlink MAC Packet size must be larger than MAC header size");
  this->dl_mac_data_length_max_ =
      this->dl_mac_packet_length_ - sizeof(MacPacketHeaderPacked);

  this->dl_mac_packets_perframe_ = frame_.NumDlDataSyms();
  this->dl_mac_data_bytes_num_perframe_ =
      this->dl_mac_data_length_max_ * this->dl_mac_packets_perframe_;
  this->dl_mac_bytes_num_perframe_ =
      this->dl_mac_packet_length_ * this->dl_mac_packets_perframe_;

  //((cb_len_bits / zc_size) - 1) * (zc_size / 8) + kProcBytes(32)
  const size_t dl_ldpc_input_min =
      (((dl_ldpc_config_.NumCbLen() / dl_ldpc_config_.ExpansionFactor()) - 1) *
           (dl_ldpc_config_.ExpansionFactor() / 8) +
       32);
  const size_t dl_ldpc_sugg_input = LdpcEncodingInputBufSize(
      dl_ldpc_config_.BaseGraph(), dl_ldpc_config_.ExpansionFactor());

  if (dl_ldpc_input_min >
      (this->dl_num_bytes_per_cb_ + this->dl_num_padding_bytes_per_cb_)) {
    // Can cause a lot of wasted space, specifically the second argument of the max
    const size_t increased_padding =
        Roundup<64>(dl_ldpc_sugg_input) - this->dl_num_bytes_per_cb_;

    AGORA_LOG_WARN(
        "LDPC required Input Buffer size exceeds downlink code block size!, "
        "Increased cb padding from %zu to %zu Downlink CB Bytes %zu, LDPC "
        "Input Min for zc 64:256: %zu\n",
        this->dl_num_padding_bytes_per_cb_, increased_padding,
        this->dl_num_bytes_per_cb_, dl_ldpc_input_min);
    this->dl_num_padding_bytes_per_cb_ = increased_padding;
  }
}

void Mcs::GenPilots() {
  if ((kUseArgos == true) || (kUseUHD == true) || (kUsePureUHD == true)) {
    std::vector<std::vector<double>> gold_ifft =
        CommsLib::GetSequence(128, CommsLib::kGoldIfft);
    std::vector<std::complex<int16_t>> gold_ifft_ci16 =
        Utils::DoubleToCint16(gold_ifft);
    for (size_t i = 0; i < 128; i++) {
      this->gold_cf32_.emplace_back(gold_ifft[0][i], gold_ifft[1][i]);
    }

    std::vector<std::vector<double>> sts_seq =
        CommsLib::GetSequence(0, CommsLib::kStsSeq);
    std::vector<std::complex<int16_t>> sts_seq_ci16 =
        Utils::DoubleToCint16(sts_seq);

    // Populate STS (stsReps repetitions)
    int sts_reps = 15;
    for (int i = 0; i < sts_reps; i++) {
      this->beacon_ci16_.insert(this->beacon_ci16_.end(), sts_seq_ci16.begin(),
                                sts_seq_ci16.end());
    }

    // Populate gold sequence (two reps, 128 each)
    int gold_reps = 2;
    for (int i = 0; i < gold_reps; i++) {
      this->beacon_ci16_.insert(this->beacon_ci16_.end(),
                                gold_ifft_ci16.begin(), gold_ifft_ci16.end());
    }

    this->beacon_len_ = this->beacon_ci16_.size();

    if (cfg_->SampsPerSymbol() < (this->beacon_len_ + cfg_->OfdmTxZeroPrefix() +
                                  cfg_->OfdmTxZeroPostfix())) {
      std::string msg = "Minimum supported symbol_size is ";
      msg += std::to_string(this->beacon_len_);
      throw std::invalid_argument(msg);
    }

    this->beacon_ = Utils::Cint16ToUint32(this->beacon_ci16_, false, "QI");
    this->coeffs_ = Utils::Cint16ToUint32(gold_ifft_ci16, true, "QI");

    // Add addition padding for beacon sent from host
    int frac_beacon = cfg_->SampsPerSymbol() % this->beacon_len_;
    std::vector<std::complex<int16_t>> pre_beacon(cfg_->OfdmTxZeroPrefix(), 0);
    std::vector<std::complex<int16_t>> post_beacon(
        cfg_->OfdmTxZeroPostfix() + frac_beacon, 0);
    this->beacon_ci16_.insert(this->beacon_ci16_.begin(), pre_beacon.begin(),
                              pre_beacon.end());
    this->beacon_ci16_.insert(this->beacon_ci16_.end(), post_beacon.begin(),
                              post_beacon.end());
  }

  // Generate common pilots based on Zadoff-Chu sequence for channel estimation
  auto zc_seq_double =
      CommsLib::GetSequence(cfg_->OfdmDataNum(), CommsLib::kLteZadoffChu);
  auto zc_seq = Utils::DoubleToCfloat(zc_seq_double);
  this->common_pilot_ =
      CommsLib::SeqCyclicShift(zc_seq, M_PI / 4);  // Used in LTE SRS

  pilots_ = static_cast<complex_float*>(Agora_memory::PaddedAlignedAlloc(
      Agora_memory::Alignment_t::kAlign64,
      cfg_->OfdmDataNum() * sizeof(complex_float)));
  pilots_sgn_ = static_cast<complex_float*>(Agora_memory::PaddedAlignedAlloc(
      Agora_memory::Alignment_t::kAlign64,
      cfg_->OfdmDataNum() * sizeof(complex_float)));  // used in CSI estimation
  for (size_t i = 0; i < cfg_->OfdmDataNum(); i++) {
    pilots_[i] = {this->common_pilot_[i].real(), this->common_pilot_[i].imag()};
    auto pilot_sgn = this->common_pilot_[i] /
                     (float)std::pow(std::abs(this->common_pilot_[i]), 2);
    pilots_sgn_[i] = {pilot_sgn.real(), pilot_sgn.imag()};
  }

  RtAssert(pilot_ifft_ == nullptr, "pilot_ifft_ should be null");
  AllocBuffer1d(&pilot_ifft_, this->ofdm_ca_num_,
                Agora_memory::Alignment_t::kAlign64, 1);

  for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
    // FFT Shift
    const size_t k = j + cfg_->OfdmDataStart() >= ofdm_ca_num_ / 2
                         ? j + cfg_->OfdmDataStart() - ofdm_ca_num_ / 2
                         : j + cfg_->OfdmDataStart() + ofdm_ca_num_ / 2;
    pilot_ifft_[k] = pilots_[j];
  }
  CommsLib::IFFT(pilot_ifft_, this->ofdm_ca_num_, false);

  // Generate UE-specific pilots based on Zadoff-Chu sequence for phase tracking
  this->ue_specific_pilot_.Malloc(cfg_->UeAntNum(), cfg_->OfdmDataNum(),
                                  Agora_memory::Alignment_t::kAlign64);
  this->ue_specific_pilot_t_.Calloc(cfg_->UeAntNum(), cfg_->SampsPerSymbol(),
                                    Agora_memory::Alignment_t::kAlign64);

  ue_pilot_ifft_.Calloc(cfg_->UeAntNum(), this->ofdm_ca_num_,
                        Agora_memory::Alignment_t::kAlign64);
  for (size_t i = 0; i < cfg_->UeAntNum(); i++) {
    auto zc_ue_pilot_i = CommsLib::SeqCyclicShift(
        zc_seq,
        (i + cfg_->UeAntOffset()) * (float)M_PI / 6);  // LTE DMRS
    for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
      this->ue_specific_pilot_[i][j] = {zc_ue_pilot_i[j].real(),
                                        zc_ue_pilot_i[j].imag()};
      // FFT Shift
      const size_t k = j + cfg_->OfdmDataStart() >= ofdm_ca_num_ / 2
                           ? j + cfg_->OfdmDataStart() - ofdm_ca_num_ / 2
                           : j + cfg_->OfdmDataStart() + ofdm_ca_num_ / 2;
      ue_pilot_ifft_[i][k] = this->ue_specific_pilot_[i][j];
    }
    CommsLib::IFFT(ue_pilot_ifft_[i], ofdm_ca_num_, false);
  }
}

void Mcs::GenData() {
  this->GenPilots();
  // Get uplink and downlink raw bits either from file or random numbers
  const size_t dl_num_bytes_per_ue_pad =
      Roundup<64>(this->dl_num_bytes_per_cb_) *
      this->dl_ldpc_config_.NumBlocksInSymbol();
  dl_bits_.Calloc(frame_.NumDLSyms(),
                  dl_num_bytes_per_ue_pad * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);
  dl_iq_f_.Calloc(frame_.NumDLSyms(), cfg_->OfdmDataNum() * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);
  dl_iq_t_.Calloc(frame_.NumDLSyms(), cfg_->SampsPerSymbol() * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);

  const size_t ul_num_bytes_per_ue_pad =
      Roundup<64>(this->ul_num_bytes_per_cb_) *
      this->ul_ldpc_config_.NumBlocksInSymbol();
  ul_bits_.Calloc(frame_.NumULSyms(),
                  ul_num_bytes_per_ue_pad * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);
  ul_iq_f_.Calloc(frame_.NumULSyms(), cfg_->OfdmDataNum() * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);
  ul_iq_t_.Calloc(frame_.NumULSyms(), cfg_->SampsPerSymbol() * cfg_->UeAntNum(),
                  Agora_memory::Alignment_t::kAlign64);

#ifdef GENERATE_DATA
  for (size_t ue_id = 0; ue_id < cfg_->UeAntNum(); ue_id++) {
    for (size_t j = 0; j < num_bytes_per_ue_pad; j++) {
      int cur_offset = j * cfg_->UeAntNum() + ue_id;
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
        std::to_string(cfg_->UeAntTotal()) + ".bin";
    AGORA_LOG_SYMBOL("mac_sched->Cfg(): Reading raw ul data from %s\n",
                     ul_data_file.c_str());
    FILE* fd = std::fopen(ul_data_file.c_str(), "rb");
    if (fd == nullptr) {
      // AGORA_LOG_ERROR("Failed to open antenna file %s. Error %s.\n",
      //                 ul_data_file.c_str(), strerror(errno));
      throw std::runtime_error("mac_sched->Cfg(): Failed to open antenna file");
    }

    for (size_t i = frame_.ClientUlPilotSymbols(); i < frame_.NumULSyms();
         i++) {
      if (std::fseek(fd, (ul_data_bytes_num_persymbol_ * cfg_->UeAntOffset()),
                     SEEK_CUR) != 0) {
        AGORA_LOG_ERROR(
            std::to_string(cfg_->UeAntTotal()),
            " *** Error: failed to seek propertly (pre) into %s file\n",
            ul_data_file.c_str());
        RtAssert(false,
                 "Failed to seek propertly into " + ul_data_file + "file\n");
      }
      for (size_t j = 0; j < cfg_->UeAntNum(); j++) {
        size_t r = std::fread(this->ul_bits_[i] + (j * ul_num_bytes_per_ue_pad),
                              sizeof(int8_t), ul_data_bytes_num_persymbol_, fd);
        if (r < ul_data_bytes_num_persymbol_) {
          AGORA_LOG_ERROR(
              " *** Error: Uplink bad read from file %s (batch %zu : %zu) "
              "%zu : %zu\n",
              ul_data_file.c_str(), i, j, r, ul_data_bytes_num_persymbol_);
        }
      }
      if (std::fseek(
              fd,
              ul_data_bytes_num_persymbol_ *
                  (cfg_->UeAntTotal() - cfg_->UeAntOffset() - cfg_->UeAntNum()),
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
        std::to_string(cfg_->UeAntTotal()) + ".bin";

    AGORA_LOG_SYMBOL("mac_sched->Cfg(): Reading raw dl data from %s\n",
                     dl_data_file.c_str());
    FILE* fd = std::fopen(dl_data_file.c_str(), "rb");
    if (fd == nullptr) {
      AGORA_LOG_ERROR("Failed to open antenna file %s. Error %s.\n",
                      dl_data_file.c_str(), strerror(errno));
      throw std::runtime_error(
          "mac_sched->Cfg(): Failed to open dl antenna file");
    }

    for (size_t i = this->frame_.ClientDlPilotSymbols();
         i < this->frame_.NumDLSyms(); i++) {
      for (size_t j = 0; j < cfg_->UeAntNum(); j++) {
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
      this->ul_ldpc_config_.NumBlocksInSymbol() * cfg_->UeAntNum();

  SimdAlignByteVector ul_scramble_buffer(
      this->ul_num_bytes_per_cb_ + ul_num_padding_bytes_per_cb_, std::byte(0));

  int8_t* ldpc_input = nullptr;
  // Encode uplink bits
  Table<int8_t> ul_encoded_bits;
  ul_encoded_bits.Malloc(frame_.NumULSyms() * ul_num_blocks_per_symbol,
                         ul_encoded_bytes_per_block,
                         Agora_memory::Alignment_t::kAlign64);
  ul_mod_bits_.Calloc(frame_.NumULSyms(),
                      Roundup<64>(cfg_->OfdmDataNum()) * cfg_->UeAntNum(),
                      Agora_memory::Alignment_t::kAlign32);
  auto* ul_temp_parity_buffer = new int8_t[LdpcEncodingParityBufSize(
      this->ul_ldpc_config_.BaseGraph(),
      this->ul_ldpc_config_.ExpansionFactor())];

  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t j = 0; j < cfg_->UeAntNum(); j++) {
      for (size_t k = 0; k < ul_ldpc_config_.NumBlocksInSymbol(); k++) {
        int8_t* coded_bits_ptr =
            ul_encoded_bits[i * ul_num_blocks_per_symbol +
                            j * ul_ldpc_config_.NumBlocksInSymbol() + k];

        if (cfg_->ScrambleEnabled()) {
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
                        ul_encoded_bytes_per_block,
                        this->current_ul_mcs_.mod_order_bits);
      }
    }
  }

  // Generate freq-domain uplink symbols
  Table<complex_float> ul_iq_ifft;
  ul_iq_ifft.Calloc(frame_.NumULSyms(), this->ofdm_ca_num_ * cfg_->UeAntNum(),
                    Agora_memory::Alignment_t::kAlign64);
  std::vector<FILE*> vec_fp_tx;
  if (kOutputUlScData) {
    for (size_t i = 0; i < cfg_->UeNum(); i++) {
      const std::string filename_ul_data_f =
          kUlDataFreqPrefix + ul_modulation_ + "_" +
          std::to_string(cfg_->OfdmDataNum()) + "_" +
          std::to_string(ofdm_ca_num_) + "_" +
          std::to_string(kOfdmSymbolPerSlot) + "_" +
          std::to_string(frame_.NumULSyms()) + "_" +
          std::to_string(kOutputFrameNum) + "_" + cfg_->UeChannel() + "_" +
          std::to_string(i) + ".bin";
      ul_tx_f_data_files_.push_back(filename_ul_data_f.substr(
          filename_ul_data_f.find_last_of("/\\") + 1));
      FILE* fp_tx_f = std::fopen(filename_ul_data_f.c_str(), "wb");
      if (fp_tx_f == nullptr) {
        AGORA_LOG_ERROR("Failed to create ul sc data file %s. Error %s.\n",
                        filename_ul_data_f.c_str(), strerror(errno));
        throw std::runtime_error(
            "mac_sched->Cfg(): Failed to create ul sc data file");
      }
      vec_fp_tx.push_back(fp_tx_f);
    }
  }
  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t u = 0; u < cfg_->UeAntNum(); u++) {
      const size_t q = u * cfg_->OfdmDataNum();

      for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
        const size_t sc = j + cfg_->OfdmDataStart();
        if (i >= frame_.ClientUlPilotSymbols()) {
          int8_t* mod_input_ptr =
              GetModBitsBuf(ul_mod_bits_, Direction::kUplink, 0, i, u, j);
          ul_iq_f_[i][q + j] = ModSingleUint8(
              *mod_input_ptr,
              modulation_tables_
                  .ul_tables[current_ul_mcs_.mod_order_bits / 2 - 1]);
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
                        ofdm_ca_num_, vec_fp_tx.at(u / cfg_->NumUeChannels()));
        if (write_status != ofdm_ca_num_) {
          AGORA_LOG_ERROR(
              "mac_sched->Cfg(): Failed to write ul sc data file\n");
        }
      }
      CommsLib::IFFT(&ul_iq_ifft[i][u * ofdm_ca_num_], ofdm_ca_num_, false);
    }
  }
  if (kOutputUlScData) {
    for (size_t i = 0; i < vec_fp_tx.size(); i++) {
      const auto close_status = std::fclose(vec_fp_tx.at(i));
      if (close_status != 0) {
        AGORA_LOG_ERROR(
            "mac_sched->Cfg(): Failed to close ul sc data file %zu\n", i);
      }
    }
  }

  // Encode downlink bits
  const size_t dl_encoded_bytes_per_block =
      BitsToBytes(this->dl_ldpc_config_.NumCbCodewLen());
  const size_t dl_num_blocks_per_symbol =
      this->dl_ldpc_config_.NumBlocksInSymbol() * cfg_->UeAntNum();

  SimdAlignByteVector dl_scramble_buffer(
      this->dl_num_bytes_per_cb_ + dl_num_padding_bytes_per_cb_, std::byte(0));

  Table<int8_t> dl_encoded_bits;
  dl_encoded_bits.Malloc(frame_.NumDLSyms() * dl_num_blocks_per_symbol,
                         dl_encoded_bytes_per_block,
                         Agora_memory::Alignment_t::kAlign64);
  dl_mod_bits_.Calloc(frame_.NumDLSyms(),
                      Roundup<64>(cfg_->GetOFDMDataNum()) * cfg_->UeAntNum(),
                      Agora_memory::Alignment_t::kAlign32);
  auto* dl_temp_parity_buffer = new int8_t[LdpcEncodingParityBufSize(
      this->dl_ldpc_config_.BaseGraph(),
      this->dl_ldpc_config_.ExpansionFactor())];

  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t j = 0; j < cfg_->UeAntNum(); j++) {
      for (size_t k = 0; k < dl_ldpc_config_.NumBlocksInSymbol(); k++) {
        int8_t* coded_bits_ptr =
            dl_encoded_bits[i * dl_num_blocks_per_symbol +
                            j * dl_ldpc_config_.NumBlocksInSymbol() + k];

        if (cfg_->ScrambleEnabled()) {
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
                        dl_encoded_bytes_per_block,
                        this->current_dl_mcs_.mod_order_bits);
      }
    }
  }

  // Generate freq-domain downlink symbols
  Table<complex_float> dl_iq_ifft;
  dl_iq_ifft.Calloc(frame_.NumDLSyms(), ofdm_ca_num_ * cfg_->UeAntNum(),
                    Agora_memory::Alignment_t::kAlign64);
  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t u = 0; u < cfg_->UeAntNum(); u++) {
      size_t q = u * cfg_->OfdmDataNum();

      for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
        size_t sc = j + cfg_->OfdmDataStart();
        if (cfg_->IsDataSubcarrier(j) == true) {
          int8_t* mod_input_ptr =
              GetModBitsBuf(dl_mod_bits_, Direction::kDownlink, 0, i, u,
                            cfg_->GetOFDMDataIndex(j));
          dl_iq_f_[i][q + j] = ModSingleUint8(
              *mod_input_ptr,
              modulation_tables_
                  .dl_tables[current_dl_mcs_.mod_order_bits / 2 - 1]);
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
  float ul_max_mag = CommsLib::FindMaxAbs(
      ul_iq_ifft, frame_.NumULSyms(), cfg_->UeAntNum() * this->ofdm_ca_num_);
  float dl_max_mag = CommsLib::FindMaxAbs(
      dl_iq_ifft, frame_.NumDLSyms(), cfg_->UeAntNum() * this->ofdm_ca_num_);
  float ue_pilot_max_mag = CommsLib::FindMaxAbs(
      ue_pilot_ifft_, cfg_->UeAntNum(), this->ofdm_ca_num_);
  float pilot_max_mag = CommsLib::FindMaxAbs(pilot_ifft_, this->ofdm_ca_num_);
  // additional 2^2 (6dB) power backoff
  this->scale_ =
      2 * std::max({ul_max_mag, dl_max_mag, ue_pilot_max_mag, pilot_max_mag});

  float dl_papr =
      dl_max_mag / CommsLib::FindMeanAbs(dl_iq_ifft, frame_.NumDLSyms(),
                                         cfg_->UeAntNum() * this->ofdm_ca_num_);
  float ul_papr =
      ul_max_mag / CommsLib::FindMeanAbs(ul_iq_ifft, frame_.NumULSyms(),
                                         cfg_->UeAntNum() * this->ofdm_ca_num_);
  std::printf("Uplink PAPR %2.2f dB, Downlink PAPR %2.2f dB\n",
              10 * std::log10(ul_papr), 10 * std::log10(dl_papr));

  // Generate time domain symbols for downlink
  for (size_t i = 0; i < frame_.NumDLSyms(); i++) {
    for (size_t u = 0; u < cfg_->UeAntNum(); u++) {
      size_t q = u * this->ofdm_ca_num_;
      size_t r = u * cfg_->SampsPerSymbol();
      CommsLib::Ifft2tx(&dl_iq_ifft[i][q], &this->dl_iq_t_[i][r],
                        this->ofdm_ca_num_, cfg_->OfdmTxZeroPrefix(),
                        cfg_->CpLen(), kDebugDownlink ? 1 : this->scale_);
    }
  }

  // Generate time domain uplink symbols
  for (size_t i = 0; i < frame_.NumULSyms(); i++) {
    for (size_t u = 0; u < cfg_->UeAntNum(); u++) {
      size_t q = u * this->ofdm_ca_num_;
      size_t r = u * cfg_->SampsPerSymbol();
      CommsLib::Ifft2tx(&ul_iq_ifft[i][q], &ul_iq_t_[i][r], this->ofdm_ca_num_,
                        cfg_->OfdmTxZeroPrefix(), cfg_->CpLen(), this->scale_);
    }
  }

  // Generate time domain ue-specific pilot symbols
  for (size_t i = 0; i < cfg_->UeAntNum(); i++) {
    CommsLib::Ifft2tx(ue_pilot_ifft_[i], this->ue_specific_pilot_t_[i],
                      this->ofdm_ca_num_, cfg_->OfdmTxZeroPrefix(),
                      cfg_->CpLen(), kDebugDownlink ? 1 : this->scale_);
  }

  this->pilot_ci16_.resize(cfg_->SampsPerSymbol(), 0);
  CommsLib::Ifft2tx(pilot_ifft_, this->pilot_ci16_.data(), ofdm_ca_num_,
                    cfg_->OfdmTxZeroPrefix(), cfg_->CpLen(), this->scale_);

  for (size_t i = 0; i < ofdm_ca_num_; i++) {
    this->pilot_cf32_.emplace_back(pilot_ifft_[i].re / this->scale_,
                                   pilot_ifft_[i].im / this->scale_);
  }
  this->pilot_cf32_.insert(this->pilot_cf32_.begin(),
                           this->pilot_cf32_.end() - cfg_->CpLen(),
                           this->pilot_cf32_.end());  // add CP

  // generate a UINT32 version to write to FPGA buffers
  this->pilot_ = Utils::Cfloat32ToUint32(this->pilot_cf32_, false, "QI");

  std::vector<uint32_t> pre_uint32(cfg_->OfdmTxZeroPrefix(), 0);
  this->pilot_.insert(this->pilot_.begin(), pre_uint32.begin(),
                      pre_uint32.end());
  this->pilot_.resize(cfg_->SampsPerSymbol());

  this->pilot_ue_sc_.resize(cfg_->UeAntNum());
  this->pilot_ue_ci16_.resize(cfg_->UeAntNum());
  for (size_t ue_id = 0; ue_id < cfg_->UeAntNum(); ue_id++) {
    this->pilot_ue_ci16_.at(ue_id).resize(frame_.NumPilotSyms());
    for (size_t pilot_idx = 0; pilot_idx < frame_.NumPilotSyms(); pilot_idx++) {
      this->pilot_ue_ci16_.at(ue_id).at(pilot_idx).resize(
          cfg_->SampsPerSymbol(), 0);
      if (cfg_->FreqOrthogonalPilot() || ue_id == pilot_idx) {
        std::vector<arma::uword> pilot_sc_list;
        for (size_t sc_id = 0; sc_id < cfg_->OfdmDataNum(); sc_id++) {
          const size_t org_sc = sc_id + cfg_->OfdmDataStart();
          const size_t center_sc = ofdm_ca_num_ / 2;
          // FFT Shift
          const size_t shifted_sc = (org_sc >= center_sc)
                                        ? (org_sc - center_sc)
                                        : (org_sc + center_sc);
          if (cfg_->FreqOrthogonalPilot() == false ||
              sc_id % cfg_->PilotScGroupSize() == ue_id) {
            pilot_ifft_[shifted_sc] = pilots_[sc_id];
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
                          ofdm_ca_num_, cfg_->OfdmTxZeroPrefix(), cfg_->CpLen(),
                          this->scale_);
      }
    }
  }

  if (kDebugPrintPilot) {
    std::cout << "Pilot data = [" << std::endl;
    for (size_t sc_id = 0; sc_id < cfg_->OfdmDataNum(); sc_id++) {
      std::cout << pilots_[sc_id].re << "+1i*" << pilots_[sc_id].im << " ";
    }
    std::cout << std::endl << "];" << std::endl;
    for (size_t ue_id = 0; ue_id < cfg_->UeAntNum(); ue_id++) {
      std::cout << "pilot_ue_sc_[" << ue_id << "] = [" << std::endl
                << pilot_ue_sc_.at(ue_id).as_row() << "];" << std::endl;
      std::cout << "ue_specific_pilot_[" << ue_id << "] = [" << std::endl;
      for (size_t sc_id = 0; sc_id < cfg_->OfdmDataNum(); sc_id++) {
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
  ue_pilot_ifft_.Free();
}

size_t Mcs::DecodeBroadcastSlots(const int16_t* const bcast_iq_samps) {
  size_t start_tsc = GetTime::WorkerRdtsc();
  size_t delay_offset = (cfg_->OfdmRxZeroPrefixClient() + cfg_->CpLen()) * 2;
  complex_float* bcast_fft_buff = static_cast<complex_float*>(
      Agora_memory::PaddedAlignedAlloc(Agora_memory::Alignment_t::kAlign64,
                                       cfg_->OfdmCaNum() * sizeof(float) * 2));
  SimdConvertShortToFloat(&bcast_iq_samps[delay_offset],
                          reinterpret_cast<float*>(bcast_fft_buff),
                          cfg_->OfdmCaNum() * 2);
  CommsLib::FFT(bcast_fft_buff, cfg_->OfdmCaNum());
  CommsLib::FFTShift(bcast_fft_buff, cfg_->OfdmCaNum());
  auto* bcast_buff_complex = reinterpret_cast<arma::cx_float*>(bcast_fft_buff);

  const size_t sc_num = cfg_->GetOFDMCtrlNum();
  const size_t ctrl_sc_num =
      dl_bcast_ldpc_config_.NumCbCodewLen() / dl_bcast_mod_order_bits_;
  std::vector<arma::cx_float> csi_buff(cfg_->OfdmDataNum());
  arma::cx_float* eq_buff =
      static_cast<arma::cx_float*>(Agora_memory::PaddedAlignedAlloc(
          Agora_memory::Alignment_t::kAlign64, sc_num * sizeof(float) * 2));

  // estimate channel from pilot subcarriers
  float phase_shift = 0;
  for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
    size_t sc_id = j + cfg_->OfdmDataStart();
    complex_float p = pilots_[j];
    if (j % cfg_->OfdmPilotSpacing() == 0) {
      csi_buff.at(j) = (bcast_buff_complex[sc_id] / arma::cx_float(p.re, p.im));
    } else {
      ///\todo not correct when 0th subcarrier is not pilot
      csi_buff.at(j) = csi_buff.at(j - 1);
      if (j % cfg_->OfdmPilotSpacing() == 1) {
        phase_shift += arg((bcast_buff_complex[sc_id] / csi_buff.at(j)) *
                           arma::cx_float(p.re, -p.im));
      }
    }
  }
  phase_shift /= cfg_->GetOFDMPilotNum();
  for (size_t j = 0; j < cfg_->OfdmDataNum(); j++) {
    size_t sc_id = j + cfg_->OfdmDataStart();
    if (cfg_->IsControlSubcarrier(j) == true) {
      eq_buff[cfg_->GetOFDMCtrlIndex(j)] =
          (bcast_buff_complex[sc_id] / csi_buff.at(j)) *
          exp(arma::cx_float(0, -phase_shift));
    }
  }
  int8_t* demod_buff_ptr = static_cast<int8_t*>(
      Agora_memory::PaddedAlignedAlloc(Agora_memory::Alignment_t::kAlign64,
                                       dl_bcast_mod_order_bits_ * ctrl_sc_num));
  Demodulate(reinterpret_cast<float*>(&eq_buff[0]), demod_buff_ptr,
             2 * ctrl_sc_num, dl_bcast_mod_order_bits_, false);

  const int num_bcast_bytes = BitsToBytes(dl_bcast_ldpc_config_.NumCbLen());
  std::vector<uint8_t> decode_buff(num_bcast_bytes, 0u);

  DataGenerator::GetDecodedData(demod_buff_ptr, &decode_buff.at(0),
                                dl_bcast_ldpc_config_, num_bcast_bytes,
                                cfg_->ScrambleEnabled());
  FreeBuffer1d(&bcast_fft_buff);
  FreeBuffer1d(&eq_buff);
  FreeBuffer1d(&demod_buff_ptr);
  const double duration =
      GetTime::CyclesToUs(GetTime::WorkerRdtsc() - start_tsc, cfg_->FreqGhz());
  if (kDebugPrintInTask) {
    std::printf("DecodeBroadcast completed in %2.2f us\n", duration);
  }
  return (reinterpret_cast<size_t*>(decode_buff.data()))[0];
}

void Mcs::GenBroadcastSlots(std::vector<std::complex<int16_t>*>& bcast_iq_samps,
                            std::vector<size_t> ctrl_msg) {
  //RtAssert(-1 == 1, "In Gen broadcast slots");

  std::cout << "In Gen broadcast slots" << std::endl << std::flush;
  ///\todo enable a vector of bytes to TX'ed in each symbol
  assert(bcast_iq_samps.size() == this->frame_.NumDlControlSyms());
  const size_t start_tsc = GetTime::WorkerRdtsc();

  int num_bcast_bytes = BitsToBytes(dl_bcast_ldpc_config_.NumCbLen());
  std::vector<int8_t> bcast_bits_buffer(num_bcast_bytes, 0);

  Table<complex_float> dl_bcast_mod_table;
  InitModulationTable(dl_bcast_mod_table, dl_bcast_mod_order_bits_);

  for (size_t i = 0; i < this->frame_.NumDlControlSyms(); i++) {
    std::memcpy(bcast_bits_buffer.data(), ctrl_msg.data(), sizeof(size_t));

    const auto coded_bits_ptr = DataGenerator::GenCodeblock(
        dl_bcast_ldpc_config_, &bcast_bits_buffer.at(0), num_bcast_bytes,
        cfg_->ScrambleEnabled());

    auto modulated_vector = DataGenerator::GetModulation(
        &coded_bits_ptr[0], dl_bcast_mod_table,
        dl_bcast_ldpc_config_.NumCbCodewLen(), cfg_->OfdmDataNum(),
        dl_bcast_mod_order_bits_);
    auto mapped_symbol = DataGenerator::MapOFDMSymbol(
        cfg_, modulated_vector, pilots_, SymbolType::kControl);
    auto ofdm_symbol = DataGenerator::BinForIfft(cfg_, mapped_symbol, true);
    CommsLib::IFFT(&ofdm_symbol[0], ofdm_ca_num_, false);
    // additional 2^2 (6dB) power backoff
    float dl_bcast_scale =
        2 * CommsLib::FindMaxAbs(&ofdm_symbol[0], ofdm_symbol.size());
    CommsLib::Ifft2tx(&ofdm_symbol[0], bcast_iq_samps[i], this->ofdm_ca_num_,
                      cfg_->OfdmTxZeroPrefix(), cfg_->CpLen(), dl_bcast_scale);
  }
  dl_bcast_mod_table.Free();
  const double duration =
      GetTime::CyclesToUs(GetTime::WorkerRdtsc() - start_tsc, cfg_->FreqGhz());
  if (kDebugPrintInTask) {
    std::printf("GenBroadcast completed in %2.2f us\n", duration);
  }
}

void Mcs::UpdateCtrlMCS() {
  if (this->frame_.NumDlControlSyms() > 0) {
    const size_t dl_bcast_mcs_index = kControlMCS;
    const size_t bcast_base_graph =
        1;  // TODO: For MCS < 5, base_graph 1 doesn't work
    dl_bcast_mod_order_bits_ = GetModOrderBits(dl_bcast_mcs_index);
    const size_t dl_bcast_code_rate = GetCodeRate(dl_bcast_mcs_index);
    std::string dl_bcast_modulation = MapModToStr(dl_bcast_mod_order_bits_);
    const int16_t max_decoder_iter = 5;
    size_t bcast_zc =
        SelectZc(bcast_base_graph, dl_bcast_code_rate, dl_bcast_mod_order_bits_,
                 cfg_->GetOFDMCtrlNum(), kCbPerSymbol, "downlink broadcast");

    // Always positive since dl_code_rate is smaller than 1
    size_t bcast_num_rows =
        static_cast<size_t>(std::round(
            1024.0 * LdpcNumInputCols(bcast_base_graph) / dl_bcast_code_rate)) -
        (LdpcNumInputCols(bcast_base_graph) - 2);

    uint32_t bcast_num_cb_len = LdpcNumInputBits(bcast_base_graph, bcast_zc);
    uint32_t bcast_num_cb_codew_len =
        LdpcNumEncodedBits(bcast_base_graph, bcast_zc, bcast_num_rows);
    dl_bcast_ldpc_config_ =
        LDPCconfig(bcast_base_graph, bcast_zc, max_decoder_iter, true,
                   bcast_num_cb_len, bcast_num_cb_codew_len, bcast_num_rows, 0);

    dl_bcast_ldpc_config_.NumBlocksInSymbol(
        (cfg_->GetOFDMCtrlNum() * dl_bcast_mod_order_bits_) /
        dl_bcast_ldpc_config_.NumCbCodewLen());
    RtAssert(dl_bcast_ldpc_config_.NumBlocksInSymbol() > 0,
             "Downlink Broadcast LDPC expansion factor is too large for number "
             "of OFDM data "
             "subcarriers.");
    AGORA_LOG_INFO(
        "Downlink Broadcast MCS Info: LDPC: Zc: %d, %zu code blocks per "
        "symbol, "
        "%d "
        "information "
        "bits per encoding, %d bits per encoded code word, decoder "
        "iterations: %d, code rate %.3f (nRows = %zu), modulation %s\n",
        dl_bcast_ldpc_config_.ExpansionFactor(),
        dl_bcast_ldpc_config_.NumBlocksInSymbol(),
        dl_bcast_ldpc_config_.NumCbLen(), dl_bcast_ldpc_config_.NumCbCodewLen(),
        dl_bcast_ldpc_config_.MaxDecoderIter(),
        1.f * LdpcNumInputCols(dl_bcast_ldpc_config_.BaseGraph()) /
            (LdpcNumInputCols(dl_bcast_ldpc_config_.BaseGraph()) - 2 +
             dl_bcast_ldpc_config_.NumRows()),
        dl_bcast_ldpc_config_.NumRows(), dl_bcast_modulation.c_str());
  }
}

void Mcs::DumpMcsInfo() {
  AGORA_LOG_INFO(
      "Uplink MCS Info: LDPC: Zc: %d, %zu code blocks per symbol, %d "
      "information "
      "bits per encoding, %d bits per encoded code word, decoder "
      "iterations: %d, code rate %.3f (nRows = %zu), modulation %s\n",
      ul_ldpc_config_.ExpansionFactor(), ul_ldpc_config_.NumBlocksInSymbol(),
      ul_ldpc_config_.NumCbLen(), ul_ldpc_config_.NumCbCodewLen(),
      ul_ldpc_config_.MaxDecoderIter(),
      1.f * LdpcNumInputCols(ul_ldpc_config_.BaseGraph()) /
          (LdpcNumInputCols(ul_ldpc_config_.BaseGraph()) - 2 +
           ul_ldpc_config_.NumRows()),
      ul_ldpc_config_.NumRows(), ul_modulation_.c_str());
  AGORA_LOG_INFO(
      "Downlink MCS Info: LDPC: Zc: %d, %zu code blocks per symbol, %d "
      "information "
      "bits per encoding, %d bits per encoded code word, decoder "
      "iterations: %d, code rate %.3f (nRows = %zu), modulation %s\n",
      dl_ldpc_config_.ExpansionFactor(), dl_ldpc_config_.NumBlocksInSymbol(),
      dl_ldpc_config_.NumCbLen(), dl_ldpc_config_.NumCbCodewLen(),
      dl_ldpc_config_.MaxDecoderIter(),
      1.f * LdpcNumInputCols(dl_ldpc_config_.BaseGraph()) /
          (LdpcNumInputCols(dl_ldpc_config_.BaseGraph()) - 2 +
           dl_ldpc_config_.NumRows()),
      dl_ldpc_config_.NumRows(), dl_modulation_.c_str());
}

LDPCconfig Mcs::UlLdpcConfig() { return ul_ldpc_config_; }

LDPCconfig Mcs::DlLdpcConfig() { return dl_ldpc_config_; }
