/**
 * @file mcs.cc
 * @brief Class implementation for mcs handling
 */

#include "mcs.h"

#include <cstddef>

#include "comms-constants.inc"
#include "comms-lib.h"
#include "data_generator.h"
#include "datatype_conversion.h"
#include "ldpc_updater.h"
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

// , size_t ofdm_ctrl_num
// size_t ofdm_data_num, size_t get_ofdm_data_num,

const size_t ofdm_data_num_ul_;
const size_t ofdm_data_num_dl_;
const size_t ofdm_ctrl_data_;

//Need the second argument for GetOFDMDataNum
Mcs::Mcs(const OfdmConfig ofdm_data, UlMcsParams ul_mcs_params,
         DlMcsParams dl_mcs_params, FrameStats frame)
    : ofdm_data_num_ul_(ofdm_data.ofdm_data_num_ul_),
      ofdm_data_num_dl_(ofdm_data.ofdm_data_num_dl_),
      ofdm_ctrl_data_(ofdm_data.ofdm_ctrl_data_),
      ul_mcs_params_(ul_mcs_params),
      dl_mcs_params_(dl_mcs_params),
      ul_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      dl_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      dl_bcast_ldpc_config_(0, 0, 0, false, 0, 0, 0, 0),
      frame_(frame) {
  ldpc_updater = std::make_unique<LdpcUpdater> CreateModulationTables();
  //Initialize UL MCS
  InitializeUlMcs(ul_mcs_params_);

  //Initialize DL MCS
  InitializeDlMcs(dl_mcs_params_);
}

Mcs::~Mcs() = default;

void Mcs::CheckUlMcs(float snr, size_t frame_id) {
  float spectral_efficiency = SpectralEffeciency(snr);
  std::cout << "checking the dl mcs" << std::endl << std::flush;
  if (abs(kMcsIndexToSpectralEffeciency.at(current_dl_mcs_.mcs_index_) -
          spectral_efficiency) > 0.1) {
    //Find the closest spectral effeciency to the measured average special
    //effeciency and update the dl mcs to the corresponding dl mcs.

    float min_spectral_effeciency_delta = MAXFLOAT;
    size_t optimal_mcs_index = 0;

    for (auto map_iter = kMcsIndexToSpectralEffeciency.begin();
         map_iter != kMcsIndexToSpectralEffeciency.end(); ++map_iter) {
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
  current_ul_mcs_.dir = Direction::kUplink;
  current_ul_mcs_.frame_number_ = 0;
  if (ul_mcs.find("mcs_index") == ul_mcs.end()) {
    ul_modulation_ = ul_mcs.value("modulation", "16QAM");
    current_ul_mcs_.mod_order_bits_ = kModulStringMap.at(ul_modulation_);

    double ul_code_rate_usr = ul_mcs.value("code_rate", 0.333);
    size_t code_rate_int =
        static_cast<size_t>(std::round(ul_code_rate_usr * 1024.0));

    current_ul_mcs_.mcs_index_ =
        CommsLib::GetMcsIndex(current_ul_mcs_.mod_order_bits_, code_rate_int);
    current_ul_mcs_.code_rate_ = GetCodeRate(current_ul_mcs_.mcs_index_);
    if (current_ul_mcs_.code_rate_ / 1024.0 != ul_code_rate_usr) {
      AGORA_LOG_WARN(
          "Rounded the user-defined uplink code rate to the closest standard "
          "rate %zu/1024.\n",
          current_ul_mcs_.code_rate_);
    }
  } else {
    current_ul_mcs_.mcs_index_ =
        ul_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_ul_mcs_.mod_order_bits_ =
        GetModOrderBits(current_ul_mcs_.mcs_index_);
    current_ul_mcs_.code_rate_ = GetCodeRate(current_ul_mcs_.mcs_index_);
  }
}

void Mcs::InitializeDlMcs(const nlohmann::json dl_mcs) {
  current_dl_mcs_.dir = Direction::kDownlink;
  current_dl_mcs_.frame_number_ = 0;
  if (dl_mcs.find("mcs_index") == dl_mcs.end()) {
    dl_modulation_ = dl_mcs.value("modulation", "16QAM");
    current_dl_mcs_.mod_order_bits_ = kModulStringMap.at(dl_modulation_);

    double dl_code_rate_usr = dl_mcs.value("code_rate", 0.333);
    size_t code_rate_int =
        static_cast<size_t>(std::round(dl_code_rate_usr * 1024.0));
    current_dl_mcs_.mcs_index_ =
        CommsLib::GetMcsIndex(current_dl_mcs_.mod_order_bits_, code_rate_int);
    current_dl_mcs_.code_rate_ = GetCodeRate(current_dl_mcs_.mcs_index_);
    if (current_dl_mcs_.code_rate_ / 1024.0 != dl_code_rate_usr) {
      AGORA_LOG_WARN(
          "Rounded the user-defined downlink code rate to the closest standard "
          "rate %zu/1024.\n",
          current_dl_mcs_.code_rate_);
    }
  } else {
    current_dl_mcs_.mcs_index_ =
        dl_mcs.value("mcs_index", 10);  // 16QAM, 340/1024
    current_dl_mcs_.mod_order_bits_ =
        GetModOrderBits(current_dl_mcs_.mcs_index_);
    current_dl_mcs_.code_rate_ = GetCodeRate(current_dl_mcs_.mcs_index_);
  }
}

// inline size_t SelectZc(size_t base_graph, size_t code_rate,
//                        size_t mod_order_bits, size_t num_sc, size_t cb_per_sym,
//                        const std::string& dir) {
//   size_t n_zc = sizeof(kZc) / sizeof(size_t);
//   std::vector<size_t> zc_vec(kZc, kZc + n_zc);
//   std::sort(zc_vec.begin(), zc_vec.end());
//   // According to cyclic_shift.cc cyclic shifter for zc
//   // larger than 256 has not been implemented, so we skip them here.
//   size_t max_zc_index =
//       (std::find(zc_vec.begin(), zc_vec.end(), kMaxSupportedZc) -
//        zc_vec.begin());
//   size_t max_uncoded_bits =
//       static_cast<size_t>(num_sc * code_rate * mod_order_bits / 1024.0);
//   size_t zc = SIZE_MAX;
//   size_t i = 0;
//   for (; i < max_zc_index; i++) {
//     if ((zc_vec.at(i) * LdpcNumInputCols(base_graph) * cb_per_sym <
//          max_uncoded_bits) &&
//         (zc_vec.at(i + 1) * LdpcNumInputCols(base_graph) * cb_per_sym >
//          max_uncoded_bits)) {
//       zc = zc_vec.at(i);
//       break;
//     }
//   }
//   if (zc == SIZE_MAX) {
//     AGORA_LOG_WARN(
//         "Exceeded possible range of LDPC lifting Zc for " + dir +
//             "! Setting lifting size to max possible value(%zu).\nThis may lead "
//             "to too many unused subcarriers. For better use of the PHY "
//             "resources, you may reduce your coding or modulation rate.\n",
//         kMaxSupportedZc);
//     zc = kMaxSupportedZc;
//   }
//   return zc;
// }

void Mcs::CreateModulationTables() {
  std::cout << "CREATING MODULATION TABLES" << std::endl << std::flush;

  for (size_t i = 0; i < kNumTables; i++) {
    std::cout << "creating ul table: " << std::endl << std::flush;
    InitModulationTable(modulation_tables_.ul_tables_[i], (i + 1) * 2);
    std::cout << modulation_tables_.ul_tables_[i][0][0].re << " "
              << modulation_tables_.ul_tables_[i][0][0].im << std::endl
              << std::flush;

    std::cout << "creating dl table" << std::endl << std::flush;
    InitModulationTable(modulation_tables_.dl_tables_[i], (i + 1) * 2);
    std::cout << "DL table: at 0, 0: "
              << modulation_tables_.dl_tables_[i][0][0].re << " "
              << modulation_tables_.dl_tables_[i][0][0].im << std::endl
              << std::flush;
  }
}

void Mcs::UpdateMcs(size_t current_frame_number) {
  std::cout << "Updating MCS." << std::endl << std::flush;
  UpdateUlMcs(current_frame_number);
  UpdateDlMcs(current_frame_number);
}

void Mcs::UpdateUlMcs(size_t current_frame_number) {
  if (current_frame_number >= next_ul_mcs_.frame_number_) {
    current_ul_mcs_.frame_number_ = current_frame_number;
    current_ul_mcs_.mcs_index_ = next_ul_mcs_.mcs_index_;
  }
  UpdateUlLdpcConfig();
}

void Mcs::UpdateDlMcs(size_t current_frame_number) {
  if (current_frame_number >= next_dl_mcs_.frame_number_) {
    current_dl_mcs_.frame_number_ = current_frame_number;
    current_dl_mcs_.mcs_index_ = next_dl_mcs_.mcs_index_;
  }
  UpdateDlLdpcConfig();
}

void Mcs::SetNextUlMcs(size_t frame_number, size_t mod_order_bits) {
  next_ul_mcs_.frame_number_ = frame_number;
  next_ul_mcs_.mcs_index_ = mod_order_bits;
}

void Mcs::SetNextDlMcs(size_t frame_number, size_t mod_order_bits) {
  next_dl_mcs_.frame_number_ = frame_number;
  next_dl_mcs_.mcs_index_ = mod_order_bits;
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
