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
#include "ldpc_generator.h"
#include "logger.h"
#include "modulation.h"
#include "scrambler.h"
#include "simd_types.h"
#include "symbols.h"
#include "utils_ldpc.h"

//Need the second argument for GetOFDMDataNum
Mcs::Mcs(const OfdmConfig ofdm_data, McsParams ul_mcs_params,
         McsParams dl_mcs_params, FrameStats frame)
    : ul_mcs_params_(ul_mcs_params),
      dl_mcs_params_(dl_mcs_params),
      frame_(frame) {
  CreateModulationTables();
  //Initialize UL MCS

  ldpc_generator_ = std::make_unique<LdpcGenerator>(
      ofdm_data, ul_mcs_params.base_graph, ul_mcs_params.early_term_,
      ul_mcs_params.max_decoder_iter, mcs_params.dl_base_graph,
      dl_mcs_params.early_term, dl_mcs_params.max_decoder_iter);
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

McsScheme Mcs::InitializeMcs(McsParams mcs_params) {
  McsScheme current_mcs;
  current_mcs.frame_number_ = 0;

  //If the string has been set.
  if (mcs_params.modulation_ != NULL) {
    current_mcs.mod_order_bits_ = kModulStringMap.at(mcs_params.modulation_);

    size_t code_rate_int =
        static_cast<size_t>(std::round(mcs_params.code_rate_usr_ * 1024.0));

    current_mcs.mcs_index_ =
        CommsLib::GetMcsIndex(current_mcs.mod_order_bits_, code_rate_int);
    current_mcs.code_rate_ = GetCodeRate(current_mcs.mcs_index_);
    if (current_mcs.code_rate_ / 1024.0 != mcs_params.code_rate_usr_) {
      AGORA_LOG_WARN(
          "Rounded the user-defined uplink code rate to the closest standard "
          "rate %zu/1024.\n",
          current_mcs.code_rate_);
    }
  } else {
    current_mcs.mcs_index_ = mcs_params.mcs_index;
    current_mcs.mod_order_bits_ = GetModOrderBits(current_mcs.mcs_index_);
    current_mcs.code_rate_ = GetCodeRate(current_mcs.mcs_index_);
  }
  return current_mcs;
}

//This function is so I can check that the next mcs_frame isn't null before
// updating the mcs so I can know if it was actually set or not.
McsScheme Mcs::InitializeNextMcs() {
  McsScheme next_mcs;
  next_mcs.frame_number_ = NULL;
  next_mcs.code_rate_ = NULL;
  next_mcs.mcs_index_ = NULL;
  next_mcs.mod_order_bits_ = NULL;
}

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

McsScheme Mcs::UpdateCurrentUlMcs(McsScheme current_ul_mcs,
                                  McsScheme next_ul_mcs,
                                  size_t current_frame_number) {
  if (current_frame_number >= next_ul_mcs.frame_number_) {
    current_ul_mcs.frame_number_ = current_frame_number;
    current_ul_mcs.mcs_index_ = next_ul_mcs.mcs_index_;
  }
  return current_ul_mcs;
}

McsScheme Mcs::UpdateCurrentDlMcs(McsScheme current_dl_mcs,
                                  McsScheme next_dl_mcs,
                                  size_t current_frame_number) {
  if (current_frame_number >= next_dl_mcs.frame_number_) {
    current_dl_mcs.frame_number_ = current_frame_number;
    current_dl_mcs.mcs_index_ = next_dl_mcs.mcs_index_;
  }
  return current_dl_mcs;
}

void Mcs::SetNextUlMcs(size_t frame_number, size_t mod_order_bits) {
  next_ul_mcs_.frame_number_ = frame_number;
  next_ul_mcs_.mcs_index_ = mod_order_bits;
}

void Mcs::SetNextDlMcs(size_t frame_number, size_t mod_order_bits) {
  next_dl_mcs_.frame_number_ = frame_number;
  next_dl_mcs_.mcs_index_ = mod_order_bits;
}

//Needs to be updated to user new ldpc struct.
// void Mcs::DumpMcsInfo() {
//   AGORA_LOG_INFO(
//       "Uplink MCS Info: LDPC: Zc: %d, %zu code blocks per symbol, %d "
//       "information "
//       "bits per encoding, %d bits per encoded code word, decoder "
//       "iterations: %d, code rate %.3f (nRows = %zu), modulation %s\n",
//       ul_ldpc_config_.ExpansionFactor(), ul_ldpc_config_.NumBlocksInSymbol(),
//       ul_ldpc_config_.NumCbLen(), ul_ldpc_config_.NumCbCodewLen(),
//       ul_ldpc_config_.MaxDecoderIter(),
//       1.f * LdpcNumInputCols(ul_ldpc_config_.BaseGraph()) /
//           (LdpcNumInputCols(ul_ldpc_config_.BaseGraph()) - 2 +
//            ul_ldpc_config_.NumRows()),
//       ul_ldpc_config_.NumRows(), ul_modulation_.c_str());
//   AGORA_LOG_INFO(
//       "Downlink MCS Info: LDPC: Zc: %d, %zu code blocks per symbol, %d "
//       "information "
//       "bits per encoding, %d bits per encoded code word, decoder "
//       "iterations: %d, code rate %.3f (nRows = %zu), modulation %s\n",
//       dl_ldpc_config_.ExpansionFactor(), dl_ldpc_config_.NumBlocksInSymbol(),
//       dl_ldpc_config_.NumCbLen(), dl_ldpc_config_.NumCbCodewLen(),
//       dl_ldpc_config_.MaxDecoderIter(),
//       1.f * LdpcNumInputCols(dl_ldpc_config_.BaseGraph()) /
//           (LdpcNumInputCols(dl_ldpc_config_.BaseGraph()) - 2 +
//            dl_ldpc_config_.NumRows()),
//       dl_ldpc_config_.NumRows(), dl_modulation_.c_str());
// }
