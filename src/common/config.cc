// Copyright (c) 2018-2022, Rice University
// RENEW OPEN SOURCE LICENSE: http://renew-wireless.org/license

/**
 * @file config.cc
 * @brief Implementation file for the configuration class which importants
 * json configuration values into class variables
 */

#include "config.h"

#include <ctime>
#include <filesystem>
#include <utility>

#include "comms-constants.inc"
#include "comms-lib.h"
#include "data_generator.h"
#include "datatype_conversion.h"
#include "fivegconfig.h"
#include "gettime.h"
#include "logger.h"
#include "mcs.h"
#include "message.h"
#include "modulation.h"
#include "phy_ldpc_decoder_5gnr.h"
#include "scrambler.h"
#include "simd_types.h"
#include "utils_ldpc.h"

using json = nlohmann::json;

static constexpr size_t kMacAlignmentBytes = 64u;
static constexpr bool kDebugPrintConfiguration = false;
static constexpr size_t kShortIdLen = 3;
static constexpr size_t kVarNodesSize = 1024 * 1024 * sizeof(int16_t);
static constexpr size_t kControlMCS = 5;  // QPSK, 379/1024

//Eventually Will have to fully migrate this to MCS.
static const std::string kLogFilepath =
    TOSTRING(PROJECT_DIRECTORY) "/files/log/";
static const std::string kExperimentFilepath =
    TOSTRING(PROJECT_DIRECTORY) "/files/experiment/";
static const std::string kUlDataFilePrefix =
    kExperimentFilepath + "LDPC_orig_ul_data_";
static const std::string kDlDataFilePrefix =
    kExperimentFilepath + "LDPC_orig_dl_data_";
static const std::string kUlDataFreqPrefix = kExperimentFilepath + "ul_data_f_";

/// Print the I/Q samples in the pilots
static constexpr bool kDebugPrintPilot = false;

Config::Config(std::string jsonfilename)
    : freq_ghz_(GetTime::MeasureRdtscFreq()),
      frame_(""),
      config_filename_(std::move(jsonfilename)) {
  auto time = std::time(nullptr);
  auto local_time = *std::localtime(&time);
  timestamp_ = std::to_string(1900 + local_time.tm_year) + "-" +
               std::to_string(1 + local_time.tm_mon) + "-" +
               std::to_string(local_time.tm_mday) + "-" +
               std::to_string(local_time.tm_hour) + "-" +
               std::to_string(local_time.tm_min) + "-" +
               std::to_string(local_time.tm_sec);

  std::string conf;
  Utils::LoadTddConfig(config_filename_, conf);
  // Allow json comments
  const auto tdd_conf = json::parse(conf, nullptr, true, true);

  // Initialize the compute configuration
  // Default exclude 1 core with id = 0
  std::vector<size_t> excluded(1, 0);
  if (tdd_conf.contains("exclude_cores")) {
    auto exclude_cores = tdd_conf.at("exclude_cores");
    excluded.resize(exclude_cores.size());
    for (size_t i = 0; i < exclude_cores.size(); i++) {
      excluded.at(i) = exclude_cores.at(i);
    }
  }
  SetCpuLayoutOnNumaNodes(true, excluded);

  num_cells_ = tdd_conf.value("cells", 1);
  num_radios_ = 0;
  ue_num_ = 0;

  std::string serials_str;
  std::string serial_file = tdd_conf.value("serial_file", "");
  if (serial_file.empty() == false) {
    Utils::LoadTddConfig(serial_file, serials_str);
  }
  if (serials_str.empty() == false) {
    const auto j_serials = json::parse(serials_str, nullptr, true, true);

    std::stringstream ss;
    json j_bs_serials;
    ss << j_serials.value("BaseStations", j_bs_serials);
    j_bs_serials = json::parse(ss);
    ss.str(std::string());
    ss.clear();

    RtAssert(j_bs_serials.size() == num_cells_, "Incorrect cells number!");
    external_ref_node_.resize(num_cells_, false);
    for (size_t i = 0; i < num_cells_; i++) {
      json serials_conf;
      std::string cell_str = "BS" + std::to_string(i);
      ss << j_bs_serials.value(cell_str, serials_conf);
      serials_conf = json::parse(ss);
      ss.str(std::string());
      ss.clear();

      auto hub_serial = serials_conf.value("hub", "");
      hub_id_.push_back(hub_serial);
      auto sdr_serials = serials_conf.value("sdr", json::array());
      RtAssert(!sdr_serials.empty(), "BS has zero sdrs!");
      radio_id_.insert(radio_id_.end(), sdr_serials.begin(), sdr_serials.end());
      num_radios_ += sdr_serials.size();
      cell_id_.resize(num_radios_, i);

      auto refnode_serial = serials_conf.value("reference", "");
      if (refnode_serial.empty()) {
        AGORA_LOG_INFO(
            "No reference node ID found in topology file! Taking the last node "
            "%s as reference node!\n",
            radio_id_.back().c_str());
        refnode_serial = radio_id_.back();
        ref_radio_.push_back(radio_id_.size() - 1);
      } else {
        auto serial_iterator =
            std::find(sdr_serials.begin(), sdr_serials.end(), refnode_serial);
        if (serial_iterator == sdr_serials.end()) {
          radio_id_.push_back(refnode_serial);
          ref_radio_.push_back(radio_id_.size() - 1);
          num_radios_++;
          cell_id_.resize(num_radios_, i);
          external_ref_node_.at(i) = true;
        } else {
          size_t index = radio_id_.size() - sdr_serials.size() +
                         serial_iterator - sdr_serials.begin();
          ref_radio_.push_back(index);
        }
      }
    }

    json j_ue_serials;
    ss << j_serials.value("Clients", j_ue_serials);
    j_ue_serials = json::parse(ss);
    ss.str(std::string());
    ss.clear();

    auto ue_serials = j_ue_serials.value("sdr", json::array());
    ue_radio_id_.assign(ue_serials.begin(), ue_serials.end());
  } else if (kUseArgos == true) {
    throw std::runtime_error(
        "Hardware is enabled but the serials files was not accessable");
  }

  if (radio_id_.empty()) {
    num_radios_ = tdd_conf.value("bs_radio_num", 8);
    external_ref_node_.resize(num_cells_, false);
    cell_id_.resize(num_radios_, 0);

    //Add in serial numbers
    for (size_t radio = 0; radio < num_radios_; radio++) {
      AGORA_LOG_TRACE("Adding BS_SIM_RADIO_%d\n", radio);
      radio_id_.emplace_back("BS_SIM_RADIO_" + std::to_string(radio));
    }
  }

  if (ue_radio_id_.empty()) {
    ue_num_ = tdd_conf.value("ue_radio_num", 8);
    for (size_t ue_radio = 0; ue_radio < ue_num_; ue_radio++) {
      std::stringstream ss;
      ss << std::setw(kShortIdLen) << std::setfill('0') << ue_radio;
      const std::string ue_name = "UE_SIM_RADIO_" + ss.str();
      AGORA_LOG_TRACE("Adding %s\n", ue_name.c_str());
      ue_radio_id_.push_back(ue_name);
    }
  }
  ue_num_ = ue_radio_id_.size();
  for (size_t i = 0; i < ue_num_; i++) {
    ue_radio_name_.push_back(
        "UE" + (ue_radio_id_.at(i).length() > kShortIdLen
                    ? ue_radio_id_.at(i).substr(ue_radio_id_.at(i).length() -
                                                kShortIdLen)
                    : ue_radio_id_.at(i)));
  }

  channel_ = tdd_conf.value("channel", "A");
  ue_channel_ = tdd_conf.value("ue_channel", channel_);
  num_channels_ = std::min(channel_.size(), kMaxChannels);
  num_ue_channels_ = std::min(ue_channel_.size(), kMaxChannels);
  bs_ant_num_ = num_channels_ * num_radios_;
  ue_ant_num_ = ue_num_ * num_ue_channels_;

  bf_ant_num_ = bs_ant_num_;
  for (size_t i = 0; i < num_cells_; i++) {
    if (external_ref_node_.at(i) == true) {
      bf_ant_num_ = bs_ant_num_ - num_channels_;
    }
  }

  if (ref_radio_.empty() == false) {
    for (size_t i = 0; i < num_cells_; i++) {
      ref_ant_.push_back(ref_radio_.at(i) * num_channels_);
    }
  }

  if ((kUseArgos == true) || (kUseUHD == true) || (kUsePureUHD == true)) {
    RtAssert(num_radios_ != 0, "Error: No radios exist in Argos mode");
  }

  /* radio configurations */
  freq_ = tdd_conf.value("frequency", 3.6e9);
  single_gain_ = tdd_conf.value("single_gain", true);
  tx_gain_a_ = tdd_conf.value("tx_gain_a", 20);
  rx_gain_a_ = tdd_conf.value("rx_gain_a", 20);
  tx_gain_b_ = tdd_conf.value("tx_gain_b", 20);
  rx_gain_b_ = tdd_conf.value("rx_gain_b", 20);
  calib_tx_gain_a_ = tdd_conf.value("calib_tx_gain_a", tx_gain_a_);
  calib_tx_gain_b_ = tdd_conf.value("calib_tx_gain_b", tx_gain_b_);

  auto gain_tx_json_a = tdd_conf.value("ue_tx_gain_a", json::array());
  if (gain_tx_json_a.empty()) {
    client_tx_gain_a_.resize(ue_num_, 20);
  } else {
    RtAssert(gain_tx_json_a.size() == ue_num_,
             "ue_tx_gain_a size must be same as the number of clients!");
    client_tx_gain_a_.assign(gain_tx_json_a.begin(), gain_tx_json_a.end());
  }
  auto gain_tx_json_b = tdd_conf.value("ue_tx_gain_b", json::array());
  if (gain_tx_json_b.empty()) {
    client_tx_gain_b_.resize(ue_num_, 0);
  } else {
    RtAssert(gain_tx_json_b.size() == ue_num_,
             "ue_tx_gain_b size must be same as the number of clients!");
    client_tx_gain_b_.assign(gain_tx_json_b.begin(), gain_tx_json_b.end());
  }
  auto gain_rx_json_a = tdd_conf.value("ue_rx_gain_a", json::array());
  if (gain_rx_json_a.empty()) {
    client_rx_gain_a_.resize(ue_num_, 20);
  } else {
    RtAssert(gain_rx_json_a.size() == ue_num_,
             "ue_rx_gain_a size must be same as the number of clients!");
    client_rx_gain_a_.assign(gain_rx_json_a.begin(), gain_rx_json_a.end());
  }
  auto gain_rx_json_b = tdd_conf.value("ue_rx_gain_b", json::array());
  if (gain_rx_json_b.empty()) {
    client_rx_gain_b_.resize(ue_num_, 0);
  } else {
    RtAssert(gain_rx_json_b.size() == ue_num_,
             "ue_rx_gain_b size must be same as the number of clients!");
    client_rx_gain_b_.assign(gain_rx_json_b.begin(), gain_rx_json_b.end());
  }

  rate_ = tdd_conf.value("sample_rate", 5e6);
  nco_ = tdd_conf.value("nco_frequency", 0.75 * rate_);
  bw_filter_ = rate_ + 2 * nco_;
  radio_rf_freq_ = freq_ - nco_;
  beacon_ant_ = tdd_conf.value("beacon_antenna", 0);
  beamsweep_ = tdd_conf.value("beamsweep", false);
  sample_cal_en_ = tdd_conf.value("calibrate_digital", false);
  imbalance_cal_en_ = tdd_conf.value("calibrate_analog", false);
  init_calib_repeat_ = tdd_conf.value("init_calib_repeat", 0);
  smooth_calib_ = tdd_conf.value("smooth_calib", false);
  beamforming_str_ = tdd_conf.value("beamforming", "ZF");
  beamforming_algo_ = kBeamformingStr.at(beamforming_str_);
  num_spatial_streams_ = tdd_conf.value("spatial_streams", ue_ant_num_);

  bs_server_addr_ = tdd_conf.value("bs_server_addr", "127.0.0.1");
  bs_rru_addr_ = tdd_conf.value("bs_rru_addr", "127.0.0.1");
  ue_server_addr_ = tdd_conf.value("ue_server_addr", "127.0.0.1");
  ue_rru_addr_ = tdd_conf.value("ue_rru_addr", "127.0.0.1");
  mac_remote_addr_ = tdd_conf.value("mac_remote_addr", "127.0.0.1");
  bs_server_port_ = tdd_conf.value("bs_server_port", 8000);
  bs_rru_port_ = tdd_conf.value("bs_rru_port", 9000);
  ue_rru_port_ = tdd_conf.value("ue_rru_port", 7000);
  ue_server_port_ = tdd_conf.value("ue_server_port", 6000);

  dpdk_num_ports_ = tdd_conf.value("dpdk_num_ports", 1);
  dpdk_port_offset_ = tdd_conf.value("dpdk_port_offset", 0);
  dpdk_mac_addrs_ = tdd_conf.value("dpdk_mac_addrs", "");

  ue_mac_tx_port_ = tdd_conf.value("ue_mac_tx_port", kMacUserRemotePort);
  ue_mac_rx_port_ = tdd_conf.value("ue_mac_rx_port", kMacUserLocalPort);
  bs_mac_tx_port_ = tdd_conf.value("bs_mac_tx_port", kMacBaseRemotePort);
  bs_mac_rx_port_ = tdd_conf.value("bs_mac_rx_port", kMacBaseLocalPort);

  log_listener_addr_ = tdd_conf.value("log_listener_addr", "");
  log_listener_port_ = tdd_conf.value("log_listener_port", 33300);

  log_sc_num_ = tdd_conf.value("log_sc_num", 4);
  log_timestamp_ = tdd_conf.value("log_timestamp", false);

  /* frame configurations */
  cp_len_ = tdd_conf.value("cp_size", 0);
  ofdm_ca_num_ = tdd_conf.value("fft_size", 2048);
  ofdm_data_num_ = tdd_conf.value("ofdm_data_num", 1200);
  ofdm_tx_zero_prefix_ = tdd_conf.value("ofdm_tx_zero_prefix", 0);
  ofdm_tx_zero_postfix_ = tdd_conf.value("ofdm_tx_zero_postfix", 0);
  ofdm_rx_zero_prefix_bs_ =
      tdd_conf.value("ofdm_rx_zero_prefix_bs", 0) + cp_len_;
  ofdm_rx_zero_prefix_client_ = tdd_conf.value("ofdm_rx_zero_prefix_client", 0);
  ofdm_rx_zero_prefix_cal_ul_ =
      tdd_conf.value("ofdm_rx_zero_prefix_cal_ul", 0) + cp_len_;
  ofdm_rx_zero_prefix_cal_dl_ =
      tdd_conf.value("ofdm_rx_zero_prefix_cal_dl", 0) + cp_len_;
  RtAssert(cp_len_ % 16 == 0,
           "cyclic prefix must be a multiple of subcarriers "
           "per cacheline.");
  RtAssert(ofdm_data_num_ % kSCsPerCacheline == 0,
           "ofdm_data_num must be a multiple of subcarriers per cacheline");
  RtAssert(ofdm_data_num_ % kTransposeBlockSize == 0,
           "Transpose block size must divide number of OFDM data subcarriers");
  ofdm_pilot_spacing_ = tdd_conf.value("ofdm_pilot_spacing", 16);
  ofdm_data_start_ = tdd_conf.value("ofdm_data_start",
                                    ((ofdm_ca_num_ - ofdm_data_num_) / 2) /
                                        kSCsPerCacheline * kSCsPerCacheline);
  RtAssert(ofdm_data_start_ % kSCsPerCacheline == 0,
           "ofdm_data_start must be a multiple of subcarriers per cacheline");
  ofdm_data_stop_ = ofdm_data_start_ + ofdm_data_num_;

  // Build subcarrier map for data ofdm symbols
  ul_symbol_map_.resize(ofdm_data_num_, SubcarrierType::kData);
  dl_symbol_map_.resize(ofdm_data_num_);
  control_symbol_map_.resize(ofdm_data_num_);
  // Maps subcarrier index to data index
  dl_symbol_data_id_.resize(ofdm_data_num_, 0);
  dl_symbol_ctrl_id_.resize(ofdm_data_num_, 0);
  size_t data_idx = 0;
  size_t ctrl_idx = 0;
  for (size_t i = 0; i < ofdm_data_num_; i++) {
    if (i % ofdm_pilot_spacing_ == 0) {  // TODO: make this index configurable
      dl_symbol_map_.at(i) = SubcarrierType::kDMRS;
      control_symbol_map_.at(i) = SubcarrierType::kDMRS;
    } else {
      dl_symbol_map_.at(i) = SubcarrierType::kData;
      dl_symbol_data_id_.at(i) = data_idx;
      data_idx++;
      if (i % ofdm_pilot_spacing_ == 1) {
        control_symbol_map_.at(i) = SubcarrierType::kPTRS;
      } else {
        control_symbol_map_.at(i) = SubcarrierType::kData;
        dl_symbol_ctrl_id_.at(i) = ctrl_idx;
        ctrl_idx++;
      }
    }
  }

  bigstation_mode_ = tdd_conf.value("bigstation_mode", false);
  freq_orthogonal_pilot_ = tdd_conf.value("freq_orthogonal_pilot", false);
  pilot_sc_group_size_ =
      tdd_conf.value("pilot_sc_group_size", kTransposeBlockSize);
  if (freq_orthogonal_pilot_) {
    RtAssert(pilot_sc_group_size_ == kTransposeBlockSize,
             "In this version, pilot_sc_group_size must be equal to Transpose "
             "Block Size " +
                 std::to_string(kTransposeBlockSize));
    RtAssert(ofdm_data_num_ % pilot_sc_group_size_ == 0,
             "ofdm_data_num must be evenly divided by pilot_sc_group_size " +
                 std::to_string(pilot_sc_group_size_));
    RtAssert(ue_ant_num_ <= pilot_sc_group_size_,
             "user antennas must be no more than pilot_sc_group_size " +
                 std::to_string(pilot_sc_group_size_));
  }
  correct_phase_shift_ = tdd_conf.value("correct_phase_shift", false);

  hw_framer_ = tdd_conf.value("hw_framer", true);
  if (kUseUHD || kUsePureUHD) {
    hw_framer_ = false;
  } else {
    RtAssert(hw_framer_ == true,
             "Base Station hardware framer (hw_framer) set to false is "
             "unsupported in this version of Agora");
  }
  ue_hw_framer_ = tdd_conf.value("ue_hw_framer", false);
  RtAssert(ue_hw_framer_ == false,
           "User equiptment hardware framer (ue_hw_framer) set to true is "
           "unsupported in this version of Agora");
  ue_resync_period_ = tdd_conf.value("ue_resync_period", 0);

  // If frames not specified explicitly, construct default based on frame_type /
  // symbol_num_perframe / pilot_num / ul_symbol_num_perframe /
  // dl_symbol_num_perframe / dl_data_symbol_start
  if (tdd_conf.find("frame_schedule") == tdd_conf.end()) {
    size_t ul_data_symbol_num_perframe = kDefaultULSymPerFrame;
    size_t ul_data_symbol_start = kDefaultULSymStart;
    size_t dl_data_symbol_num_perframe = kDefaultDLSymPerFrame;
    size_t dl_data_symbol_start = kDefaultDLSymStart;

    size_t symbol_num_perframe =
        tdd_conf.value("symbol_num_perframe", kDefaultSymbolNumPerFrame);
    size_t pilot_symbol_num_perframe = tdd_conf.value(
        "pilot_num",
        freq_orthogonal_pilot_ ? kDefaultFreqOrthPilotSymbolNum : ue_ant_num_);

    size_t beacon_symbol_position = tdd_conf.value("beacon_position", SIZE_MAX);

    ul_data_symbol_num_perframe =
        tdd_conf.value("ul_symbol_num_perframe", ul_data_symbol_num_perframe);

    if (ul_data_symbol_num_perframe == 0) {
      ul_data_symbol_start = 0;
    } else {
      // Start position of the first UL symbol
      ul_data_symbol_start =
          tdd_conf.value("ul_data_symbol_start", ul_data_symbol_start);
    }
    const size_t ul_data_symbol_stop =
        ul_data_symbol_start + ul_data_symbol_num_perframe;

    //Dl symbols
    dl_data_symbol_num_perframe =
        tdd_conf.value("dl_symbol_num_perframe", dl_data_symbol_num_perframe);

    if (dl_data_symbol_num_perframe == 0) {
      dl_data_symbol_start = 0;
    } else {
      // Start position of the first DL symbol
      dl_data_symbol_start =
          tdd_conf.value("dl_data_symbol_start", dl_data_symbol_start);
    }
    const size_t dl_data_symbol_stop =
        dl_data_symbol_start + dl_data_symbol_num_perframe;

    if ((ul_data_symbol_num_perframe + dl_data_symbol_num_perframe +
         pilot_symbol_num_perframe) > symbol_num_perframe) {
      AGORA_LOG_ERROR(
          "!!!!! Invalid configuration pilot + ul + dl exceeds total symbols "
          "!!!!!\n");
      AGORA_LOG_ERROR(
          "Uplink symbols: %zu, Downlink Symbols :%zu, Pilot Symbols: %zu, "
          "Total Symbols: %zu\n",
          ul_data_symbol_num_perframe, dl_data_symbol_num_perframe,
          pilot_symbol_num_perframe, symbol_num_perframe);
      throw std::runtime_error("Invalid Frame Configuration");
    } else if (((ul_data_symbol_num_perframe > 0) &&
                (dl_data_symbol_num_perframe > 0)) &&
               (((ul_data_symbol_start >= dl_data_symbol_start) &&
                 (ul_data_symbol_start < dl_data_symbol_stop)) ||
                ((ul_data_symbol_stop > dl_data_symbol_start) &&
                 (ul_data_symbol_stop <= dl_data_symbol_stop)))) {
      AGORA_LOG_ERROR(
          "!!!!! Invalid configuration ul and dl symbol overlap detected "
          "!!!!!\n");
      AGORA_LOG_ERROR(
          "Uplink - start: %zu - stop :%zu, Downlink - start: %zu - stop %zu\n",
          ul_data_symbol_start, ul_data_symbol_stop, dl_data_symbol_start,
          dl_data_symbol_stop);
      throw std::runtime_error("Invalid Frame Configuration");
    }

    char first_sym;
    char second_sym;
    size_t first_sym_start;
    size_t first_sym_count;
    size_t second_sym_start;
    size_t second_sym_count;
    if ((dl_data_symbol_num_perframe > 0) &&
        (dl_data_symbol_start <= ul_data_symbol_start)) {
      first_sym = 'D';
      first_sym_start = dl_data_symbol_start;
      first_sym_count = dl_data_symbol_num_perframe;
      second_sym = 'U';
      second_sym_start = ul_data_symbol_start;
      second_sym_count = ul_data_symbol_num_perframe;
    } else {
      first_sym = 'U';
      first_sym_start = ul_data_symbol_start;
      first_sym_count = ul_data_symbol_num_perframe;
      second_sym = 'D';
      second_sym_start = dl_data_symbol_start;
      second_sym_count = dl_data_symbol_num_perframe;
    }
    AGORA_LOG_SYMBOL(
        "Symbol %c, start %zu, count %zu. Symbol %c, start %zu, count %zu. "
        "Total Symbols: %zu\n",
        first_sym, first_sym_start, first_sym_start, second_sym,
        second_sym_start, second_sym_start, symbol_num_perframe);

    std::string sched = "";
    // Offset the pilots, if the beacon comes first
    if (beacon_symbol_position == 0) {
      sched = "G";
    }
    sched.append(pilot_symbol_num_perframe, 'P');
    // ( )PGGGG1111111111GGGG2222222222GGGG
    if (first_sym_start > 0) {
      const int guard_symbols = first_sym_start - sched.length();
      if (guard_symbols > 0) {
        sched.append(guard_symbols, 'G');
      }
      if (first_sym_count > 0) {
        sched.append(first_sym_count, first_sym);
      }
    }
    if (second_sym_start > 0) {
      const int guard_symbols = second_sym_start - sched.length();
      if (guard_symbols > 0) {
        sched.append(guard_symbols, 'G');
      }
      if (second_sym_count > 0) {
        sched.append(second_sym_count, second_sym);
      }
    }
    const int guard_symbols = symbol_num_perframe - sched.length();
    if (guard_symbols > 0) {
      sched.append(guard_symbols, 'G');
    }

    // Add the beacon
    if (beacon_symbol_position < sched.length()) {
      if (sched.at(beacon_symbol_position) != 'G') {
        AGORA_LOG_ERROR("Invalid beacon location %zu replacing %c\n",
                        beacon_symbol_position,
                        sched.at(beacon_symbol_position));
        throw std::runtime_error("Invalid Frame Configuration");
      }
      sched.replace(beacon_symbol_position, 1, "B");
    }
    frame_ = FrameStats(sched);
  } else {
    json jframes = tdd_conf.value("frame_schedule", json::array());

    // Only allow 1 unique frame type
    assert(jframes.size() == 1);
    std::string frame = jframes.at(0).get<std::string>();
    /*
    If an apostrophe delimiter is found in the frame string, execute logic to
    convert a subframe formated frame into the symbol formated frame that Agora
    is designed to handle.
    */
    if (frame.find(",") != std::string::npos) {
      std::vector<std::string> flex_formats =
          tdd_conf.value("flex_formats", json::array());
      FiveGConfig fivegconfig = FiveGConfig(tdd_conf);
      frame = fivegconfig.FiveGFormat();
      rate_ = fivegconfig.SamplingRate();
      ofdm_data_start_ = fivegconfig.OfdmDataStart();
    }
    frame_ = FrameStats(frame);
  }
  AGORA_LOG_INFO("Config: Frame schedule %s (%zu symbols)\n",
                 frame_.FrameIdentifier().c_str(), frame_.NumTotalSyms());

  if (frame_.IsRecCalEnabled()) {
    RtAssert(bf_ant_num_ >= frame_.NumDLCalSyms(),
             "Too many DL Cal symbols for the number of base station antennas");

    RtAssert(((bf_ant_num_ % frame_.NumDLCalSyms()) == 0),
             "Number of Downlink calibration symbols per frame must complete "
             "calibration on frame boundary!");
  }
  // Check for frame validity.
  // We should remove the restriction of the beacon symbol placement when tested
  // more thoroughly
  if (((frame_.NumBeaconSyms() > 1)) ||
      ((frame_.NumBeaconSyms() == 1) && (frame_.GetBeaconSymbolLast() > 1))) {
    AGORA_LOG_ERROR("Invalid beacon symbol placement\n");
    throw std::runtime_error("Invalid beacon symbol placement");
  }

  // client_dl_pilot_sym uses the first x 'D' symbols for downlink channel
  // estimation for each user.
  size_t client_dl_pilot_syms = tdd_conf.value("client_dl_pilot_syms", 0);
  RtAssert(client_dl_pilot_syms <= frame_.NumDLSyms(),
           "Number of ul_iq_t_DL pilot symbol exceeds number of DL symbols!");
  // client_ul_pilot_sym uses the first x 'U' symbols for downlink channel
  // estimation for each user.
  size_t client_ul_pilot_syms = tdd_conf.value("client_ul_pilot_syms", 0);
  RtAssert(client_ul_pilot_syms <= frame_.NumULSyms(),
           "Number of UL pilot symbol exceeds number of UL symbols!");

  frame_.SetClientPilotSyms(client_ul_pilot_syms, client_dl_pilot_syms);

  if ((freq_orthogonal_pilot_ == false) &&
      (ue_ant_num_ != frame_.NumPilotSyms())) {
    RtAssert(
        false,
        "Number of pilot symbols: " + std::to_string(frame_.NumPilotSyms()) +
            " does not match number of UEs: " + std::to_string(ue_ant_num_));
  }
  if ((freq_orthogonal_pilot_ == false) && (ue_radio_id_.empty() == true) &&
      (tdd_conf.find("ue_radio_num") == tdd_conf.end())) {
    ue_num_ = frame_.NumPilotSyms();
    ue_ant_num_ = ue_num_ * num_ue_channels_;
  }
  ue_ant_offset_ = tdd_conf.value("ue_ant_offset", 0);
  ue_ant_total_ = tdd_conf.value("ue_ant_total", ue_ant_num_);

  auto tx_advance = tdd_conf.value("tx_advance", json::array());
  if (tx_advance.empty()) {
    cl_tx_advance_.resize(ue_num_, 0);
  } else {
    RtAssert(tx_advance.size() == ue_num_,
             "tx_advance size must be same as the number of clients!");
    cl_tx_advance_.assign(tx_advance.begin(), tx_advance.end());
  }

  auto corr_scale = tdd_conf.value("corr_scale", json::array());
  if (corr_scale.empty()) {
    cl_corr_scale_.resize(ue_num_, 1.f);
  } else {
    RtAssert(corr_scale.size() == ue_num_,
             "corr_scale size must be same as the number of clients!");
    cl_corr_scale_.assign(corr_scale.begin(), corr_scale.end());
  }

  if (std::filesystem::is_directory(kExperimentFilepath) == false) {
    std::filesystem::create_directory(kExperimentFilepath);
  }

  // set trace file path
  const std::string ul_present_str = (frame_.NumULSyms() > 0 ? "uplink-" : "");
  const std::string dl_present_str =
      (frame_.NumDLSyms() > 0 ? "downlink-" : "");
  std::string filename =
      kLogFilepath + "trace-" + ul_present_str + dl_present_str + timestamp_ +
      "_" + std::to_string(num_cells_) + "_" + std::to_string(BsAntNum()) +
      "x" + std::to_string(UeAntTotal()) + ".hdf5";
  trace_file_ = tdd_conf.value("trace_file", filename);

  // Agora configurations
  frames_to_test_ = tdd_conf.value("max_frame", 9600);
  core_offset_ = tdd_conf.value("core_offset", 0);
  worker_thread_num_ = tdd_conf.value("worker_thread_num", 25);
  socket_thread_num_ = tdd_conf.value("socket_thread_num", 4);
  ue_core_offset_ = tdd_conf.value("ue_core_offset", 0);
  ue_worker_thread_num_ = tdd_conf.value("ue_worker_thread_num", 25);
  ue_socket_thread_num_ = tdd_conf.value("ue_socket_thread_num", 4);
  fft_thread_num_ = tdd_conf.value("fft_thread_num", 5);
  demul_thread_num_ = tdd_conf.value("demul_thread_num", 5);
  decode_thread_num_ = tdd_conf.value("decode_thread_num", 10);
  beam_thread_num_ = worker_thread_num_ - fft_thread_num_ - demul_thread_num_ -
                     decode_thread_num_;

  demul_block_size_ = tdd_conf.value("demul_block_size", 48);
  RtAssert(demul_block_size_ % kSCsPerCacheline == 0,
           "Demodulation block size must be a multiple of subcarriers per "
           "cacheline");
  RtAssert(
      demul_block_size_ % kTransposeBlockSize == 0,
      "Demodulation block size must be a multiple of transpose block size");
  demul_events_per_symbol_ = 1 + (ofdm_data_num_ - 1) / demul_block_size_;

  beam_block_size_ = tdd_conf.value("beam_block_size", 1);
  if (freq_orthogonal_pilot_) {
    if (beam_block_size_ == 1) {
      AGORA_LOG_INFO("Setting beam_block_size to pilot_sc_group_size %zu\n",
                     pilot_sc_group_size_);
      beam_block_size_ = pilot_sc_group_size_;
    }

    //Set beam block size to the pilot sc group size so events arn't generated for the redundant sc
    if ((beam_block_size_ % pilot_sc_group_size_) != 0) {
      AGORA_LOG_WARN(
          "beam_block_size(%zu) is not a multiple of pilot_sc_group_size(%zu). "
          "Efficiency will be decreased.  Please consider updating your "
          "settings\n",
          beam_block_size_, pilot_sc_group_size_);
    }
  }
  beam_events_per_symbol_ = 1 + (ofdm_data_num_ - 1) / beam_block_size_;

  fft_block_size_ = tdd_conf.value("fft_block_size", 1);
  fft_block_size_ = std::max(fft_block_size_, num_channels_);
  RtAssert(bs_ant_num_ % fft_block_size_ == 0,
           "FFT block size is set to an invalid value - all rx symbols per "
           "frame must fit inside an fft block");

  encode_block_size_ = tdd_conf.value("encode_block_size", 1);

  noise_level_ = tdd_conf.value("noise_level", 0.03);  // default: 30 dB
  AGORA_LOG_SYMBOL("Noise level: %.2f\n", noise_level_);

  // Scrambler and descrambler configurations
  scramble_enabled_ = tdd_conf.value("wlan_scrambler", true);

  // LDPC Coding and Modulation configurations
  ul_mcs_params_ = this->Parse(tdd_conf, "ul_mcs");

  dl_mcs_params_ = this->Parse(tdd_conf, "dl_mcs");



  fft_in_rru_ = tdd_conf.value("fft_in_rru", false);

  samps_per_symbol_ =
      ofdm_tx_zero_prefix_ + ofdm_ca_num_ + cp_len_ + ofdm_tx_zero_postfix_;
  packet_length_ =
      Packet::kOffsetOfData + ((kUse12BitIQ ? 3 : 4) * samps_per_symbol_);
  dl_packet_length_ = Packet::kOffsetOfData + (samps_per_symbol_ * 4);

  //Don't check for jumbo frames when using the hardware, this might be temp
  if (!kUseArgos) {
    RtAssert(packet_length_ < 9000,
             "Packet size must be smaller than jumbo frame");
  }

  this->running_.store(true);
  /* 12 bit samples x2 for I + Q */
  static const size_t kBitsPerSample = 12 * 2;
  const double bit_rate_mbps = (rate_ * kBitsPerSample) / 1e6;
  //For framer mode, we can ignore the Beacon
  //Double count the UlCal and DLCal to simplify things
  //Peak network traffic is the bit rate for 1 symbol, for non-hardware framer mode
  //the device can generate 2*rate_ traffic (for each tx symbol)
  const size_t bs_tx_symbols =
      frame_.NumDLSyms() + frame_.NumDLCalSyms() + frame_.NumULCalSyms();
  const size_t bs_rx_symbols = frame_.NumPilotSyms() + frame_.NumULSyms() +
                               frame_.NumDLCalSyms() + frame_.NumULCalSyms();
  const double per_bs_radio_traffic =
      ((static_cast<double>(bs_tx_symbols + bs_rx_symbols)) /
       frame_.NumTotalSyms()) *
      bit_rate_mbps;

  const size_t ue_tx_symbols = frame_.NumULSyms() + frame_.NumPilotSyms();
  //Rx all symbols, Tx the tx symbols (ul + pilots)
  const double per_ue_radio_traffic =
      (bit_rate_mbps *
       (static_cast<double>(ue_tx_symbols) / frame_.NumTotalSyms())) +
      bit_rate_mbps;

  AGORA_LOG_INFO(
      "Config: %zu BS antennas, %zu UE antennas, %zu pilot symbols per "
      "frame,\n"
      "\t%zu uplink data symbols per frame, %zu downlink data symbols "
      "per frame,\n"
      "\t%zu OFDM subcarriers (%zu data subcarriers),\n"
      "\tBeamforming %s, \n"
      //"\t%zu UL codeblocks per symbol, "
      
      "\tFrame time %.3f usec\n"
      "Radio Network Traffic Peak (Mbps): %.3f\n"
      "Radio Network Traffic Avg  (Mbps): %.3f\n"
      "Basestation Network Traffic Peak (Mbps): %.3f\n"
      "Basestation Network Traffic Avg  (Mbps): %.3f\n"
      "UE Network Traffic Peak (Mbps): %.3f\n"
      "UE Network Traffic Avg  (Mbps): %.3f\n"
      "All UEs Network Traffic Avg (Mbps): %.3f\n"
      "All UEs Network Traffic Avg (Mbps): %.3f\n",
      bs_ant_num_, ue_ant_num_, frame_.NumPilotSyms(), frame_.NumULSyms(),
      frame_.NumDLSyms(), ofdm_ca_num_, ofdm_data_num_, beamforming_str_.c_str(),
      this->GetFrameDurationSec() * 1e6,
      bit_rate_mbps, per_bs_radio_traffic, bit_rate_mbps * bs_ant_num_,
      per_bs_radio_traffic * bs_ant_num_, 2 * bit_rate_mbps,
      per_ue_radio_traffic, 2 * bit_rate_mbps * ue_ant_num_,
      per_ue_radio_traffic * ue_ant_num_);

  if (frame_.IsRecCalEnabled()) {
    AGORA_LOG_INFO(
        "Reciprical Calibration Enabled.  Full calibration data ready every "
        "%zu frame(s) using %zu symbols per frame\n",
        RecipCalFrameCnt(), frame_.NumDLCalSyms());
  }

  //Print();
}

json Config::Parse(const json& in_json, const std::string& json_handle) {
  json out_json;
  std::stringstream ss;
  ss << in_json.value(json_handle, out_json);
  out_json = json::parse(ss);
  if (out_json == nullptr) {
    out_json = json::object();
  }
  ss.str(std::string());
  ss.clear();
  return out_json;
}

Config::~Config() = default;

/* TODO Inspect and document */
size_t Config::GetSymbolId(size_t input_id) const {
  size_t symbol_id = SIZE_MAX;

  if (input_id < this->frame_.NumPilotSyms()) {
    symbol_id = this->Frame().GetPilotSymbol(input_id);
  } else {
    int new_idx = input_id - this->frame_.NumPilotSyms();

    // std::printf("\n*****GetSymbolId %d %zu\n", new_idx, input_id);
    if ((new_idx >= 0) &&
        (static_cast<size_t>(new_idx) < this->frame_.NumULSyms())) {
      symbol_id = this->Frame().GetULSymbol(new_idx);
    }
  }
  return symbol_id;
}

/* Returns True if symbol is valid index and is of symbol type 'P'
   False otherwise */
bool Config::IsPilot(size_t /*frame_id*/, size_t symbol_id) const {
  bool is_pilot = false;
  assert(symbol_id < this->frame_.NumTotalSyms());
  char s = frame_.FrameIdentifier().at(symbol_id);
#ifdef DEBUG3
  std::printf("IsPilot(%zu, %zu) = %c\n", frame_id, symbol_id, s);
#endif
  /* TODO should use the symbol type here */
  is_pilot = (s == 'P');
  return is_pilot;
}

/* Returns True if user equiptment and is a client dl pilot_
 * False otherwise */
bool Config::IsDlPilot(size_t /*frame_id*/, size_t symbol_id) const {
  bool is_pilot = false;
  assert(symbol_id < this->frame_.NumTotalSyms());
  char s = frame_.FrameIdentifier().at(symbol_id);
#ifdef DEBUG3
  std::printf("IsDlPilot(%zu, %zu) = %c\n", frame_id, symbol_id, s);
#endif
  if ((s == 'D') && (this->frame_.ClientDlPilotSymbols() > 0)) {
    size_t dl_index = this->frame_.GetDLSymbolIdx(symbol_id);
    is_pilot = (this->frame_.ClientDlPilotSymbols() > dl_index);
  }
  return is_pilot;
}

bool Config::IsCalDlPilot(size_t /*frame_id*/, size_t symbol_id) const {
  bool is_cal_dl_pilot = false;
  assert(symbol_id < this->frame_.NumTotalSyms());
  is_cal_dl_pilot = (this->frame_.FrameIdentifier().at(symbol_id) == 'C');
  return is_cal_dl_pilot;
}

bool Config::IsCalUlPilot(size_t /*frame_id*/, size_t symbol_id) const {
  bool is_cal_ul_pilot = false;
  assert(symbol_id < this->frame_.NumTotalSyms());
  is_cal_ul_pilot = (this->frame_.FrameIdentifier().at(symbol_id) == 'L');
  return is_cal_ul_pilot;
}

bool Config::IsUplink(size_t /*frame_id*/, size_t symbol_id) const {
  assert(symbol_id < this->frame_.NumTotalSyms());
  char s = frame_.FrameIdentifier().at(symbol_id);
#ifdef DEBUG3
  std::printf("IsUplink(%zu, %zu) = %c\n", frame_id, symbol_id, s);
#endif
  return (s == 'U');
}

bool Config::IsDownlinkBroadcast(size_t frame_id, size_t symbol_id) const {
  assert(symbol_id < this->frame_.NumTotalSyms());
  char s = frame_.FrameIdentifier().at(symbol_id);
#ifdef DEBUG3
  std::printf("IsDownlinkBroadcast(%zu, %zu) = %c\n", frame_id, symbol_id, s);
#else
  unused(frame_id);
#endif
  return (s == 'S');
}

bool Config::IsDownlink(size_t frame_id, size_t symbol_id) const {
  assert(symbol_id < this->frame_.NumTotalSyms());
  char s = frame_.FrameIdentifier().at(symbol_id);
#ifdef DEBUG3
  std::printf("IsDownlink(%zu, %zu) = %c\n", frame_id, symbol_id, s);
#else
  unused(frame_id);
#endif
  return (s == 'D');
}

SymbolType Config::GetSymbolType(size_t symbol_id) const {
  return kSymbolMap.at(this->frame_.FrameIdentifier().at(symbol_id));
}

//FrameStats Config::Frame() { return this->frame_; }

void Config::Print() const {
  if (kDebugPrintConfiguration == true) {
    std::cout << "Freq Ghz: " << freq_ghz_ << std::endl
              << "BaseStation ant num: " << bs_ant_num_ << std::endl
              << "BeamForming ant num: " << bf_ant_num_ << std::endl
              << "Ue num: " << ue_num_ << std::endl
              << "Ue ant num: " << ue_ant_num_ << std::endl
              << "Ue ant total: " << ue_ant_total_ << std::endl
              << "Ue ant offset: " << ue_ant_offset_ << std::endl
              << "OFDM Ca num: " << ofdm_ca_num_ << std::endl
              << "Cp Len: " << cp_len_ << std::endl
              << "Ofdm data num: " << ofdm_data_num_ << std::endl
              << "Ofdm data start: " << ofdm_data_start_ << std::endl
              << "Ofdm data stop: " << ofdm_data_stop_ << std::endl
              << "Ofdm pilot spacing: " << ofdm_pilot_spacing_ << std::endl
              << "Hardware framer: " << hw_framer_ << std::endl
              << "Ue Hardware framer: " << ue_hw_framer_ << std::endl
              << "Freq: " << freq_ << std::endl
              << "Rate: " << rate_ << std::endl
              << "NCO: " << nco_ << std::endl
              << "Scrambler Enabled: " << scramble_enabled_ << std::endl
              << "Radio Rf Freq: " << radio_rf_freq_ << std::endl
              << "Bw filter: " << bw_filter_ << std::endl
              << "Single Gain: " << single_gain_ << std::endl
              << "Tx Gain A: " << tx_gain_a_ << std::endl
              << "Rx Gain A: " << rx_gain_a_ << std::endl
              << "Tx Gain B: " << tx_gain_b_ << std::endl
              << "Rx Gain B: " << rx_gain_b_ << std::endl
              << "Calib Tx Gain A: " << calib_tx_gain_a_ << std::endl
              << "Calib Tx Gain B: " << calib_tx_gain_b_ << std::endl
              << "Num Cells: " << num_cells_ << std::endl
              << "Num Bs Radios: " << num_radios_ << std::endl
              << "Num Bs Channels: " << num_channels_ << std::endl
              << "Num Ue Channels: " << num_ue_channels_ << std::endl
              << "Beacon Ant: " << beacon_ant_ << std::endl
              << "Calib init repeat: " << init_calib_repeat_ << std::endl
              << "Beamsweep " << beamsweep_ << std::endl
              << "Sample Cal En: " << sample_cal_en_ << std::endl
              << "Imbalance Cal: " << imbalance_cal_en_ << std::endl
              << "Beamforming: " << beamforming_str_ << std::endl
              << "Bs Channel: " << channel_ << std::endl
              << "Ue Channel: " << ue_channel_ << std::endl
              << "Max Frames: " << frames_to_test_ << std::endl
              << "Transport Block Size: " << transport_block_size_ << std::endl
              << "Noise Level: " << noise_level_
              << std::endl
              << std::endl
              << "FFT in rru: " << fft_in_rru_ << std::endl;
  }
}
extern "C" {
__attribute__((visibility("default"))) Config* ConfigNew(char* filename) {
  auto* cfg = new Config(filename);
  //cfg->GenData();
  return cfg;
}
}