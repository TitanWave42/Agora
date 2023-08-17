#include "ldpc_updater.h"

LdpcUpdater::LdpcUpdater(const OfdmConfig ofdm_data, uint16_t ul_base_graph, bool ul_early_term,
                         int16_t ul_max_decoder_iter, uint16_t dl_base_graph,
                         bool dl_early_term, int16_t dl_max_decoder_iter)
    : ofdm_data_num_ul_(ofdm_data.ofdm_data_num_ul_),
      ofdm_data_num_dl_(ofdm_data.ofdm_data_num_dl_),
      ofdm_ctrl_data_(ofdm_data.ofdm_ctrl_data_),
      ul_base_graph_(ul_base_graph),
      ul_early_term_(ul_early_term),
      ul_max_decoder_iter_(ul_max_decoder_iter),
      dl_base_graph_(dl_base_graph),
      dl_early_term_(dl_early_term),
      dl_max_decoder_iter_(dl_max_decoder_iter) {

        GenerateLdpcTable();
      }

LdpcUpdater::~LdpcUpdater() = default;

void LdpcUpdater::GenerateLdpcTable() {
  for (size_t mcs_index = 0; mcs_index < NumMcsIndices; mcs_index++) {
    std::pair<UlMcsParams*, DlMcsParams*> ldpc_info;
    UlMcsParams* ul_mcs_params;
    DlMcsParams* dl_mcs_params;

    ldpc_info.ul_mcs_params->ul_ldpc_config_ = UpdateUlLdpcConfig(mcs_index);
    CalculateUlLdpcProperties(ldpc_info.ul_mcs_params);
    ldpc_info.dl_mcs_params->dl_ldpc_config_ = UpdateDlLdpcConfig(mcs_index);
    CalculateDlLdpcProperties(ldpc_info.dl_mcs_params)
        ul_dl_ldpc_table.push_back(ldpc_info);
  }
}

LDPCconfig LdpcUpdater::UpdateUlLdpcConfig(size_t ul_mcs_index) {
  uint16_t base_graph = this->ul_base_graph_;

  size_t ul_mod_order_bits = GetModOrderBits(ul_mcs_index);
  size_t ul_code_rate = GetCodeRate(ul_mcs_index);

  size_t zc = SelectZc(base_graph, ul_code_rate, ul_mod_order_bits,
                       ofdm_data_num_ul_, Mcs::kCbPerSymbol, "uplink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows = static_cast<size_t>(std::round(
                        1024.0 * LdpcNumInputCols(base_graph) / ul_code_rate)) -
                    (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  LDPCconfig ul_ldpc_config_ =
      LDPCconfig(base_graph, zc, this->ul_max_decoder_iter, this->ul_early_term,
                 num_cb_len, num_cb_codew_len, num_rows, 0);

  ul_ldpc_config_.NumBlocksInSymbol((ofdm_data_num_ul_ * ul_mod_order_bits) /
                                    ul_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumULSyms() == 0) || (ul_ldpc_config_.NumBlocksInSymbol() > 0),
      "Uplink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");
  return ul_ldpc_config_;
}

LDPCconfig LdpcUpdater::UpdateDlLdpcConfig(size_t dl_mcs_index) {
  uint16_t base_graph = this->dl_base_graph_;

  size_t dl_mod_order_bits = GetModOrderBits(dl_mcs_index);
  size_t dl_code_rate = GetCodeRate(dl_mcs_index);

  size_t zc = SelectZc(base_graph, dl_code_rate, dl_mod_order_bits,
                       ofdm_data_num_dl_, Mcs::kCbPerSymbol, "downlink");

  // Always positive since ul_code_rate is smaller than 1024
  size_t num_rows = static_cast<size_t>(std::round(
                        1024.0 * LdpcNumInputCols(base_graph) / dl_code_rate)) -
                    (LdpcNumInputCols(base_graph) - 2);

  uint32_t num_cb_len = LdpcNumInputBits(base_graph, zc);
  uint32_t num_cb_codew_len = LdpcNumEncodedBits(base_graph, zc, num_rows);
  LDPCconfig dl_ldpc_config_ = LDPCconfig(
      base_graph, zc, this->dl_max_decoder_iter_, this->dl_early_term_,
      num_cb_len, num_cb_codew_len, num_rows, 0);

  dl_ldpc_config_.NumBlocksInSymbol((ofdm_data_num_dl_ * dl_mod_order_bits) /
                                    dl_ldpc_config_.NumCbCodewLen());
  RtAssert(
      (frame_.NumDLSyms() == 0) || (dl_ldpc_config_.NumBlocksInSymbol() > 0),
      "Downlink LDPC expansion factor is too large for number of OFDM data "
      "subcarriers.");

  return dl_ldpc_config_;
}

LDPCconfig LdpcUpdater::UpdateCtrlMCS(size_t control_mcs_index) {
  if (this->frame_.NumDlControlSyms() > 0) {
    const size_t dl_bcast_mcs_index = control_mcs_index;
    const size_t bcast_base_graph =
        1;  // TODO: For MCS < 5, base_graph 1 doesn't work
    dl_bcast_mod_order_bits_ = GetModOrderBits(dl_bcast_mcs_index);
    const size_t dl_bcast_code_rate = GetCodeRate(dl_bcast_mcs_index);
    std::string dl_bcast_modulation = MapModToStr(dl_bcast_mod_order_bits_);
    const int16_t max_decoder_iter = 5;
    size_t bcast_zc =
        SelectZc(bcast_base_graph, dl_bcast_code_rate, dl_bcast_mod_order_bits_,
                 ofdm_ctrl_data_, kCbPerSymbol, "downlink broadcast");

    // Always positive since dl_code_rate is smaller than 1
    size_t bcast_num_rows =
        static_cast<size_t>(std::round(
            1024.0 * LdpcNumInputCols(bcast_base_graph) / dl_bcast_code_rate)) -
        (LdpcNumInputCols(bcast_base_graph) - 2);

    uint32_t bcast_num_cb_len = LdpcNumInputBits(bcast_base_graph, bcast_zc);
    uint32_t bcast_num_cb_codew_len =
        LdpcNumEncodedBits(bcast_base_graph, bcast_zc, bcast_num_rows);
    LDPCconfig dl_bcast_ldpc_config_ =
        LDPCconfig(bcast_base_graph, bcast_zc, max_decoder_iter, true,
                   bcast_num_cb_len, bcast_num_cb_codew_len, bcast_num_rows, 0);

    dl_bcast_ldpc_config_.NumBlocksInSymbol(
        (ofdm_ctrl_data_ * dl_bcast_mod_order_bits_) /
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
    return dl_bcast_ldpc_config_;
  }
}

void LdpcUpdater::CalculateUlLdpcProperties(UlMcsParams* ul_mcs_params) {
  ul_mcs_params->ul_num_bytes_per_cb_ = ul_ldpc_config_.NumCbLen() / 8;
  ul_mcs_params->ul_num_padding_bytes_per_cb_ =
      Roundup<64>(ul_mcs_params->ul_num_bytes_per_cb_) -
      ul_mcs_params->ul_num_bytes_per_cb_;
  ul_mcs_params->ul_data_bytes_num_persymbol_ =
      ul_mcs_params->ul_num_bytes_per_cb_ * ul_mcs_params->NumBlocksInSymbol();
  uul_mcs_params->ul_mac_packet_length_ =
      ul_mcs_params->ul_data_bytes_num_persymbol_;

  const size_t ul_ldpc_input_min =
      (((ul_mcs_params->ul_ldpc_config_.NumCbLen() /
         ul_mcs_params->ul_ldpc_config_.ExpansionFactor()) -
        1) *
           (ul_mcs_params->ul_ldpc_config_.ExpansionFactor() / 8) +
       32);
  const size_t ul_ldpc_sugg_input = LdpcEncodingInputBufSize(
      ul_mcs_params->ul_ldpc_config_.BaseGraph(),
      ul_mcs_params->ul_ldpc_config_.ExpansionFactor());

  if (ul_ldpc_input_min > (ul_mcs_params->ul_num_bytes_per_cb_ +
                           ul_mcs_params->ul_num_padding_bytes_per_cb_)) {
    // Can cause a lot of wasted space, specifically the second argument of the max
    const size_t increased_padding =
        Roundup<64>(ul_ldpc_sugg_input) - ul_mcs_params->ul_num_bytes_per_cb_;

    AGORA_LOG_WARN(
        "LDPC required Input Buffer size exceeds uplink code block size!, "
        "Increased cb padding from %zu to %zu uplink CB Bytes %zu, LDPC "
        "Input Min for zc 64:256: %zu\n",
        ul_num_padding_bytes_per_cb_, increased_padding,
        ul_mcs_params->ul_num_bytes_per_cb_, ul_ldpc_input_min);
    ul_mcs_params->ul_num_padding_bytes_per_cb_ = increased_padding;
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

  return ul_ldpc_properties;
}

void LdpcUpdater::CalculateDlLdpcProperties(DlMcsParams* dl_mcs_params) {
  dl_mcs_params->dl_num_bytes_per_cb_ =
      dl_mcs_params->dl_ldpc_config_.NumCbLen() / 8;
  dl_mcs_params->dl_num_padding_bytes_per_cb_ =
      Roundup<64>(dl_mcs_params->dl_num_bytes_per_cb_) -
      dl_mcs_params->dl_num_bytes_per_cb_;
  dl_mcs_params->dl_data_bytes_num_persymbol_ =
      dl_mcs_params->dl_num_bytes_per_cb_ * dl_ldpc_config_.NumBlocksInSymbol();
  dl_mcs_params->dl_mac_packet_length_ =
      dl_mcs_params->dl_data_bytes_num_persymbol_;
  // Smallest over the air packet structure
  RtAssert(frame_.NumDLSyms() == 0 || dl_mcs_params->dl_mac_packet_length_ >
                                          sizeof(MacPacketHeaderPacked),
           "Downlink MAC Packet size must be larger than MAC header size");
  dl_mcs_params->dl_mac_data_length_max_ =
      dl_mcs_params->dl_mac_packet_length_ - sizeof(MacPacketHeaderPacked);

  dl_mcs_params->dl_mac_packets_perframe_ = frame_.NumDlDataSyms();
  dl_mcs_params->dl_mac_data_bytes_num_perframe_ =
      dl_mcs_params->dl_mac_data_length_max_ *
      dl_mcs_params->dl_mac_packets_perframe_;
  dl_mcs_params->dl_mac_bytes_num_perframe_ =
      dl_mcs_params->dl_mac_packet_length_ *
      dl_mcs_params->dl_mac_packets_perframe_;

  //((cb_len_bits / zc_size) - 1) * (zc_size / 8) + kProcBytes(32)
  const size_t dl_ldpc_input_min =
      (((dl_mcs_params->dl_ldpc_config_.NumCbLen() /
         dl_mcs_params->dl_ldpc_config_.ExpansionFactor()) -
        1) *
           (dl_mcs_params->dl_ldpc_config_.ExpansionFactor() / 8) +
       32);
  const size_t dl_ldpc_sugg_input = LdpcEncodingInputBufSize(
      dl_mcs_params->dl_ldpc_config_.BaseGraph(),
      dl_mcs_params->dl_ldpc_config_.ExpansionFactor());

  if (dl_ldpc_input_min > (dl_mcs_params->dl_num_bytes_per_cb_ +
                           dl_mcs_params->dl_num_padding_bytes_per_cb_)) {
    // Can cause a lot of wasted space, specifically the second argument of the max
    const size_t increased_padding =
        Roundup<64>(dl_ldpc_sugg_input) - dl_mcs_params->dl_num_bytes_per_cb_;

    AGORA_LOG_WARN(
        "LDPC required Input Buffer size exceeds downlink code block size!, "
        "Increased cb padding from %zu to %zu Downlink CB Bytes %zu, LDPC "
        "Input Min for zc 64:256: %zu\n",
        dl_mcs_params->dl_num_padding_bytes_per_cb_, increased_padding,
        dl_mcs_params->dl_num_bytes_per_cb_, dl_ldpc_input_min);
    dl_mcs_params->dl_num_padding_bytes_per_cb_ = increased_padding;
  }
}