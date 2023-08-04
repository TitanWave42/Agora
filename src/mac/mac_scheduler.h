/**
 * @file mac_scheduler.h
 * @brief Declaration file for the simple MAC scheduler class
 */
#ifndef MAC_SCHEDULER_H_
#define MAC_SCHEDULER_H_

#include "mcs.h"

class MacScheduler {
 public:
  explicit MacScheduler(Config* const cfg);
  ~MacScheduler();

  inline void Running(bool value) { this->mcs_->Running(value); }
  inline bool Running() const { return mcs_->Running(); }

  inline Table<complex_float>& UeSpecificPilot() {
    return this->mcs_->UeSpecificPilot();
  }
  inline Table<std::complex<int16_t>>& UeSpecificPilotT() {
    return this->mcs_->UeSpecificPilotT();
  }


  inline void InitializeUlMcs(const nlohmann::json ul_mcs) {
    this->mcs_->InitializeUlMcs(ul_mcs);
  }
  inline void InitializeDlMcs(const nlohmann::json dl_mcs) {
    this->mcs_->InitializeDlMcs(dl_mcs);
  }

  //inline size_t Mcs::MacBytesNumPerframe(Direction dir) const
  inline size_t DecodeBroadcastSlots(const int16_t* const bcast_iq_samps) {
    return this->mcs_->DecodeBroadcastSlots(bcast_iq_samps);
  }
  inline void GenBroadcastSlots(
      std::vector<std::complex<int16_t>*>& bcast_iq_samps,
      std::vector<size_t> ctrl_msg) {
    this->mcs_->GenBroadcastSlots(bcast_iq_samps, ctrl_msg);
  }
  inline void GenData() { this->mcs_->GenData(); }
  inline void GenPilots() { this->mcs_->GenPilots(); }
  inline Table<int8_t>& DlBits() { return mcs_->DlBits(); }
  inline Table<int8_t>& UlBits() { return mcs_->UlBits(); }

  inline Config* Cfg() { return this->cfg_; }
  //inline Mcs* GetMcs() { return this->mcs_; }

  inline float Scale() const { return this->mcs_->Scale(); }
  inline size_t MacBytesNumPerframe(Direction dir) const {
    return mcs_->MacBytesNumPerframe(dir);
  }

  inline const LDPCconfig& LdpcConfig(Direction dir) const {
    return mcs_->LdpcConfig(dir);
  }
  inline size_t NumBytesPerCb(Direction dir) const {
    return mcs_->NumBytesPerCb(dir);
  }

  inline size_t ModOrderBits(Direction dir) const {
    return mcs_->ModOrderBits(dir);
  }

  inline size_t NumPaddingBytesPerCb(Direction dir) const {
    return mcs_->NumPaddingBytesPerCb(dir);
  }

  inline size_t MacPacketLength(Direction dir) const {
    return mcs_->MacPacketLength(dir);
  }

  inline size_t MacPayloadMaxLength(Direction dir) const {
    return mcs_->MacPayloadMaxLength(dir);
  }

  inline size_t SubcarrierPerCodeBlock(Direction dir) const {
    return mcs_->SubcarrierPerCodeBlock(dir);
  }
  inline size_t MacDataBytesNumPerframe(Direction dir) const {
    return mcs_->MacBytesNumPerframe(dir);
  }
  inline size_t MacPacketsPerframe(Direction dir) const {
    return mcs_->MacPacketsPerframe(dir);
  }
  inline size_t McsIndex(Direction dir) const { return mcs_->McsIndex(dir); }
  inline const complex_float* PilotsSgn() const { return mcs_->PilotsSgn(); };
  inline Table<complex_float>& UlIqF() { return mcs_->UlIqF(); }
  inline Table<complex_float>& DlIqF() { return mcs_->DlIqF(); }
  inline Table<std::complex<int16_t>>& DlIqT() { return mcs_->DlIqT(); }
  inline Table<complex_float>& ModTable(Direction dir) {
    return mcs_->ModTable(dir);
  }
  inline const std::vector<std::complex<float>>& PilotCf32() const {
    return mcs_->PilotCf32();
  }
  inline std::vector<std::complex<int16_t>>& PilotUeCi16(size_t ue_id,
                                                         size_t pilot_idx) {
    return mcs_->PilotUeCi16(ue_id, pilot_idx);
  }
  inline const nlohmann::json& McsParams(Direction dir) const {
    return mcs_->McsParams(dir);
  }
  inline Table<int8_t>& UlModBits() { return mcs_->UlModBits(); }
  inline Table<int8_t>& DlModBits() { return mcs_->DlModBits(); }
  inline const arma::uvec& PilotUeSc(size_t ue_id) const {
    return mcs_->PilotUeSc(ue_id);
  }
  //inline void GenData() { mcs_->GenData(); }
  inline std::string Modulation(Direction dir) const {
    return mcs_->Modulation(dir);
  }
  /// Get info bits for this symbol, user and code block ID
  inline int8_t* GetInfoBits(Table<int8_t>& info_bits, Direction dir,
                             size_t symbol_id, size_t ue_id,
                             size_t cb_id) const {
    return mcs_->GetInfoBits(info_bits, dir, symbol_id, ue_id, cb_id);
  }
  /// Get encoded_buffer for this frame, symbol, user and code block ID
  inline int8_t* GetModBitsBuf(Table<int8_t>& mod_bits_buffer, Direction dir,
                               size_t frame_id, size_t symbol_id, size_t ue_id,
                               size_t sc_id) const {
    return mcs_->GetModBitsBuf(mod_bits_buffer, dir, frame_id, symbol_id, ue_id,
                               sc_id);
  }
  //EVENTUALLY UPDATE THIS FUNCTION TO TAKE num_byte_per_cb as an argument
  //because the num_bytes_per_cb is determined by the MCS and is effectively
  //How much data we can read from the date buffer.
  /// Get info bits for this symbol, user and code block ID
  /// Get mac bits for this frame, symbol, user and code block ID
  inline int8_t* GetMacBits(Table<int8_t>& info_bits, Direction dir,
                            size_t frame_id, size_t symbol_id, size_t ue_id,
                            size_t cb_id) const {
    return mcs_->GetMacBits(info_bits, dir, frame_id, symbol_id, ue_id, cb_id);
  }

  bool IsUeScheduled(size_t frame_id, size_t sc_id, size_t ue_id);
  size_t ScheduledUeIndex(size_t frame_id, size_t sc_id, size_t sched_ue_id);
  arma::uvec ScheduledUeList(size_t frame_id, size_t sc_id);
  arma::uvec ScheduledUeMap(size_t frame_id, size_t sc_id);
  size_t ScheduledUeUlMcs(size_t frame_id, size_t ue_id);
  size_t ScheduledUeDlMcs(size_t frame_id, size_t ue_id);

 private:
  size_t num_groups_;
  Table<int> schedule_buffer_;
  Table<size_t> schedule_buffer_index_;
  Table<size_t> ul_mcs_buffer_;
  Table<size_t> dl_mcs_buffer_;
  Config* const cfg_;
  std::unique_ptr<Mcs> mcs_;
};

#endif  // MAC_SCHEDULER_H_
