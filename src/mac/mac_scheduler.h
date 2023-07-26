/**
 * @file mac_scheduler.h
 * @brief Declaration file for the simple MAC scheduler class
 */
#ifndef MAC_SCHEDULER_H_
#define MAC_SCHEDULER_H_

#include "config.h"
#include "mcs.h"

class MacScheduler {
 public:
  explicit MacScheduler(Config* const cfg);
  ~MacScheduler();

  bool IsUeScheduled(size_t frame_id, size_t sc_id, size_t ue_id);
  size_t ScheduledUeIndex(size_t frame_id, size_t sc_id, size_t sched_ue_id);
  arma::uvec ScheduledUeList(size_t frame_id, size_t sc_id);
  arma::uvec ScheduledUeMap(size_t frame_id, size_t sc_id);
  size_t ScheduledUeUlMcs(size_t frame_id, size_t ue_id);
  size_t ScheduledUeDlMcs(size_t frame_id, size_t ue_id);

  inline Config* Cfg() { return this->cfg_; }

 private:
  size_t num_groups_;
  Table<int> schedule_buffer_;
  Table<size_t> schedule_buffer_index_;
  Table<size_t> ul_mcs_buffer_;
  Table<size_t> dl_mcs_buffer_;
  Config* const cfg_;
  Mcs* const mcs_;
};

#endif  // MAC_SCHEDULER_H_
