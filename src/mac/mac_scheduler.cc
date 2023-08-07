/**
 * @file  mac_scheduler.cc
 * @brief Declaration file for the simple MAC scheduler
 */
#include "mac_scheduler.h"

#include "comms-lib.h"

MacScheduler::MacScheduler(Config* const cfg)
    : cfg_(cfg), mcs_(std::make_unique<Mcs>(cfg)) {
  //mcs_ = std::make_unique<Mcs>(cfg);
  num_groups_ =
      (cfg_->SpatialStreamsNum() == cfg_->UeAntNum()) ? 1 : cfg_->UeAntNum();
  schedule_buffer_.Calloc(num_groups_, cfg_->UeAntNum() * cfg_->OfdmDataNum(),
                          Agora_memory::Alignment_t::kAlign64);
  schedule_buffer_index_.Calloc(num_groups_,
                                cfg_->SpatialStreamsNum() * cfg_->OfdmDataNum(),
                                Agora_memory::Alignment_t::kAlign64);
  // Create round-robin schedule
  for (size_t gp = 0u; gp < num_groups_; gp++) {
    for (size_t sc = 0; sc < cfg_->OfdmDataNum(); sc++) {
      for (size_t ue = gp; ue < gp + cfg_->SpatialStreamsNum(); ue++) {
        size_t cur_ue = ue % cfg_->UeAntNum();
        // for now all SCs are allocated to scheduled UEs
        schedule_buffer_[gp][cur_ue + cfg_->UeAntNum() * sc] = 1;
        schedule_buffer_index_[gp][(ue - gp) + cfg_->SpatialStreamsNum() * sc] =
            cur_ue;
      }
    }
  }

  ul_mcs_buffer_.Calloc(num_groups_, cfg_->UeAntNum(),
                        Agora_memory::Alignment_t::kAlign64);
  dl_mcs_buffer_.Calloc(num_groups_, cfg_->UeAntNum(),
                        Agora_memory::Alignment_t::kAlign64);
  for (size_t gp = 0u; gp < num_groups_; gp++) {
    for (size_t ue = 0; ue < cfg_->UeAntNum(); ue++) {
      ul_mcs_buffer_[gp][ue] = mcs_->McsIndex(Direction::kUplink);
      dl_mcs_buffer_[gp][ue] = mcs_->McsIndex(Direction::kDownlink);
    }
  }
}

MacScheduler::~MacScheduler() {
  schedule_buffer_.Free();
  schedule_buffer_index_.Free();
  ul_mcs_buffer_.Free();
  dl_mcs_buffer_.Free();
}

void MacScheduler::CheckDlMcs(float snr, size_t frame_id, size_t ant_id) {
  float spectral_efficiency = this->mcs_->SpectralEffeciency(snr);
  std::cout << "checking the dl mcs" << std::endl << std::flush;
  if (abs(kmcs_index_to_spectral_effeciency.at(mcs_->CurrentDlMcs().mcs_index) -
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
    mcs_->SetNextDlMcs(frame_id, GetModOrderBits(optimal_mcs_index));
  }
}

bool MacScheduler::IsUeScheduled(size_t frame_id, size_t sc_id, size_t ue_id) {
  size_t gp = frame_id % num_groups_;
  return (schedule_buffer_[gp][ue_id + cfg_->UeAntNum() * sc_id] != 0);
}

size_t MacScheduler::ScheduledUeIndex(size_t frame_id, size_t sc_id,
                                      size_t sched_ue_id) {
  return (size_t)this->ScheduledUeList(frame_id, sc_id)[sched_ue_id];
}

arma::uvec MacScheduler::ScheduledUeMap(size_t frame_id, size_t sc_id) {
  size_t gp = frame_id % num_groups_;
  return arma::uvec(reinterpret_cast<unsigned long long*>(
                        &schedule_buffer_[gp][cfg_->UeAntNum() * sc_id]),
                    cfg_->UeAntNum(), false);
}

arma::uvec MacScheduler::ScheduledUeList(size_t frame_id, size_t sc_id) {
  size_t gp = frame_id % num_groups_;
  return sort(arma::uvec(
      reinterpret_cast<unsigned long long*>(
          &schedule_buffer_index_[gp][cfg_->SpatialStreamsNum() * sc_id]),
      cfg_->SpatialStreamsNum(), false));
}

size_t MacScheduler::ScheduledUeUlMcs(size_t frame_id, size_t ue_id) {
  size_t gp = frame_id % num_groups_;
  return ul_mcs_buffer_[gp][ue_id];
}

size_t MacScheduler::ScheduledUeDlMcs(size_t frame_id, size_t ue_id) {
  size_t gp = frame_id % num_groups_;
  return dl_mcs_buffer_[gp][ue_id];
}
