/**
 * @file mcs.h
 * @brief Class definition for mcs.
 */
#ifndef MCS_H_
#define MCS_H_

#include "config.h"
#include "framestats.h"
#include "ldpc_config.h"
#include "ldpc_generator.h"
#include "memory_manage.h"
#include "nlohmann/json.hpp"
#include "utils.h"

//#include "comms-constants.inc"

constexpr size_t kNumTables = 4;

//For determining if we need to change mcs
constexpr float kSnrDeltaTolerance = 0.1;

static const std::map<size_t, float> kMcsIndexToSpectralEffeciency = {
    {0, 0.2344},  {1, 0.3066},  {2, 0.3770},  {3, 0.4902},  {4, 0.6016},
    {5, 0.7402},  {6, 0.8770},  {7, 1.0273},  {8, 1.1758},  {9, 1.3262},
    {10, 1.3281}, {11, 1.4766}, {12, 1.6953}, {13, 1.9141}, {14, 2.1602},
    {15, 2.4063}, {16, 2.5703}, {17, 2.5664}, {18, 2.7305}, {19, 3.0293},
    {20, 3.3223}, {21, 3.6094}, {22, 3.9023}, {23, 4.2129}, {24, 4.5234},
    {25, 4.8164}, {26, 5.1152}, {27, 5.3320}, {28, 5.5547}, {29, 5.8906},
    {30, 6.2266}, {31, 6.5703}};

struct McsScheme {
  size_t frame_number_;
  size_t mcs_index_;
  size_t mod_order_bits_;
  size_t code_rate_;
};

struct UserMcs {
  McsScheme current_ul_mcs_;
  McsScheme current_dl_mcs_;
  McsScheme next_ul_mcs_;
  McsScheme next_dl_mcs_;
}

struct ModulationTables {
  Table<complex_float> ul_tables_[kNumTables];
  Table<complex_float> dl_tables_[kNumTables];
};

class Mcs {
 public:
  explicit Mcs();
  ~Mcs();

  void CheckUlMcs(float snr, size_t frame_id);

  //Example accessors in the mcs for values that are requested from the config.
  inline LDPCconfig ul_ldpc_config_(size_t mcs_index) {
    return this->ldpc_generator_.LdpcTable().at(mcs_index).first;
  }
  inline LDPCconfig dl_ldpc_config_(size_t mcs_index) {
    return this->ldpc_generator_.LdpcTable().at(mcs_index).second;
  }

 private:
  std::string ul_modulation_;
  std::string dl_modulation_;
  FrameStats frame_;
  std::unique_ptr<LdpcGenerator> ldpc_generator_;

  void UpdateUlMcs(size_t current_frame_number);
  void UpdateDlMcs(size_t current_frame_number);
  void CalculateLdpcProperties();
  void CreateModulationTables();

  void SetNextUlMcs(size_t frame_number, size_t mod_order_bits);
  void SetNextDlMcs(size_t frame_number, size_t mod_order_bits);
  void InitializeUlMcs(const nlohmann::json ul_mcs);
  void InitializeDlMcs(const nlohmann::json dl_mcs);
  void UpdateUlLdpcConfig();
  void UpdateDlLdpcConfig();

  void UpdateMcs(size_t current_frame_number);
  void DumpMcsInfo();
  void UpdateCtrlMCS();
};

#endif /* MCS_H_ */
