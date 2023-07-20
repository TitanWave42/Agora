/**
 * @file mcs.h
 * @brief Class definition for mcs.
 */
#ifndef MCS_H_
#define MCS_H_

#include "utils.h"


 struct Modulation_Scheme {
    size_t frame_number;
    size_t modulation_type;
  };



class Mcs {

  public:
    explicit Mcs(size_t initial_modulation);
    ~Mcs();

    void Create_Modulation_Tables();
    void Update_Modulation_Scheme();
    void Set_Next_Modulation_Scheme();


  private:
    Modulation_Scheme current_modulation_;
    Modulation_Scheme next_modulation_;

}


#endif  // MCS_H_
