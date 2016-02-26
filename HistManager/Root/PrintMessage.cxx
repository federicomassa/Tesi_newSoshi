#include "HistManager/PrintMessage.h"
#include <iostream>

void PrintMessage(const TString& message, const Verbosity& message_verbosity, const Verbosity& verbosity_level) {
  const int& message_verbosity_int = static_cast<int>(message_verbosity);
  const int& verbosity_level_int = static_cast<int>(verbosity_level);

  if (message_verbosity_int >= verbosity_level_int)
    std::cout << message << std::endl;

}
