#ifndef PRINTMESSAGE_H
#define PRINTMESSAGE_H

#include <TString.h>

enum class Verbosity {
  VERBOSE,
    INFO,
    DEBUG,
    WARNING,
    ERROR,
    FATAL};

void PrintMessage(const TString&, const Verbosity&, const Verbosity&);

#endif //PRINTMESSAGE_H
