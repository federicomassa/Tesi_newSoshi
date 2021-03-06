#include <string>
#include <iostream>

void Assert(const std::string& message, bool condition) {
  try {
    if (condition) return;
    else {
      const std::string errStr = "ERROR:\t";
      throw (errStr + message);
    }
  }
  catch (std::string& e) {
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
}

