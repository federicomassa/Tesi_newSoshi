#include <iostream>
#include <string>

void CheckPtr(void* ptr) {
  try {
    if (ptr == NULL) {
      const std::string error = "ERROR:\t null pointer exception";
      throw error;
    }
  }
  catch(std::string& e) {
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
}
