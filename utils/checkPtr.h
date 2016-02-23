#include <iostream>

void CheckPtr(void* ptr) {
  if (ptr == 0) {
    std::cout << "ERROR: null pointer exception" << std::endl;
    exit(EXIT_FAILURE);
  }
}
