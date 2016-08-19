#include "SB.cxx"

#include <fstream>

void SB_all() {

  ofstream f("latexSB.txt");

  SB("HZZ4mu_ExtBrl4/hist-pu200.root", "ZZ4mu_ExtBrl4/hist-pu200_gamma.root", "", "ExtBrl4",&f);
  SB("HZZ4mu_IExtBrl4/hist-pu200.root", "ZZ4mu_IExtBrl4/hist-pu200.root", "", "IExtBrl4",&f);
  SB("HZZ4mu_InclBrl4/hist-pu200.root", "ZZ4mu_InclBrl4/hist-pu200_gamma.root", "", "InclBrl4",&f);

  SB("HZZ4mu_ExtBrl4_27/hist-pu200.root", "ZZ4mu_ExtBrl4_27/hist-pu200_gamma.root", "", "ExtBrl27",&f);
  SB("HZZ4mu_IExtBrl4_27/hist-pu200.root", "ZZ4mu_IExtBrl4_27/hist-pu200.root", "", "IExtBrl27",&f);
  SB("HZZ4mu_InclBrl4_27/hist-pu200.root", "ZZ4mu_InclBrl4_27/hist-pu200_gamma.root", "", "InclBrl27",&f);

  SB("HZZ4mu_ExtBrl4_32/hist-pu200.root", "ZZ4mu_ExtBrl4_32/hist-pu200_gamma.root", "", "ExtBrl32",&f);
  SB("HZZ4mu_IExtBrl4_32/hist-pu200.root", "ZZ4mu_IExtBrl4_32/hist-pu200.root", "", "IExtBrl32",&f);
  SB("HZZ4mu_InclBrl4_32/hist-pu200.root", "ZZ4mu_InclBrl4_32/hist-pu200_gamma.root", "", "InclBrl32",&f);


}
