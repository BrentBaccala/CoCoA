#include <iostream>
#include <string>
#include "OnlineHelp.H"

std::string packageDir = "packages";

int main()
{
  CoCoA::OnlineHelp::PrintWordlist(std::cout);
}
