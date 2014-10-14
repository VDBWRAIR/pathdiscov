// regex_match example
#include <iostream>
#include <string>
#include <regex>

int main ()
{

  if (std::regex_match ("subject", std::regex("(sub)(.*)") ))
    std::cout << "string literal matched\n";

  std::string s ("subject");
  std::regex e ("(sub)(.*)");
  if (std::regex_match (s,e))
    std::cout << "string object matched\n";

  return 0;
}
