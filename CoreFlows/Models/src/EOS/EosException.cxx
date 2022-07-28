#include "EosException.hxx"

using namespace std;

EosException::EosException(std::string reason):_reason(reason)
{
}

EosException::EosException(std::string reason, std::string file, int line):_reason(reason)
{
}

EosException::~EosException() throw ()
{
}

const char *EosException::what() const throw()
{
  return _reason.c_str();
}
