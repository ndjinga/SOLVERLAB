#include "CdmathException.hxx"

using namespace std;

CdmathException::CdmathException(std::string reason):_reason(reason)
{
}

CdmathException::CdmathException(std::string reason, std::string file, int line):_reason(reason)
{
}

CdmathException::~CdmathException() throw ()
{
}

const char *CdmathException::what() const throw()
{
  return _reason.c_str();
}
