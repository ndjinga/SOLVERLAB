#ifndef __CDMATHEXCEPTION_HXX__
#define __CDMATHEXCEPTION_HXX__

#include <string>
#include <exception>


class CdmathException : public std::exception
{
	public:
	CdmathException(std::string reason);
	CdmathException(std::string reason, std::string file, int line);
    ~CdmathException() throw ();
    const char *what() const throw();
  protected:
    std::string _reason;
  };

#endif
