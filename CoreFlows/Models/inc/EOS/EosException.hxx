#ifndef __EOSEXCEPTION_HXX__
#define __EOSEXCEPTION_HXX__

#include <string>
#include <exception>


class EosException : public std::exception
{
	public:
	EosException(std::string reason);
	EosException(std::string reason, std::string file, int line);
    ~EosException() throw ();
    const char *what() const throw();
  protected:
    std::string _reason;
  };

#endif
