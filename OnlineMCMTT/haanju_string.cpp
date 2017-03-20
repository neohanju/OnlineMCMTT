#include <sstream>
#include <vector>
#include <memory>    // For std::unique_ptr
#include <stdarg.h>  // For va_start, etc.
#include "haanju_string.hpp"

namespace hj
{

std::string FormattedString(const std::string fmt_str, ...)
{
	int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
	std::string str;
	std::unique_ptr<char[]> formatted;
	va_list ap;
	while (1)
	{
		formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
		strcpy_s(formatted.get(), n, fmt_str.c_str());
		va_start(ap, fmt_str);
		final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
		va_end(ap);
		if (final_n < 0 || final_n >= n)
		{
			n += abs(final_n - n + 1);
		}
		else
		{
			break;
		}
	}
	return std::string(formatted.get());
}

long GetTimeFromFileName(
	const std::string _frameFileName, 
	const std::string _strPrefix, 
	const std::string _strPostfix,
	const std::string _strFormat)
{
	std::string strTime = _frameFileName;	
	long readTime     = -1;
	
	// remove prefix
	size_t foundPos = strTime.find(_strPrefix);
	if (std::string::npos == foundPos) { return readTime; }
	strTime.erase(foundPos, _strPrefix.length());

	// remove postfix
	foundPos = strTime.find(_strPostfix);
	if (std::string::npos == foundPos) { return readTime; }
	strTime.erase(foundPos, _strPostfix.length());

	// convert string to integer
	if (0 == _strFormat.compare("MMDDHHMMSSmmm"))
	{
		strTime.erase(strTime.begin(), strTime.begin() + 4);
		readTime = std::stol(strTime);
	}

	return readTime;
}

}
