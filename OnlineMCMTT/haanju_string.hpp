/******************************************************************************
* Title        : haanju_string
* Author       : Haanju Yoo
* Initial Date : 2016.08.28 (ver. 0.9)
* Version Num. : 0.9 
* Description  : contains string related operation.
******************************************************************************
*               .__                           __.
*                \ `\~~---..---~~~~~~--.---~~| /
*                 `~-.   `                   .~         _____
*                     ~.                .--~~    .---~~~    /
*                      / .-.      .-.      |  <~~        __/
*                     |  |_|      |_|       \  \     .--'
*                    /-.      -       .-.    |  \_   \_
*                    \-'   -..-..-    `-'    |    \__  \_
*                     `.                     |     _/  _/
*                     ~-                .,-\   _/  _/
*                      /                 -~~~~\ /_  /_
*                     |               /   |    \  \_  \_
*                     |   /          /   /      | _/  _/
*                     |  |          |   /    .,-|/  _/
*                     )__/           \_/    -~~~| _/
*                       \                      /  \
*                        |           |        /_---`
*                        \    .______|      ./
*                        (   /        \    /
*                        `--'          /__/
*
******************************************************************************/

#ifndef __HAANJU_STRING_HPP__
#define __HAANJU_STRING_HPP__

#include <sstream>
#include <vector>
#include <memory>    // For std::unique_ptr
#include <stdarg.h>  // For va_start, etc.

namespace hj
{

std::string FormattedString(const std::string fmt_str, ...);
long GetTimeFromFileName(
	const std::string _frameFileName, 
	const std::string _strPrefix, 
	const std::string _strPostfix,
	const std::string _strFormat = "MMDDHHMMSSmmm");

template<typename _Tp>
int StringDelimite(std::string inputString, const char deliminator, std::vector<_Tp> &resultVector)
{
	resultVector.clear();
	std::stringstream ss(inputString);

	_Tp element;
	while (ss >> element)
	{
		resultVector.push_back(element);
		if (deliminator == ss.peek()) { ss.ignore(); }
	}
	return (int)resultVector.size();
}

}


#endif

//()()
//('')HAANJU.YOO


