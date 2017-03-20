#include <stdarg.h>
#include "haanju_misc.hpp"

#define HJ_DEBUG

/************************************************************************
 Method Name: MakeTrackIDList
 Description:
	- make a string filled with track IDs in the input track set
 Input Arguments:
	- tracks: target track set
 Return Values:
	- string filled with track IDs.
************************************************************************/
std::string hj::MakeTrackIDList(hj::TrackSet *tracks)
{
	std::string strResult("{");
	for (TrackSet::iterator trackIter = tracks->begin();
		trackIter != tracks->end();
		trackIter++)
	{
		strResult = strResult + std::to_string((*trackIter)->id);
		if (trackIter < tracks->end() - 1) { strResult += ","; }
	}
	strResult += "}";

	return strResult;
}


/************************************************************************
 Method Name: printf_debug
 Description:
	- print debugging messages to the console when 'HJ_DEBUG' is defined
 Input Arguments:
	- format: formatted string
	- ...   : arguments in formatted string
 Return Values:
	- none
************************************************************************/
void hj::printf_debug(const char *_format, ...)
{
#ifdef HJ_DEBUG
	va_list args;
	va_start(args, _format);
	vprintf(_format, args);
	va_end(args);
#endif
}


//()()
//('')HAANJU.YOO


