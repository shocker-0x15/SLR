//  debugPrintf.cpp
//
//  Created by 渡部 心 on 2016/08/08.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "defines.h"

#ifdef SLR_Defs_MSVC
SLR_API void debugPrintf(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	char str[1024];
	vsprintf(str, fmt, args);
	va_end(args);
	OutputDebugString(str);
}
#endif
