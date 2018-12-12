//
//  bmp_exporter.h
//
//  Created by 渡部 心 on 2014/07/28.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_bmp_exporter__
#define __SLR_bmp_exporter__

#include <cstdint>

void saveBMP(const char* filename, const void* pixels, uint32_t width, uint32_t height);

#endif /* __SLR_bmp_exporter__ */
