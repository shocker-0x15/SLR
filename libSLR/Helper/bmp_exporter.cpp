//
//  BMPExporter.cpp
//  OpenCL_PathTracer
//
//  Created by 渡部 心 on 2014/07/29.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#include "bmp_exporter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void saveBMP(const char* filename, const void* pixels, uint32_t width, uint32_t height) {
    uint32_t byteWidth = 3 * width + width % 4;
    const uint32_t FILE_HEADER_SIZE = 14;
    const uint32_t INFO_HEADER_SIZE = 40;
    const uint32_t HEADER_SIZE = FILE_HEADER_SIZE + INFO_HEADER_SIZE;
    uint8_t header[HEADER_SIZE];
    uint32_t fileSize = byteWidth * height + HEADER_SIZE;
    uint32_t dataOffset = HEADER_SIZE;
    uint16_t numPlanes = 1;
    uint16_t numBits = 24;
    uint32_t compress = 0;
    uint32_t dataSize = byteWidth * height;
    uint32_t xppm = 1;
    uint32_t yppm = 1;
    
    FILE* fp = fopen(filename, "wb");
    if (fp == nullptr)
        return;
    
    memset(header, 0, HEADER_SIZE);
    header[0] = 'B';
    header[1] = 'M';
    memcpy(header + 2, &fileSize, sizeof(uint32_t));
    memcpy(header + 10, &dataOffset, sizeof(uint32_t));
    
    memcpy(header + 14, &INFO_HEADER_SIZE, sizeof(uint32_t));
    memcpy(header + 18, &width, sizeof(uint32_t));
    memcpy(header + 22, &height, sizeof(uint32_t));
    memcpy(header + 26, &numPlanes, sizeof(uint16_t));
    memcpy(header + 28, &numBits, sizeof(uint16_t));
    memcpy(header + 30, &compress, sizeof(uint32_t));
    memcpy(header + 34, &dataSize, sizeof(uint32_t));
    memcpy(header + 38, &xppm, sizeof(uint32_t));
    memcpy(header + 42, &yppm, sizeof(uint32_t));
    
    fwrite(header, sizeof(uint8_t), HEADER_SIZE, fp);
    fwrite(pixels, sizeof(uint8_t), dataSize, fp);
    
    fclose(fp);
}
