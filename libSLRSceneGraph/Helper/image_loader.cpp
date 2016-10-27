//
//  ImageLoader.cpp
//
//  Created by 渡部 心 on 2014/05/05.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#include "image_loader.h"
#include <libpng16/png.h>
#include <ImfInputFile.h>
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <cstdlib>
#include <string>
#include <cassert>

enum EXRType {
    Plane = 0,
    LatitudeLongitude,
    CubeMap
};

bool getEXRInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize, EXRType* exrtype) {
    using namespace Imf;
    using namespace Imath;
    RgbaInputFile file(filePath.c_str());
    Imf::Header header = file.header();
    
    Box2i dw = file.dataWindow();
    *width = dw.max.x - dw.min.x + 1;
    *height = dw.max.y - dw.min.y + 1;
    
    *requiredSize = *height * *width * sizeof(Rgba);
    
    return true;
}

bool loadEXR(const std::string &filePath, uint8_t* storage) {
    using namespace Imf;
    using namespace Imath;
    RgbaInputFile file(filePath.c_str());
    Imf::Header header = file.header();
    
    Box2i dw = file.dataWindow();
    long width = dw.max.x - dw.min.x + 1;
    long height = dw.max.y - dw.min.y + 1;
    Array2D<Rgba> pixels{height, width};
    pixels.resizeErase(height, width);
    file.setFrameBuffer(&pixels[0][0] - dw.min.x - dw.min.y * width, 1, width);
    file.readPixels(dw.min.y, dw.max.y);
    
    long rowSize = width * sizeof(Rgba);
    uint8_t* orgDataHead = storage;
    uint8_t* curDataHead = orgDataHead;
    for (int i = 0; i < height; ++i) {
        memcpy(curDataHead, pixels[i], rowSize);
        curDataHead += rowSize;
    }
    
    return true;
}


bool getJPEGInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize) {
    return false;
}

bool loadJPEG(const std::string &filePath, uint8_t* storage) {
    return false;
}


struct PNGStructures {
    png_structp pngStruct;
    png_infop pngInfo;
    png_infop pngEndInfo;
};

int libPNGreadChunkCallback(png_structp pngStruct, png_unknown_chunkp chunk) {
    return 0;
}

void libPNGreadRowCallback(png_structp pngStruct, png_uint_32 row, int pass) {
    
}

bool preloadPNG(const std::string &filePath, FILE** fpRet, PNGStructures* pngStructures) {
    FILE* &fp = *fpRet;
    png_structp &pngStruct = pngStructures->pngStruct;
    png_infop &pngInfo = pngStructures->pngInfo;
    png_infop &pngEndInfo = pngStructures->pngEndInfo;
    
    fp = fopen(filePath.c_str(), "rb");
    if (fp == nullptr) {
        printf("Failed to open the file.\n%s\n", filePath.c_str());
        return false;
    }
    
    uint8_t headSignature[8];
    fread(headSignature, 1, sizeof(headSignature), fp);
    if (png_sig_cmp(headSignature, 0, sizeof(headSignature))) {
        fclose(fp);
        printf("This is not a png file.\n");
        return false;
    }
    
    pngStruct = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!pngStruct) {
        fclose(fp);
        printf("libpng: Failed to create a read struct.\n");
        return false;
    }
    
    // 画像データ前のチャンク用の構造体を確保。
    pngInfo = png_create_info_struct(pngStruct);
    if (!pngInfo) {
        png_destroy_read_struct(&pngStruct, nullptr, nullptr);
        fclose(fp);
        printf("libpng: Failed to create an info struct.\n");
        return false;
    }
    
    // 画像データ後のチャンク用の構造体を確保。
    pngEndInfo = png_create_info_struct(pngStruct);
    if (!pngEndInfo) {
        png_destroy_read_struct(&pngStruct, &pngInfo, nullptr);
        fclose(fp);
        printf("libpng: Failed to create an end info struct.\n");
        return false;
    }
    
    if (setjmp(png_jmpbuf(pngStruct))) {
        png_destroy_read_struct(&pngStruct, &pngInfo, &pngEndInfo);
        fclose(fp);
        printf("libpng: Error.\n");
        return false;
    }
    
    png_init_io(pngStruct, fp);
    png_set_sig_bytes(pngStruct, sizeof(headSignature));
    png_set_read_user_chunk_fn(pngStruct, png_get_user_chunk_ptr(pngStruct), libPNGreadChunkCallback);
    //    png_set_read_status_fn(pngStruct, libPNGreadRowCallback);
//#ifdef PNG_UNKNOWN_CHUNKS_SUPPORTED
//    png_byte chunkList[0] = {};
//    png_set_keep_unknown_chunks(pngStruct, PNG_HANDLE_CHUNK_AS_DEFAULT, chunkList, 0);
//#endif
    png_set_user_limits(pngStruct, 5120, 5120);// 5120x5120より大きい画像は却下する。
    //    png_set_chunk_cache_max(pngStruct, 0x7FFFFFFF);// 補助チャンクsPLT, tEXt, iTXt, zTXtの限界数を設定する。
    
    png_read_info(pngStruct, pngInfo);
    
    return true;
}

bool getPNGInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize, ColorFormat* color) {
    FILE* fp;
    PNGStructures pngStructures;
    
    preloadPNG(filePath, &fp, &pngStructures);
    
    int bitDepth, colorType, interlaceMethod, compressionMethod, filterMethod;
    png_get_IHDR(pngStructures.pngStruct, pngStructures.pngInfo, width, height, &bitDepth, &colorType, &interlaceMethod, &compressionMethod, &filterMethod);
    assert(colorType != PNG_COLOR_TYPE_GA);
    if (colorType == PNG_COLOR_TYPE_RGB)
        *color = ColorFormat::RGB_8x4;
    else if (colorType == PNG_COLOR_TYPE_RGBA)
        *color = ColorFormat::RGBA8x4;
    else if (colorType == PNG_COLOR_TYPE_GRAY)
        *color = ColorFormat::Gray8;
    else if (colorType == PNG_COLOR_TYPE_PALETTE)
        *color = ColorFormat::RGB_8x4;
    
    uint32_t stride;
    if (*color == ColorFormat::RGB_8x4 || *color == ColorFormat::RGBA8x4)
        stride = 4;
    else
        stride = 1;
    size_t rowSize = *width * stride * sizeof(uint8_t);
    *requiredSize = *height * rowSize;
    
    fclose(fp);
    
    return true;
}

bool loadPNG(const std::string &filePath, uint8_t* storage, bool gammaCorrection) {
    FILE* fp;
    PNGStructures pngStructures;
    
    preloadPNG(filePath, &fp, &pngStructures);
    
    uint32_t width, height;
    
    int bitDepth, colorType, interlaceMethod, compressionMethod, filterMethod;
    png_get_IHDR(pngStructures.pngStruct, pngStructures.pngInfo, &width, &height, &bitDepth, &colorType, &interlaceMethod, &compressionMethod, &filterMethod);
    assert(colorType != PNG_COLOR_TYPE_GA);
    size_t channels = png_get_channels(pngStructures.pngStruct, pngStructures.pngInfo);
    (void)channels;
    size_t rowBytes = png_get_rowbytes(pngStructures.pngStruct, pngStructures.pngInfo);
    (void)rowBytes;
    
    if (bitDepth == 16)
        png_set_strip_16(pngStructures.pngStruct);// 16ビットの深度を8ビットに変換する。
//    png_set_invert_alpha(pngStruct);// アルファを透明度とする。
    if (bitDepth < 8)
        png_set_packing(pngStructures.pngStruct);// 1, 2, 4ビットを1バイトに詰めずに読み込みたい場合。関数名とは逆の印象。
    png_color_8p significantBit;
    if (png_get_sBIT(pngStructures.pngStruct, pngStructures.pngInfo, &significantBit))
        png_set_shift(pngStructures.pngStruct, significantBit);// 本来のビット深度が2のべき乗ではない場合に、本来の深度へと戻す。
    if (colorType == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(pngStructures.pngStruct);
//    if (colorType == PNG_COLOR_TYPE_RGB || colorType == PNG_COLOR_TYPE_RGB_ALPHA)
//        png_set_bgr(pngStruct);// RGBをBGRに変換する。
    if (colorType == PNG_COLOR_TYPE_RGB || colorType == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(pngStructures.pngStruct, 0xFFFF, PNG_FILLER_AFTER);// 3チャンネルの場合に1バイト0xFFを後ろ側に埋める。カラータイプは変わらない。
//    if (colorType == PNG_COLOR_TYPE_RGB || colorType == PNG_COLOR_TYPE_GRAY)
//        png_set_add_alpha(pngStruct, 0xFFFF, PNG_FILLER_AFTER);// アルファチャンネルを足す。
//    if (colorType == PNG_COLOR_TYPE_RGB_ALPHA)
//        png_set_swap_alpha(pngStruct);// RGBAのファイルをARGBで読み込みたい場合。
//    if (colorType == PNG_COLOR_TYPE_GRAY || colorType == PNG_COLOR_TYPE_GRAY_ALPHA)
//        png_set_gray_to_rgb(pngStruct);// グレイスケール画像をRGBで表現させて読み込みたい場合。
//    if (colorType == PNG_COLOR_TYPE_RGB || colorType == PNG_COLOR_TYPE_RGB_ALPHA)
//        png_set_rgb_to_gray_fixed(pngStruct, 1, 21268, 71510);RGB画像をグレイスケールで読み込みたい場合。RとGのウェイトx100000を指定する。
//    png_color_16 myBackground;
//    png_color_16p imageBackground;
//    if (png_get_bKGD(pngStruct, pngInfo, &imageBackground))
//        png_set_background(pngStruct, imageBackground, PNG_BACKGROUND_GAMMA_FILE, 1, 1.0);
//    else
//        png_set_background(pngStruct, &myBackground, PNG_BACKGROUND_GAMMA_SCREEN, 0, 1.0);
    double gamma, screenGamma;
    const char* gammaStr;
    if (gammaCorrection) {
        if ((gammaStr = getenv("SCREEN_GAMMA")) != nullptr)
            screenGamma = (double)atof(gammaStr);
        else
            screenGamma = 2.2;
    }
    else{
        screenGamma = 1.0f;
    }
    if (png_get_gAMA(pngStructures.pngStruct, pngStructures.pngInfo, &gamma))
        png_set_gamma(pngStructures.pngStruct, screenGamma, gamma);// ファイルにガンマ値がある場合。
    else
        png_set_gamma(pngStructures.pngStruct, screenGamma, 0.45455);// ファイルにガンマ値が無い場合。
//    if (bitDepth == 1 && colorType == PNG_COLOR_TYPE_GRAY)
//        png_set_invert_mono(pngStruct);// 2値画像の解釈を反転する。
//    if (colorType == PNG_COLOR_TYPE_GRAY || colorType == PNG_COLOR_TYPE_GRAY_ALPHA)
//        png_set_invert_mono(pngStruct);// グレイスケール全般の解釈を反転する。
//    if (bitDepth == 16)
//        png_set_swap(pngStruct);// 16ビット画像のエンディアンを変更する。
//    if (bitDepth < 8)
//        png_set_packswap(pngStruct);// 4ビット以下の深度の場合のパッキングの順序を変更する。
//    png_set_read_user_transform_fn(pngStruct, );// 自作変換関数を登録する。
//    int numPasses = png_set_interlace_handling(pngStruct);インターレース処理の登録。
    png_read_update_info(pngStructures.pngStruct, pngStructures.pngInfo);
    
    uint32_t stride;
    if (colorType == PNG_COLOR_TYPE_RGB || colorType == PNG_COLOR_TYPE_RGBA || colorType == PNG_COLOR_TYPE_PALETTE)
        stride = 4;
    else
        stride = 1;
    size_t rowSize = width * stride * sizeof(uint8_t);
    uint8_t* dataHead = storage;
    for (int i = 0; i < height; ++i) {
        png_read_row(pngStructures.pngStruct, dataHead, nullptr);
        dataHead += rowSize;
    }
    
    png_read_end(pngStructures.pngStruct, pngStructures.pngEndInfo);
    png_destroy_read_struct(&pngStructures.pngStruct, &pngStructures.pngInfo, &pngStructures.pngEndInfo);
    
    fclose(fp);
    
    return true;
}

bool getImageInfo(const std::string &filePath, uint32_t* width, uint32_t* height, uint64_t* requiredSize, ColorFormat* color) {
    size_t extPos = filePath.find_last_of(".");
    if (extPos == std::string::npos)
        return false;
    
    std::string ext = filePath.substr(extPos + 1);
    if (!ext.compare("jpg") || !ext.compare("jpeg")) {
        *color = ColorFormat::RGB8x3;
        return getJPEGInfo(filePath, width, height, requiredSize);
    }
    else if (!ext.compare("png")) {
        return getPNGInfo(filePath, width, height, requiredSize, color);
    }
    else if (!ext.compare("exr")) {
        EXRType exrType;
        *color = ColorFormat::RGBA16Fx4;
        return getEXRInfo(filePath, width, height, requiredSize, &exrType);
    }
    
    return false;
}

bool loadImage(const std::string &filePath, uint8_t* storage, bool gammaCorrection) {
    size_t extPos = filePath.find_last_of(".");
    if (extPos == std::string::npos)
        return false;
    
    std::string ext = filePath.substr(extPos + 1);
    if (!ext.compare("jpg") || !ext.compare("jpeg")) {
        return loadJPEG(filePath, storage);
    }
    else if (!ext.compare("png")) {
        return loadPNG(filePath, storage, gammaCorrection);
    }
    else if (!ext.compare("exr")) {
        return loadEXR(filePath, storage);
    }
    
    return false;
}
