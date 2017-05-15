//
//  images.cpp
//
//  Created by 渡部 心 on 2017/05/16.
//  Copyright © 2017年 渡部 心. All rights reserved.
//

#include "images.h"

#include <libSLR/Core/image_2d.h>

#include "Helper/image_loader.h"

namespace SLRSceneGraph {
    std::map<Image2DLoadParams, Image2DRef> s_imageDB;
    
    SLR_SCENEGRAPH_API Image2DRef createImage2D(const std::string &filepath, SLR::ImageStoreMode mode, SLR::SpectrumType spType, bool gammaCorrection) {
        Image2DLoadParams params{filepath, mode, spType, gammaCorrection};
        if (s_imageDB.count(params) > 0) {
            return s_imageDB[params];
        }
        else {
            Image2DRef ret = createShared<TiledImage2D>(filepath, mode, spType, gammaCorrection);
            s_imageDB[params] = ret;
            return ret;
        }
    };
    
    
    
    Image2D::~Image2D() {
        delete m_rawData;
    }
    
    
    
    TiledImage2D::TiledImage2D(const std::string &filePath, SLR::ImageStoreMode storeMode, SLR::SpectrumType spectrumType, bool gammaCorrection) : 
    m_filePath(filePath), m_storeMode(storeMode), m_spectrumType(spectrumType), m_gammaCorrection(gammaCorrection) {
        uint64_t requiredSize;
        bool imgSuccess;
        uint32_t width, height;
        ::ColorFormat colorFormat;
        imgSuccess = getImageInfo(m_filePath, &width, &height, &requiredSize, &colorFormat);
        SLRAssert(imgSuccess, "Error occured during getting image information.\n%s", m_filePath.c_str());
        
        void* linearData = malloc(requiredSize);
        imgSuccess = loadImage(m_filePath, (uint8_t*)linearData, m_gammaCorrection);
        SLRAssert(imgSuccess, "failed to load the image\n%s", m_filePath.c_str());
        
        SLR::ColorFormat internalFormat = (SLR::ColorFormat)colorFormat;
        
        // TODO: ?? make a memory allocator selectable.
        SLR::DefaultAllocator &defMem = SLR::DefaultAllocator::instance();
        m_rawData = new SLR::TiledImage2D(linearData, width, height, internalFormat, &defMem, storeMode, spectrumType);
        free(linearData);
    }
}
