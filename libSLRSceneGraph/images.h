//
//  images.h
//
//  Created by 渡部 心 on 2017/05/16.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_images__
#define __SLRSceneGraph_images__

#include <libSLR/defines.h>
#include "declarations.h"

#include <libSLR/Core/image_2d.h>

namespace SLRSceneGraph {
    struct Image2DLoadParams {
        std::string filePath;
        SLR::ImageStoreMode storeMode;
        SLR::SpectrumType spectrumType;
        bool gammaCorrection;
        
        bool operator<(const Image2DLoadParams &params) const {
            if (filePath < params.filePath) {
                return true;
            }
            else if (filePath == params.filePath) {
                if (storeMode < params.storeMode) {
                    return true;
                }
                else if (storeMode == params.storeMode) {
                    if (spectrumType < params.spectrumType) {
                        return true;
                    }
                    else if (spectrumType == params.spectrumType) {
                        if (gammaCorrection < params.gammaCorrection)
                            return true;
                    }
                }
            }
            return false;
        }
    };
    
    extern std::map<Image2DLoadParams, Image2DRef> s_imageDB;
    
    SLR_SCENEGRAPH_API Image2DRef createImage2D(const std::string &filepath, SLR::ImageStoreMode mode, SLR::SpectrumType spType, bool gammaCorrection);
    
    
    
    class SLR_SCENEGRAPH_API Image2D {
    protected:
        SLR::Image2D* m_rawData;
    public:
        virtual ~Image2D();
        const SLR::Image2D* getRaw() const {
            return m_rawData;
        };
    };
    
    
    
    class SLR_SCENEGRAPH_API TiledImage2D : public Image2D {
        std::string m_filePath;
        SLR::ImageStoreMode m_storeMode;
        SLR::SpectrumType m_spectrumType;
        bool m_gammaCorrection;
    public:
        TiledImage2D(const std::string &filePath, SLR::ImageStoreMode storeMode, SLR::SpectrumType spectrumType, bool gammaCorrection);
    };
}

#endif /* __SLRSceneGraph_images__ */
