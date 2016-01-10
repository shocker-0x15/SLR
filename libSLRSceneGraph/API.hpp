//
//  API.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_API_hpp
#define SLRSceneGraph_API_hpp

#include <libSLR/defines.h>
#include "references.h"

#include "Scene.h"
#include "Parser/SceneParser.hpp"
#include <libSLR/BasicTypes/Spectrum.h>

namespace SLRSceneGraph {
    bool readScene(const std::string &filePath, const SceneRef &scene, RenderingContext* context);
    
    Element SetRenderer(const ParameterList &params, RenderingContext* context, ErrorMessage* err);
    
    namespace Spectrum {
        using namespace SLR;
        
        InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2);
        InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples);
        InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples);
    }
    
    namespace Image {
        extern std::map<std::string, Image2DRef> s_imageDB;
        
        std::shared_ptr<SLR::TiledImage2D> createTiledImage(const std::string &filepath, SLR::Allocator *mem, SLR::SpectrumType spType, bool gammaCorrection = false);
    }
}

#endif /* API_hpp */
