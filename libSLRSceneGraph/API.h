//
//  API.h
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_API__
#define __SLRSceneGraph_API__

#include <libSLR/defines.h>
#include <libSLR/BasicTypes/spectrum_base.h>
#include <libSLR/Core/image_2d.h>
#include "declarations.h"
#include "Scene/Scene.h"
#include "Parser/SceneParser.h"

namespace SLRSceneGraph {
    SLR_SCENEGRAPH_API bool readScene(const std::string &filePath, const SceneRef &scene, RenderingContext* context);
    
    namespace Spectrum {
        using namespace SLR;
        
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2);
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples);
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples);
    }
}

#endif /* API_hpp */
