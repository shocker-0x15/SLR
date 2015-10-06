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

#include <libSLR/BasicTypes/Spectrum.h>

namespace SLRSceneGraph {
    namespace API {
        enum class Type : uint32_t {
            Integer,
            Float,
            String,
            Array,
            Image,
            Node
        };
        
        Image2DRef LoadImage(const std::string &filepath);
        
        NodeRef Set3DModel(const std::string &filepath);
        
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
        
        namespace SurfaceMaterial {
            SurfaceMaterialRef createMatte(const SpectrumTextureRef &reflectance, const FloatTextureRef &sigma);
            SurfaceMaterialRef createMetal(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &eta, const SpectrumTextureRef &k);
            SurfaceMaterialRef createGlass(const SpectrumTextureRef &coeffR, const SpectrumTextureRef &coeffT, const SpectrumTextureRef &etaExt, const SpectrumTextureRef &etaInt);
            SurfaceMaterialRef createModifiedWardDur(const SpectrumTextureRef &reflectance, const FloatTextureRef &anisoX, const FloatTextureRef &anisoY);
            SurfaceMaterialRef createAshikhminShirley(const SpectrumTextureRef &Rd, const SpectrumTextureRef &Rs, const FloatTextureRef &nx, const FloatTextureRef &ny);
            SurfaceMaterialRef createSummedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1);
            SurfaceMaterialRef createMixedMaterial(const SurfaceMaterialRef &mat0, const SurfaceMaterialRef &mat1, const FloatTextureRef &factor);
            EmitterSurfacePropertyRef createDiffuseEmitter(const SpectrumTextureRef &emittance);
            SurfaceMaterialRef createEmitterSurfaceMaterial(const SurfaceMaterialRef &mat, const EmitterSurfacePropertyRef &emit);
        }
    }
}

#endif /* API_hpp */
