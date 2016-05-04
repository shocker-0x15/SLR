//
//  surface_material.cpp
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_material.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/directional_distribution_functions.h"
#include "../Core/textures.h"
#include "../Core/XORShift.h"

namespace SLR {
    SubSurfaceScatteringSurfaceMaterial::SubSurfaceScatteringSurfaceMaterial(const SurfaceMaterial* surfMat,
                                                                             const InputSpectrum* l_sigma_a, const InputSpectrum* l_sigma_s, float l_g,
                                                                             const InputSpectrum* u_sigma_a, const InputSpectrum* u_sigma_s, float u_g) :
    m_surfMat(surfMat),
    m_l_sigma_a(l_sigma_a), m_l_sigma_s(l_sigma_s), m_l_g(l_g),
    m_u_sigma_a(u_sigma_a), m_u_sigma_s(u_sigma_s), m_u_g(u_g) {
        ArenaAllocator tempMem;
        
#ifdef Use_Spectral_Representation
        const uint32_t numPrecomputeWavelengths = 64;
        float lambdas[numPrecomputeWavelengths];
        for (int i = 0; i < numPrecomputeWavelengths; ++i)
            lambdas[i] = WavelengthLowBound + (WavelengthHighBound - WavelengthLowBound) / (numPrecomputeWavelengths - 1) * i;
        WavelengthSamples wls = WavelengthSamples(lambdas);
#else
        WavelengthSamples wls;
#endif
        
        SurfacePoint surfPt;
        surfPt.texCoord = TexCoord2D(0, 0);
        BSDF* innerBSDF = surfMat->getBSDF(surfPt, wls, tempMem);
        
        const uint32_t numSamples = 128;
        BSDFSample* bsdfSamples = tempMem.alloc<BSDFSample>(numSamples);
        float* uDir0 = tempMem.alloc<float>(numSamples);
        float* uDir1 = tempMem.alloc<float>(numSamples);
        float* uWl = tempMem.alloc<float>(numSamples);
        XORShift rng{394891282};
        for (int i = 0; i < numSamples; ++i) {
            new (&bsdfSamples[i]) BSDFSample(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), rng.getFloat0cTo1o());
            uDir0[i] = rng.getFloat0cTo1o();
            uDir1[i] = rng.getFloat0cTo1o();
            uWl[i] = rng.getFloat0cTo1o();
        }
        
        SampledSpectrum innerHHReflectance = innerBSDF->rho(numSamples, bsdfSamples, uDir0, uDir1, uWl);
        
#ifdef Use_Spectral_Representation
//        RegularContinuousSpectrumTemplate(RealType minWL, RealType maxWL, const RealType* vals, uint32_t numVals)
        m_innerHHReflectance = new RegularContinuousSpectrum(WavelengthLowBound, WavelengthHighBound, innerHHReflectance.values, numPrecomputeWavelengths);
#else
        m_innerHHReflectance = new RGBInputSpectrum(innerHHReflectance);
#endif
    }
    
    BSSRDF* SubSurfaceScatteringSurfaceMaterial::getBSSRDF(bool lowerHemisphere, const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        SampledSpectrum internalHHReflectance;
        if (lowerHemisphere && m_l_sigma_a != nullptr)
            return mem.create<BSSRDF>(m_l_sigma_a->evaluate(wls), m_l_sigma_s->evaluate(wls), m_l_g, internalHHReflectance);
        else if (m_u_sigma_a != nullptr)
            return mem.create<BSSRDF>(m_u_sigma_a->evaluate(wls), m_u_sigma_s->evaluate(wls), m_u_g, internalHHReflectance);
        return nullptr;
    }
    
    
    
    Fresnel* SVFresnelNoOp::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelNoOp>();
    }
    
    Fresnel* SVFresnelConductor::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelConductor>(m_eta->evaluate(surfPt.texCoord, wls), m_k->evaluate(surfPt.texCoord, wls));
    }
    
    Fresnel* SVFresnelDielectric::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelDielectric>(m_etaExt->evaluate(surfPt.texCoord, wls), m_etaInt->evaluate(surfPt.texCoord, wls));
    }
    
    
    
    MicrofacetDistribution* SVGGX::getMicrofacetDistribution(const SLR::SurfacePoint &surfPt, SLR::ArenaAllocator &mem) const {
        return mem.create<GGX>(m_alpha_g->evaluate(surfPt.texCoord));
    }
}
