//
//  BidirectionalPathTracingRenderer.h
//
//  Created by 渡部 心 on 2016/01/31.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef BidirectionalPathTracingRenderer_h
#define BidirectionalPathTracingRenderer_h

#include "../defines.h"
#include "../references.h"
#include "../Core/Renderer.h"

#include "../Core/geometry.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API BidirectionalPathTracingRenderer : public Renderer {
        struct SLR_API DDFQuery {
            Vector3D dir_sn;
            Normal3D gNormal_sn;
            uint8_t heroIndex;
            bool adjoint;
        };
        
        struct DDFProxy {
            virtual const void* getDDF() const = 0;
            virtual SampledSpectrum evaluate(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const = 0;
            virtual SampledSpectrum evaluatePDF(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal = nullptr) const = 0;
        };
        struct EDFProxy : public DDFProxy {
            const EDF* edf;
            EDFProxy(const EDF* _edf) : edf(_edf) {}
            const void* getDDF() const override { return edf; }
            SampledSpectrum evaluate(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                EDFQuery edfQuery(query.heroIndex);
                return edf->evaluate(edfQuery, dir_sn);
            }
            SampledSpectrum evaluatePDF(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                EDFQuery edfQuery(query.heroIndex);
                return edf->evaluatePDF(edfQuery, dir_sn);
            }
        };
        struct BSDFProxy : public DDFProxy {
            const BSDF* bsdf;
            BSDFProxy(const BSDF* _bsdf) : bsdf(_bsdf) {}
            const void* getDDF() const override { return bsdf; }
            SampledSpectrum evaluate(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                BSDFQuery bsdfQuery(query.heroIndex, query.dir_sn, query.gNormal_sn, DirectionType::All, query.adjoint);
                return bsdf->evaluate(bsdfQuery, dir_sn, revVal);
            }
            SampledSpectrum evaluatePDF(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                BSDFQuery bsdfQuery(query.heroIndex, query.dir_sn, query.gNormal_sn, DirectionType::All, query.adjoint);
                return bsdf->evaluatePDF(bsdfQuery, dir_sn, revVal);
            }
        };
        struct IDFProxy : public DDFProxy {
            const IDF* idf;
            IDFProxy(const IDF* _idf) : idf(_idf) {}
            const void* getDDF() const override { return idf; }
            SampledSpectrum evaluate(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return idf->evaluate(dir_sn);
            }
            SampledSpectrum evaluatePDF(const DDFQuery &query, const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return idf->evaluatePDF(dir_sn);
            }
        };
        
        struct BPTVertex {
            SurfacePoint surfPt;
            Vector3D dirIn_sn;
            Normal3D gNormal_sn;
            const DDFProxy* ddf;
            SampledSpectrum wlPDFs;
            SampledSpectrum alpha;
            SampledSpectrum areaPDF;
            SampledSpectrum RRProb;
            SampledSpectrum revAreaPDF;
            SampledSpectrum revRRProb;
            DirectionType sampledType;
            BPTVertex(const SurfacePoint &_surfPt, const Vector3D &_dirIn_sn, const Normal3D &_gNormal_sn, const DDFProxy* _ddf,
                      const SampledSpectrum &_wlPDFs, const SampledSpectrum &_alpha, const SampledSpectrum &_areaPDF, const SampledSpectrum &_RRProb, DirectionType _sampledType) :
            surfPt(_surfPt), dirIn_sn(_dirIn_sn), gNormal_sn(_gNormal_sn), ddf(_ddf),
            wlPDFs(_wlPDFs), alpha(_alpha), areaPDF(_areaPDF), RRProb(_RRProb), revAreaPDF(NAN), revRRProb(NAN), sampledType(_sampledType) {}
        };
        
        struct Job {
            const Scene* scene;
            
            ArenaAllocator* mems;
            IndependentLightPathSampler** pathSamplers;
            
            const Camera* camera;
            float timeStart;
            float timeEnd;
            
            ImageSensor* sensor;
            uint32_t imageWidth;
            uint32_t imageHeight;
            uint32_t numPixelX;
            uint32_t numPixelY;
            uint32_t basePixelX;
            uint32_t basePixelY;
            
            // working area
            uint8_t heroIndex;
            SampledSpectrum wlPDFs;
            float curPx, curPy;
            std::vector<BPTVertex> lightVertices;
            std::vector<BPTVertex> eyeVertices;
            
            void kernel(uint32_t threadID);
            void generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initWLPDFs,
                                 const SampledSpectrum &initAlpha, const SLR::Ray &initRay, const SampledSpectrum &initDirPDF, DirectionType sampledType,
                                 float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, SLR::ArenaAllocator &mem);
            float calculateMISWeight(const SampledSpectrum &lExtend1stAreaPDF, const SampledSpectrum &lExtend1stRRProb,
                                     const SampledSpectrum &lExtend2ndAreaPDF, const SampledSpectrum &lExtend2ndRRProb,
                                     const SampledSpectrum &eExtend1stAreaPDF, const SampledSpectrum &eExtend1stRRProb,
                                     const SampledSpectrum &eExtend2ndAreaPDF, const SampledSpectrum &eExtend2ndRRProb,
                                     uint32_t numLVtx, uint32_t numEVtx) const;
        };
        
        uint32_t m_samplesPerPixel;
    public:
        BidirectionalPathTracingRenderer(uint32_t spp);
        void render(const Scene &scene, const RenderSettings &settings) const override;
    };
}

#endif /* BidirectionalPathTracingRenderer_hpp */
