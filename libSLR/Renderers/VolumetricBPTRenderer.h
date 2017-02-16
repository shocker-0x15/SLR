//
//  VolumetricBPTRenderer.h
//
//  Created by 渡部 心 on 2017/02/16.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_VolumetricBPTRenderer__
#define __SLR_VolumetricBPTRenderer__

#include "../defines.h"
#include "../references.h"
#include "../Core/Renderer.h"

#include "../Core/geometry.h"
#include "../Core/directional_distribution_functions.h"

namespace SLR {
    class SLR_API VolumetricBPTRenderer : public Renderer {
        struct SLR_API DDFQuery {
            Vector3D dir_sn;
            Normal3D gNormal_sn;
            int16_t wlHint;
            bool adjoint;
        };
        
        struct DDFProxy {
            virtual const void* getDDF() const = 0;
            virtual SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const = 0;
            virtual float evaluatePDF(const Vector3D &dir_sn, float* revVal = nullptr) const = 0;
        };
        struct EDFProxy : public DDFProxy {
            const EDF* edf;
            EDFQuery query;
            
            EDFProxy(const EDF* _edf, const EDFQuery &_query) : edf(_edf), query(_query) {}
            const void* getDDF() const override { return edf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return edf->evaluate(query, dir_sn);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return edf->evaluatePDF(query, dir_sn);
            }
        };
        struct ABDFProxy : public DDFProxy {
            const AbstractBDF* abdf;
            const ABDFQuery* query;
            
            ABDFProxy(const AbstractBDF* _abdf, const ABDFQuery* _query) : abdf(_abdf), query(_query) {}
            const void* getDDF() const override { return abdf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return abdf->evaluate(query, dir_sn, revVal);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return abdf->evaluatePDF(query, dir_sn, revVal);
            }
        };
        struct IDFProxy : public DDFProxy {
            const IDF* idf;
            
            IDFProxy(const IDF* _idf) : idf(_idf) {}
            const void* getDDF() const override { return idf; }
            SampledSpectrum evaluate(const Vector3D &dir_sn, SampledSpectrum* revVal) const override {
                return idf->evaluate(dir_sn);
            }
            float evaluatePDF(const Vector3D &dir_sn, float* revVal) const override {
                return idf->evaluatePDF(dir_sn);
            }
        };
        
        struct VBPTVertex {
            const InteractionPoint* interPt;
            const DDFProxy* ddf;
            SampledSpectrum alpha;
            float cosIn;
            float spatialPDF;
            float RRProb;
            float revSpatialPDF;
            float revRRProb;
            DirectionType sampledType;
            bool lambdaSelected;
            VBPTVertex(const InteractionPoint* _interPt, const DDFProxy* _ddf,
                      const SampledSpectrum &_alpha, float _cosIn, float _spatialPDF, float _RRProb, DirectionType _sampledType, bool _lambdaSelected) :
            interPt(_interPt), ddf(_ddf),
            alpha(_alpha), cosIn(_cosIn), spatialPDF(_spatialPDF), RRProb(_RRProb), revSpatialPDF(NAN), revRRProb(NAN), sampledType(_sampledType), lambdaSelected(_lambdaSelected) {}
        };
        
        struct Job {
            const Scene* scene;
            
            ArenaAllocator* mems;
            IndependentLightPathSampler* pathSamplers;
            
            const Camera* camera;
            ImageSensor* sensor;
            float timeStart;
            float timeEnd;
            uint32_t imageWidth;
            uint32_t imageHeight;
            uint32_t numPixelX;
            uint32_t numPixelY;
            uint32_t basePixelX;
            uint32_t basePixelY;
            
            // working area
            float curPx, curPy;
            int16_t wlHint;
            std::vector<VBPTVertex> lightVertices;
            std::vector<VBPTVertex> eyeVertices;
            
            ProgressReporter* reporter;
            
            void kernel(uint32_t threadID);
            void generateSubPath(const WavelengthSamples &initWLs, const SampledSpectrum &initAlpha, const Ray &initRay, float dirPDF, DirectionType sampledType,
                                 float cosLast, bool adjoint, IndependentLightPathSampler &pathSampler, ArenaAllocator &mem);
            float calculateMISWeight(float lExtend1stSpatialPDF, float lExtend1stRRProb, float lExtend2ndSpatialPDF, float lExtend2ndRRProb,
                                     float eExtend1stSpatialPDF, float eExtend1stRRProb, float eExtend2ndSpatialPDF, float eExtend2ndRRProb,
                                     uint32_t numLVtx, uint32_t numEVtx) const;
        };
        
        uint32_t m_samplesPerPixel;
    public:
        VolumetricBPTRenderer(uint32_t spp);
        void render(const Scene &scene, const RenderSettings &settings) const override;
    };
}

#endif /* __SLR_VolumetricBPTRenderer__ */
