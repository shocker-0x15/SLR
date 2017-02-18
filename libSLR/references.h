//
//  references.h
//
//  Created by 渡部 心 on 2015/04/27.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_references_h
#define SLR_references_h

#include <vector>
#include <memory>
#include <stdint.h>

namespace SLR {
    // Memory Allocators
    class Allocator;
    class DefaultAllocator;
    class StackAllocator;
    class ArenaAllocator;
    class MSpaceAllocator;
    
    // RGB & Spectrum
    template <typename RealType> struct CompensatedSum;
    template <typename RealType> struct RGBSamplesTemplate;
    template <typename RealType> struct RGBTemplate;
    template <typename RealType> class RGBStorageTemplate;
    template <typename RealType, uint32_t N> class ContinuousSpectrumTemplate;
    template <typename RealType, uint32_t N> class RegularContinuousSpectrumTemplate;
    template <typename RealType, uint32_t N> class IrregularContinuousSpectrumTemplate;
    template <typename RealType, uint32_t N> class UpsampledContinuousSpectrumTemplate;
    template <typename RealType, uint32_t N> class ScaledAndOffsetUpsampledContinuousSpectrumTemplate;
    template <typename RealType, uint32_t N> struct WavelengthSamplesTemplate;
    template <typename RealType, uint32_t N> struct SampledSpectrumTemplate;
    template <typename RealType, uint32_t N> struct DiscretizedSpectrumTemplate;
    template <typename RealType, uint32_t numStrata> class SpectrumStorageTemplate;
    // FIXME: Current code is inconsistent with respect to float precision.
    typedef float SpectrumFloat;
    typedef RGBTemplate<SpectrumFloat> RGBInputSpectrum;
    const static uint32_t NumSpectralSamples = 16;
    const static uint32_t NumStrataForStorage = 16;
    typedef ContinuousSpectrumTemplate<SpectrumFloat, NumSpectralSamples> ContinuousSpectrum;
    typedef RegularContinuousSpectrumTemplate<SpectrumFloat, NumSpectralSamples> RegularContinuousSpectrum;
    typedef IrregularContinuousSpectrumTemplate<SpectrumFloat, NumSpectralSamples> IrregularContinuousSpectrum;
    typedef UpsampledContinuousSpectrumTemplate<SpectrumFloat, NumSpectralSamples> UpsampledContinuousSpectrum;
    typedef ScaledAndOffsetUpsampledContinuousSpectrumTemplate<SpectrumFloat, NumSpectralSamples> ScaledAndOffsetUpsampledContinuousSpectrum;
#ifdef Use_Spectral_Representation
    typedef WavelengthSamplesTemplate<SpectrumFloat, NumSpectralSamples> WavelengthSamples;
    typedef SampledSpectrumTemplate<SpectrumFloat, NumSpectralSamples> SampledSpectrum;
    typedef DiscretizedSpectrumTemplate<SpectrumFloat, NumStrataForStorage> DiscretizedSpectrum;
    typedef SpectrumStorageTemplate<SpectrumFloat, NumStrataForStorage> SpectrumStorage;
    
    typedef ContinuousSpectrum InputSpectrum;
#else
    typedef RGBSamplesTemplate<SpectrumFloat> WavelengthSamples;
    typedef RGBTemplate<SpectrumFloat> SampledSpectrum;
    typedef RGBTemplate<SpectrumFloat> DiscretizedSpectrum;
    typedef RGBStorageTemplate<SpectrumFloat> SpectrumStorage;
    
    typedef RGBInputSpectrum InputSpectrum;
#endif
    typedef CompensatedSum<SampledSpectrum> SampledSpectrumSum;
    
    // Basic Types
    template <typename RealType> struct Vector3Template;
    template <typename RealType> struct Vector4Template;
    template <typename RealType> struct Normal3Template;
    template <typename RealType> struct Point3Template;
    template <typename RealType> struct Matrix3x3Template;
    template <typename RealType> struct Matrix4x4Template;
    template <typename RealType> struct QuaternionTemplate;
    template <typename RealType> struct TexCoord2Template;
    template <typename RealType> struct RayTemplate;
    template <typename RealType> struct BoundingBox3Template;
    typedef Vector3Template<float> Vector3D;
    typedef Vector3Template<float> Tangent3D;
    typedef Vector3Template<float> Bitangent3D;
    typedef Vector4Template<float> Vector4D;
    typedef Normal3Template<float> Normal3D;
    typedef Point3Template<float> Point3D;
    typedef Matrix3x3Template<float> Matrix3x3;
    typedef Matrix4x4Template<float> Matrix4x4;
    typedef QuaternionTemplate<float> Quaternion;
    typedef TexCoord2Template<float> TexCoord2D;
    typedef RayTemplate<float> Ray;
    typedef BoundingBox3Template<float> BoundingBox3D;
    typedef CompensatedSum<float> FloatSum;
    
    // Transforms
    class Transform;
    class StaticTransform;
    class AnimatedTransform;
    class ChainedTransform;
    
    // Random Number Generators
    struct Types32bit;
    struct Types64bit;
    template <typename TypeSet> class RandomNumberGeneratorTemplate;
    template <typename TypeSet> class XORShiftRNGTemplate;
    template <typename TypeSet> class LinearCongruentialRNGTemplate;
    
    typedef RandomNumberGeneratorTemplate<Types32bit> RandomNumberGenerator;
    typedef XORShiftRNGTemplate<Types32bit> XORShiftRNG;
    typedef LinearCongruentialRNGTemplate<Types32bit> LinearCongruentialRNG;
    
    // Distributions
    template <typename RealType> class RegularConstantDiscrete1DTemplate;
    template <typename RealType> class RegularConstantContinuous1DTemplate;
    template <typename RealType> class RegularConstantContinuous2DTemplate;
    typedef RegularConstantDiscrete1DTemplate<float> RegularConstantDiscrete1D;
    typedef RegularConstantContinuous1DTemplate<float> RegularConstantContinuous1D;
    typedef RegularConstantContinuous2DTemplate<float> RegularConstantContinuous2D;
    
    // Image & Tiled Image
    class Image2D;
    template <uint32_t log2_tileWidth = 3> class TiledImage2DTemplate;
    
    typedef TiledImage2DTemplate<> TiledImage2D;
    
    //
    class ImageSensor;
    
    // geometry
    class Interaction;
    class SurfaceInteraction;
    class MediumInteraction;
    struct ReferenceFrame;
    class InteractionPoint;
    class SurfacePoint;
    class MediumPoint;
    
    // Surfaces
    class Surface;
    struct Vertex;
    class Triangle;
    class InfiniteSphere;
    class InfinitesimalPoint;
    
    // Media
    class Medium;
    class HomogeneousMedium;
    class AchromaticExtinctionGridMedium;
    class GridMedium;
    class SubGridMedium;
    class DensityGridMedium;
    
    class Object;
    
    struct LightPosQuery;
    struct LightPosQueryResult;
    struct SurfaceLightPosSample;
    struct SurfaceLightPosQueryResult;
    struct VolumetricLightPosSample;
    struct VolumetricLightPosQueryResult;
    class Light;
    class SurfaceLight;
    class VolumetricLight;
    
    // Surface Objects
    class SurfaceObject;
    class SingleSurfaceObject;
    class BumpSingleSurfaceObject;
    class InfiniteSphereSurfaceObject;
    class TransformedSurfaceObject;
    class SurfaceObjectAggregate;
    
    // Medium Objects
    class MediumObject;
    class SingleMediumObject;
    class TransformedMediumObject;
    class EnclosedMediumObject;
    class MediumObjectAggregate;
    
    // Cameras
    class Camera;
    class PerspectiveCamera;
    class EquirectangularCamera;
    
    // Path Samplers
    class FreePathSampler;
    class LightPathSampler;
    class IndependentLightPathSampler;
    
    // Accelerators
    class Accelerator;
    class StandardBVH;
    class SBVH;
    class QBVH;
    
    // Directional Distribution Functions
    struct DirectionType;
    struct EDFQuery;
    struct EDFSample;
    struct EDFQueryResult;
    struct ABDFQuery;
    struct ABDFSample;
    struct ABDFReverseInfo;
    struct ABDFQueryResult;
    struct BSDFQuery;
    struct BSDFSample;
    struct BSDFReverseInfo;
    struct BSDFQueryResult;
    struct VolumetricBSDFQuery;
    struct VolumetricBSDFSample;
    struct VolumetricBSDFReverseInfo;
    struct VolumetricBSDFQueryResult;
    struct IDFQuery;
    struct IDFSample;
    struct IDFQueryResult;
    struct PFSample;
    struct PFQueryResult;
    class EDF;
    class AbstractBDF;
    class BSDF;
    class PhaseFunction;
    class VolumetricBSDF;
    class IDF;
    class DiffuseEDF;
    class IdealDirectionalEDF;
    class IBLEDF;
    class MultiEDF;
    class Fresnel;
    class FresnelNoOp;
    class FresnelConductor;
    class FresnelDielectric;
    class MicrofacetDistribution;
    class GGX;
    class LambertianBRDF;
    class SpecularBRDF;
    class SpecularBSDF;
    class InverseBSDF;
    class NullBSDF;
    class OrenNayerBRDF;
    class ModifiedWardDurBRDF;
    class AshikhminShirleyBRDF;
    class MicrofacetBRDF;
    class MicrofacetBSDF;
    class MultiBSDF;
    class IsotropicPhaseFunction;
    class HenyeyGreensteinPhaseFunction;
    class SchlickPhaseFunction;
    class PerspectiveIDF;
    class EquirectangularIDF;
    
    // Textures & Mappings
    class Texture2DMapping;
    class Texture3DMapping;
    class OffsetAndScale2DMapping;
    class OffsetAndScale3DMapping;
    class WorldPosition3DMapping;
    class SpectrumTexture;
    class Normal3DTexture;
    class FloatTexture;
    class ConstantSpectrumTexture;
    class ConstantFloatTexture;
    class ImageSpectrumTexture;
    class ImageNormal3DTexture;
    class ImageFloatTexture;
    class CheckerBoardSpectrumTexture;
    class CheckerBoardNormal3DTexture;
    class CheckerBoardFloatTexture;
    class VoronoiSpectrumTexture;
    class VoronoiNormal3DTexture;
    class VoronoiFloatTexture;
    
    // Surface Materials
    class SVFresnel;
    class SVFresnelNoOp;
    class SVFresnelConductor;
    class SVFresnelDielectric;
    class SVMicrofacetDistribution;
    class SVGGX;
    class SurfaceMaterial;
    class DiffuseReflection;
    class SpecularReflection;
    class SpecularScattering;
    class ModifiedWardDurReflection;
    class AshikhminShirleyReflection;
    class MicrofacetReflection;
    class MicrofacetScattering;
    class SummedSurfaceMaterial;
    class MixedSurfaceMaterial;
    class EmitterSurfaceProperty;
    class DiffuseEmission;
    class IdealDirectionalEmission;
    class IBLEmission;
    class EmitterSurfaceMaterial;
    
    // Medium Materials
    class MediumMaterial;
    class IsotropicScattering;
    class EmitterMediumProperty;
    class EmitterMediumMaterial;
    
    // Scene
    class Node;
    class InternalNode;
    class ReferenceNode;
    class InfiniteSphereNode;
    class InfinitesimalPointNode;
    class PerspectiveCameraNode;
    class EquirectangularCameraNode;
    class SurfaceNode;
    class TriangleMeshNode;
    class MediumNode;
    class HomogeneousMediumNode;
    class GridMediumNode;
    class DensityGridMediumNode;
    class Scene;
    
    class RenderSettings;
    class Renderer;
    
    // Renderers
    class PathTracingRenderer;
    class BidirectionalPathTracingRenderer;
    class AMCMCPPMRenderer;
    class VolumetricPathTracingRenderer;
    class VolumetricBPTRenderer;
    
    class ProgressReporter;
}

#endif
