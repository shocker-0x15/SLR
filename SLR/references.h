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

// Memory Allocators
class Allocator;
class DefaultAllocator;
class StackAllocator;
class ArenaAllocator;
class MSpaceAllocator;

// Basic Types
template <typename RealType> struct SpectrumTemplate;
template <typename RealType> struct Vector3Template;
template <typename RealType> struct Vector4Template;
template <typename RealType> struct Normal3Template;
template <typename RealType> struct Point3Template;
template <typename RealType> struct Matrix4x4Template;
template <typename RealType> struct QuaternionTemplate;
template <typename RealType> struct TexCoord2Template;
typedef SpectrumTemplate<float> Spectrum;
typedef Vector3Template<float> Vector3D;
typedef Vector3Template<float> Tangent3D;
typedef Vector3Template<float> Bitangent3D;
typedef Vector4Template<float> Vector4D;
typedef Normal3Template<float> Normal3D;
typedef Point3Template<float> Point3D;
typedef Matrix4x4Template<float> Matrix4x4;
typedef QuaternionTemplate<float> Quaternion;
typedef TexCoord2Template<float> TexCoord2D;

//
template<typename RealType> struct CompensatedSum;
typedef CompensatedSum<float> FloatSum;
typedef CompensatedSum<Spectrum> SpectrumSum;

// Transforms
class Transform;
class StaticTransform;
class AnimatedTransform;
class ChainedTransform;

typedef std::shared_ptr<Transform> TransformRef;

// Random Number Generators
struct Types32bit;
struct Types64bit;
template <typename TypeSet> class RandomNumberGeneratorTemplate;
template <typename TypeSet> class XORShiftTemplate;

typedef RandomNumberGeneratorTemplate<Types32bit> RandomNumberGenerator;
typedef XORShiftTemplate<Types32bit> XORShift;

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

typedef std::shared_ptr<Image2D> Image2DRef;
typedef std::shared_ptr<TiledImage2D> TiledImage2DRef;

//
class ImageSensor;

// 
struct Ray;
struct BoundingBox3D;
struct Vertex;
struct ReferenceFrame;
struct Intersection;
struct SurfacePoint;

// Surface Objects
class SurfaceObject;
class SingleSurfaceObject;
class InfiniteSphereSurfaceObject;
class SurfaceObjectAggregate;
class TransformedSurfaceObject;

struct LightPosQuery;
struct LightPosQueryResult;
class Light;

struct LensPosQuery;
struct LensPosQueryResult;

// Cameras
class Camera;
class PerspectiveCamera;
class EquirectangularCamera;

// Surfaces
class Surface;
class Triangle;
class InfiniteSphere;

// Nodes
class Node;
class InternalNode;
class SurfaceObjectNode;
class ReferenceNode;
class CameraNode;
class TriangleMeshNode;
class InfiniteSphereNode;
typedef std::shared_ptr<Node> NodeRef;
typedef std::shared_ptr<InternalNode> InternalNodeRef;
typedef std::shared_ptr<SurfaceObjectNode> SurfaceObjectNodeRef;
typedef std::shared_ptr<ReferenceNode> ReferenceNodeRef;
typedef std::shared_ptr<CameraNode> CameraNodeRef;
typedef std::shared_ptr<TriangleMeshNode> TriangleMeshNodeRef;
typedef std::shared_ptr<InfiniteSphereNode> InfiniteSphereNodeRef;

struct RenderingData;

// Accelerators
class BBVH;

// Textures
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
typedef std::shared_ptr<SpectrumTexture> SpectrumTextureRef;
typedef std::shared_ptr<Normal3DTexture> Normal3DTextureRef;
typedef std::shared_ptr<FloatTexture> FloatTextureRef;

// Materials
class SpatialFresnel;
class SpatialFresnelNoOp;
class SpatialFresnelConductor;
class SpatialFresnelDielectric;
class SurfaceMaterial;
class EmitterSurfaceProperty;
class DiffuseReflection;
class SpecularReflection;
class SpecularTransmission;
class ModifiedWardDurReflection;
class AshikhminSpecularReflection;
class AshikhminDiffuseReflection;
class AddedSurfaceMaterial;
class MixedSurfaceMaterial;
class EmitterSurfaceMaterial;
class DiffuseEmission;
class IBLEmission;
typedef std::shared_ptr<SpatialFresnel> SpatialFresnelRef;
typedef std::shared_ptr<SpatialFresnelDielectric> SpatialFresnelDielectricRef;
typedef std::shared_ptr<SurfaceMaterial> SurfaceMaterialRef;
typedef std::shared_ptr<EmitterSurfaceProperty> EmitterSurfacePropertyRef;

// Directional Distribution Functions
struct DirectionType;
struct EDFQuery;
struct EDFSample;
struct EDFQueryResult;
struct BSDFQuery;
struct BSDFSample;
struct BSDFQueryResult;
struct IDFQuery;
struct IDFSample;
struct IDFQueryResult;
class EDF;
class BSDF;
class IDF;
class DiffuseEDF;
class IBLEDF;
class MultiEDF;
class Fresnel;
class FresnelNoOp;
class FresnelConductor;
class FresnelDielectric;
class LambertianBRDF;
class SpecularBRDF;
class SpecularBTDF;
class ModifiedWardDurBRDF;
class AshikhminSpecularBRDF;
class AshikhminDiffuseBRDF;
class MultiBSDF;
class PerspectiveIDF;
class EquirectangularIDF;

class Scene;

class RenderSettings;
class Renderer;

// Renderers
class PathTracingRenderer;
class AMCMCPPMRenderer;

#endif
