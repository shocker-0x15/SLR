//
//  declarations.h
//
//  Created by 渡部 心 on 2015/09/27.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_references__
#define __SLRSceneGraph_references__

#include <libSLR/declarations.h>

#ifdef SLR_Platform_Windows_MSVC
#   ifdef SLR_SCENEGRAPH_API_EXPORTS
#       define SLR_SCENEGRAPH_API __declspec(dllexport)
#   else
#       define SLR_SCENEGRAPH_API __declspec(dllimport)
#   endif
#else
#   define SLR_SCENEGRAPH_API
#endif

namespace SLRSceneGraph {
    struct SceneParsingDriver;
    
    class Statement;
    typedef std::shared_ptr<Statement> StatementRef;
    typedef std::shared_ptr<std::vector<StatementRef>> StatementsRef;
    class Expression;
    typedef std::shared_ptr<Expression> ExpressionRef;
    class Term;
    typedef std::shared_ptr<Term> TermRef;
    class SingleTerm;
    typedef std::shared_ptr<SingleTerm> SingleTermRef;
    class Value;
    typedef std::shared_ptr<Value> ValueRef;
    class ArgumentDefinition;
    typedef std::shared_ptr<ArgumentDefinition> ArgumentDefinitionRef;
    typedef std::shared_ptr<std::vector<ArgumentDefinitionRef>> ArgumentDefinitionVecRef;
    class Parameter;
    typedef std::shared_ptr<Parameter> ParameterRef;
    typedef std::shared_ptr<std::vector<ParameterRef>> ParameterVecRef;
    
    struct ErrorMessage;
    
    struct Element;
    struct ParameterList;
    typedef std::shared_ptr<ParameterList> ParameterListRef;
    class Function;
    struct TypeInfo;
    
    class LocalVariables;
    struct ExecuteContext;
    
    
    
    typedef std::shared_ptr<SLR::AssetSpectrum> AssetSpectrumRef;
    
    typedef std::shared_ptr<SLR::Transform> TransformRef;
    
    typedef std::shared_ptr<SLR::Image2D> Image2DRef;
    typedef std::shared_ptr<SLR::TiledImage2D> TiledImage2DRef;
    
    class Texture2DMapping;
    class Texture3DMapping;
    class SpectrumTexture;
    class NormalTexture;
    class FloatTexture;
    
    typedef std::shared_ptr<Texture2DMapping> Texture2DMappingRef;
    typedef std::shared_ptr<Texture3DMapping> Texture3DMappingRef;
    typedef std::shared_ptr<SpectrumTexture> SpectrumTextureRef;
    typedef std::shared_ptr<NormalTexture> NormalTextureRef;
    typedef std::shared_ptr<FloatTexture> FloatTextureRef;
    
    class SVFresnel;
    class SVMicrofacetDistribution;
    class SurfaceMaterial;
    class EmitterSurfaceProperty;
    class MediumMaterial;
    class EmitterMediumProperty;
    
    typedef std::shared_ptr<SVFresnel> SVFresnelRef;
    typedef std::shared_ptr<SVMicrofacetDistribution> SVMicrofacetDistributionRef;
    typedef std::shared_ptr<SurfaceMaterial> SurfaceMaterialRef;
    typedef std::shared_ptr<EmitterSurfaceProperty> EmitterSurfacePropertyRef;
    typedef std::shared_ptr<MediumMaterial> MediumMaterialRef;
    typedef std::shared_ptr<EmitterMediumProperty> EmitterMediumPropertyRef;
    
    class IBLEmitterSurfaceProperty;
    typedef std::shared_ptr<IBLEmitterSurfaceProperty> IBLEmitterSurfacePropertyRef;
    
    // Nodes
    class Node;
    class InternalNode;
    class ReferenceNode;
    class SurfaceNode;
    class MediumNode;
    class InfinitesimalPointNode;
    class InfiniteSphereNode;
    class PerspectiveCameraNode;
    class EquirectangularCameraNode;
    class TriangleMeshNode;
    class HomogeneousMediumNode;
    class GridMediumNode;
    class DensityGridMediumNode;
    typedef std::shared_ptr<Node> NodeRef;
    typedef std::shared_ptr<InternalNode> InternalNodeRef;
    typedef std::shared_ptr<ReferenceNode> ReferenceNodeRef;
    typedef std::shared_ptr<SurfaceNode> SurfaceNodeRef;
    typedef std::shared_ptr<MediumNode> MediumNodeRef;
    typedef std::shared_ptr<InfiniteSphereNode> InfiniteSphereNodeRef;
    typedef std::shared_ptr<TriangleMeshNode> TriangleMeshNodeRef;
    
    class Scene;
    typedef std::shared_ptr<Scene> SceneRef;
    typedef std::weak_ptr<Scene> SceneWRef;
    
    struct RenderingContext;
}

#endif /* __SLRSceneGraph_references__ */
