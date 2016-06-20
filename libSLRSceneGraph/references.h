//
//  references.h
//
//  Created by 渡部 心 on 2015/09/27.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_references_h
#define SLRSceneGraph_references_h

#include <libSLR/references.h>

#ifdef SLR_Defs_Windows
#   ifdef SLR_SCENEGRAPH_API_EXPORTS
#       define SLR_SCENEGRAPH_API __declspec(dllexport)
#   else
#       define SLR_SCENEGRAPH_API __declspec(dllimport)
#   endif
#else
#   define SLR_SCENEGRAPH_API
#endif

namespace SLRSceneGraph {
    typedef std::shared_ptr<SLR::Transform> TransformRef;
    
    typedef std::shared_ptr<SLR::InputSpectrum> InputSpectrumRef;
    
    typedef std::shared_ptr<SLR::Image2D> Image2DRef;
    typedef std::shared_ptr<SLR::TiledImage2D> TiledImage2DRef;
    
    class Texture2DMapping;
    class Texture3DMapping;
    class SpectrumTexture;
    class Normal3DTexture;
    class FloatTexture;
    
    typedef std::shared_ptr<Texture2DMapping> Texture2DMappingRef;
    typedef std::shared_ptr<Texture3DMapping> Texture3DMappingRef;
    typedef std::shared_ptr<SpectrumTexture> SpectrumTextureRef;
    typedef std::shared_ptr<Normal3DTexture> Normal3DTextureRef;
    typedef std::shared_ptr<FloatTexture> FloatTextureRef;
    
    class SurfaceMaterial;
    class EmitterSurfaceProperty;
    class SpatialFresnel;
    
    typedef std::shared_ptr<SpatialFresnel> SpatialFresnelRef;
    typedef std::shared_ptr<SurfaceMaterial> SurfaceMaterialRef;
    typedef std::shared_ptr<EmitterSurfaceProperty> EmitterSurfacePropertyRef;
    
    class IBLEmission;
    typedef std::shared_ptr<IBLEmission> IBLEmissionRef;
    
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
    
    class Scene;
    typedef std::shared_ptr<Scene> SceneRef;
    typedef std::weak_ptr<Scene> SceneWRef;
    
    struct RenderingContext;
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
}

#endif
