//
//  references.h
//
//  Created by 渡部 心 on 2015/09/27.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_references_h
#define SLRSceneGraph_references_h

#include <libSLR/references.h>

namespace SLRSceneGraph {
    typedef std::shared_ptr<SLR::Transform> TransformRef;
    
    typedef std::shared_ptr<SLR::InputSpectrum> InputSpectrumRef;
    
    typedef std::shared_ptr<SLR::Image2D> Image2DRef;
    typedef std::shared_ptr<SLR::TiledImage2D> TiledImage2DRef;
    
    class SpectrumTexture;
    class Normal3DTexture;
    class FloatTexture;
    
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
    typedef std::weak_ptr<Scene> SceneWRef;
}

#endif
