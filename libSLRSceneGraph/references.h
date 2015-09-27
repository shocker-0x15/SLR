//
//  references.h
//
//  Created by 渡部 心 on 2015/09/27.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_references_h
#define SLRSceneGraph_references_h

#include <vector>
#include <memory>
#include <stdint.h>

namespace SLR {
    class Transform;
}

namespace SLRSceneGraph {
    typedef std::shared_ptr<SLR::Transform> TransformRef;
    
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
}

#endif
