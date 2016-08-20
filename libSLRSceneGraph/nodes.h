//
//  nodes.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__entities__
#define __SLRSceneGraph__entities__

#include <libSLR/defines.h>
#include "references.h"

namespace SLRSceneGraph {
    struct SLR_SCENEGRAPH_API RenderingData {
        std::vector<SLR::SurfaceObject*> surfObjs;
        SLR::Camera* camera;
        const SLR::Transform* camTransform;
        RenderingData() : camera(nullptr), camTransform(nullptr) { }
    };
    
    class SLR_SCENEGRAPH_API Node {
    protected:
        std::string m_name;
    public:
        Node() { }
        virtual ~Node() { }
        
        void setName(const std::string &name) { m_name = name; }
        std::string getName() const { return m_name; }
        
        virtual bool contains(const NodeRef &obj) const { return this == obj.get(); }
        virtual bool hasChildren() const { return false; }
        virtual bool isInstanced() const { return false; }
        virtual NodeRef copy() const { SLRAssert_NotImplemented(); return nullptr; }
        
        virtual void applyTransform() { SLRAssert_NotImplemented(); }
        virtual void applyTransform(const SLR::StaticTransform &tf) { SLRAssert_NotImplemented(); }
        virtual void applyTransformToLeaf(const SLR::StaticTransform &tf) { SLRAssert_NotImplemented(); }
        
        virtual void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) { SLRAssert_NotImplemented(); }
    };
    
    class SLR_SCENEGRAPH_API InternalNode : public Node {
        std::vector<NodeRef> m_childNodes;
        TransformRef m_localToWorld;
    public:
        InternalNode() { }
        ~InternalNode() { }
        
        bool addChildNode(const NodeRef &node);
        const NodeRef &childNodeAt(int i) const;
        NodeRef &childNodeAt(int i);
        
        void setTransform(const TransformRef &tf);
        const TransformRef getTransform() const;
        
        bool contains(const NodeRef &obj) const override;
        bool hasChildren() const override;
        NodeRef copy() const override;
        
        void applyTransform() override;
        void applyTransform(const SLR::StaticTransform &tf) override;
        void applyTransformToLeaf(const SLR::StaticTransform &tf) override;

        void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) override;
    };
    
    class SLR_SCENEGRAPH_API SurfaceObjectNode : public Node {
    protected:
        bool m_ready;
        SLR::SingleSurfaceObject** m_refinedObjs;
        size_t m_numRefinedObjs;
    public:
        SurfaceObjectNode() : m_ready(false) { }
        virtual ~SurfaceObjectNode();
        
        virtual void applyTransformForRendering(const SLR::StaticTransform &tf) = 0;
        virtual void createSurfaceObjects() = 0;
        void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) final;
    };
    
    class SLR_SCENEGRAPH_API ReferenceNode : public Node {
        NodeRef m_node;
        bool m_ready;
        RenderingData m_subData;
        SLR::SurfaceObject* m_surfObj;
    public:
        ReferenceNode(const NodeRef node) : m_node(node), m_ready(false) { }
        bool isInstanced() const override { return true; }
        
        void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) final;
    };
    
    class SLR_SCENEGRAPH_API CameraNode : public Node {
    protected:
        bool m_ready;
        SLR::Camera* m_camera;
    public:
        CameraNode() : m_ready(false) { }
        virtual ~CameraNode();
        
        virtual void createCamera() = 0;
        void getRenderingData(SLR::ArenaAllocator &mem, const SLR::Transform* subTF, RenderingData *data) final;
    };
    
}

#endif
