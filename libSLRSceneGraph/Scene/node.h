//
//  node.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_node__
#define __SLRSceneGraph_node__

#include <libSLR/defines.h>
#include <libSLR/Core/geometry.h>
#include "../declarations.h"

namespace SLRSceneGraph {    
    class SLR_SCENEGRAPH_API Node {
    protected:
        SLR::Node* m_rawData;
        bool m_setup;
        std::string m_name;
        
        virtual void allocateRawData() = 0;
        virtual void setupRawData() = 0;
        virtual void terminateRawData() = 0;
    public:
        Node() : m_rawData(nullptr), m_setup(false) { }
        virtual ~Node();

        SLR::Node* getRaw() const {
            return m_rawData;
        }
        
        void setName(const std::string &name) { m_name = name; }
        std::string getName() const { return m_name; }
        
        virtual bool isUniqueInTree() const { return true; }
        virtual bool contains(const NodeRef &obj) const { return this == obj.get(); }
        virtual bool hasChildren() const { return false; }
        virtual NodeRef copy() const = 0;
        
        virtual void applyTransform(const SLR::StaticTransform &tf) { SLRAssert_NotImplemented(); }
        virtual void applyTransformToLeaf(const SLR::StaticTransform &tf) { SLRAssert_NotImplemented(); }
        
        virtual void prepareForRendering() = 0;
    };
    
    
    
    class SLR_SCENEGRAPH_API InternalNode : public Node {
        std::vector<NodeRef> m_childNodes;
        TransformRef m_localToWorld;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        InternalNode(const TransformRef &localToWorld);
        
        bool addChildNode(const NodeRef &node);
        const NodeRef &childNodeAt(int i) const;
        NodeRef &childNodeAt(int i);
        
        void setTransform(const TransformRef &tf);
        const TransformRef getTransform() const;
        
        bool contains(const NodeRef &obj) const override;
        bool hasChildren() const override;
        NodeRef copy() const override;

        void applyTransform(const SLR::StaticTransform &tf) override;
        void applyTransformToLeaf(const SLR::StaticTransform &tf) override;
        
        void prepareForRendering() override;
        
        void propagateTransform();
    };
    
    
    
    class SLR_SCENEGRAPH_API ReferenceNode : public Node {
        NodeRef m_node;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        ReferenceNode(const NodeRef &node);
        
        bool isUniqueInTree() const override { return false; }
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    
    class SLR_SCENEGRAPH_API SurfaceNode : public Node {
    protected:
        MediumNodeRef m_enclosedMediumNode;
    public:
        void setInternalMedium(const MediumNodeRef &medium) {
            m_enclosedMediumNode = medium;
        }
    };
    
    
    
    class SLR_SCENEGRAPH_API MediumNode : public Node {
    public:
    };
    
    
    
    class SLR_SCENEGRAPH_API InfinitesimalPointNode : public SurfaceNode {
        SLR::Point3D m_position;
        SLR::Vector3D m_direction;
        SurfaceMaterialRef m_material;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        InfinitesimalPointNode(const SLR::Point3D &p, const SLR::Vector3D &d, const SurfaceMaterialRef &mat);
        
        NodeRef copy() const override;
        
        void prepareForRendering() override;
    };
    
    
    
    class SLR_SCENEGRAPH_API InfiniteSphereNode : public SurfaceNode {
        SceneWRef m_scene;
        SpectrumTextureRef m_IBLTex;
        float m_scale;
        
        void allocateRawData() override;
        void setupRawData() override;
        void terminateRawData() override;
    public:
        InfiniteSphereNode(const SceneWRef &scene, const SpectrumTextureRef &IBLTex, float scale);
        
        NodeRef copy() const override { SLRAssert_ShouldNotBeCalled(); return nullptr; }
        
        void prepareForRendering() override;
    };
}

#endif /* __SLRSceneGraph_node__ */
