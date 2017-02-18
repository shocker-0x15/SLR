//
//  nodes.h
//
//  Created by 渡部 心 on 2017/01/05.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_nodes__
#define __SLR_nodes__

#include "../defines.h"
#include "../references.h"
#include "../BasicTypes/Point3.h"
#include "../BasicTypes/Vector3.h"

namespace SLR {
    struct SLR_API RenderingData {
        Scene* const scene;
        std::vector<SurfaceObject*> surfObjs;
        std::vector<MediumObject*> medObjs;
        Camera* camera;
        const Transform* camTransform;
        InfiniteSphereSurfaceObject* envObj;
        
        RenderingData(Scene* sc) :
        scene(sc), camera(nullptr), camTransform(nullptr), envObj(nullptr) { }
    };
    
    
    
    class SLR_API Node {
    public:
        virtual ~Node() { }
        
        virtual bool isDirectlyTransformable() const = 0;
        virtual void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) = 0;
        virtual void destroyRenderingData(Allocator* mem) = 0;
    };
    
    
    
    class SLR_API InternalNode : public Node {
        std::vector<Node*> m_childNodes;
        const Transform* m_localToWorld;
        
        Transform* m_appliedTransform;
        SurfaceObjectAggregate* m_subSurfObj;
        MediumObjectAggregate* m_subMedObj;
        std::vector<TransformedSurfaceObject*> m_TFSurfObjs;
        std::vector<TransformedMediumObject*> m_TFMedObjs;
    public:
        InternalNode(const Transform* localToWorld) :
        m_localToWorld(localToWorld),
        m_appliedTransform(nullptr), m_subSurfObj(nullptr), m_subMedObj(nullptr) { }
        
        bool isDirectlyTransformable() const override { return true; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
        
        void addChildNode(Node* node);
        void setTransform(const Transform* localToWorld) {
            m_localToWorld = localToWorld;
        }
    };
    
    
    
    class SLR_API ReferenceNode : public Node {
        Node* m_node;
        
        bool m_ready;
        SurfaceObject* m_obj;
        bool m_isAggregate;
    public:
        ReferenceNode(Node* node) :
        m_node(node),
        m_ready(false), m_obj(nullptr) { }
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API SurfaceNode : public Node {
    protected:
        MediumNode* m_enclosedMediumNode;
        SurfaceObject* m_boundarySurfObj;
        MediumObject* m_enclosedMedObj;
        Transform* m_mediumTransform;
        TransformedMediumObject* m_TFMedObj;
    public:
        SurfaceNode() : m_enclosedMediumNode(nullptr), m_boundarySurfObj(nullptr), m_enclosedMedObj(nullptr), m_mediumTransform(nullptr), m_TFMedObj(nullptr) { }
        
        void setInternalMedium(MediumNode* medium) {
            m_enclosedMediumNode = medium;
        }
    };
    
    
    
    class SLR_API MediumNode : public Node {
    public:
    };
    
    
    
    class SLR_API InfiniteSphereNode : public Node {
        const Scene* m_scene;
        const SpectrumTexture* m_IBLTex;
        float m_scale;
        
        IBLEmitterSurfaceProperty* m_emission;
        InfiniteSphereSurfaceObject* m_obj;
    public:
        InfiniteSphereNode(const Scene* scene, const SpectrumTexture* IBLTex, float scale) :
        m_scene(scene), m_IBLTex(IBLTex), m_scale(scale),
        m_emission(nullptr), m_obj(nullptr) { }
        
        bool isDirectlyTransformable() const override { return false; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
    
    
    
    class SLR_API InfinitesimalPointNode : public SurfaceNode {
        const Point3D m_position;
        const Vector3D m_direction;
        const SurfaceMaterial* m_material;
        InfinitesimalPointSurfaceShape* m_surface;
        SingleSurfaceObject* m_obj;
    public:
        InfinitesimalPointNode(const Point3D &p, const Vector3D &d, const SurfaceMaterial* mat) : m_position(p), m_direction(d), m_material(mat), m_surface(nullptr), m_obj(nullptr) { }
        
        bool isDirectlyTransformable() const override { return true; }
        void createRenderingData(Allocator* mem, const Transform* subTF, RenderingData *data) override;
        void destroyRenderingData(Allocator* mem) override;
    };
}

#endif /* __SLR_nodes__ */
