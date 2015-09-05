//
//  nodes.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__entities__
#define __SLR__entities__

#include "../defines.h"
#include "../references.h"
#include "../Core/Transform.h"

struct RenderingData {
    std::vector<SurfaceObject*> surfObjs;
    Camera* camera;
    InternalNode* camParent;
    RenderingData() : camera(nullptr), camParent(nullptr) { };
};

class Node {
protected:
    std::string m_name;
    InternalNode* m_parent;
public:
    Node() : m_parent(nullptr) { };
    virtual ~Node() { };
    
    void setName(const std::string &name) { m_name = name; };
    std::string getName() const { return m_name; };
    
    void setParent(InternalNode* parent) { m_parent = parent; };
    InternalNode* getParent() const { return m_parent; };
    
    virtual bool contains(const NodeRef &obj) const { return this == obj.get(); };
    virtual bool hasChildren() const { return false; };
    virtual bool isInstanced() const { return false; };
    
    virtual void applyTransform() { SLRAssert(false, "Not implemented."); };
    virtual void applyTransform(const StaticTransform &tf) { SLRAssert(false, "Not implemented."); };
    virtual void applyTransformToLeaf(const StaticTransform &tf) { SLRAssert(false, "Not implemented."); };
    
    virtual void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) { SLRAssert(false, "Not implemented."); };
};

class InternalNode : public Node {
    std::vector<NodeRef> m_childNodes;
    TransformRef m_localToWorld;
public:
    InternalNode() { };
    ~InternalNode() { };
    
    bool addChildNode(const NodeRef &node);
    const NodeRef &childNodeAt(int i) const;
    NodeRef &childNodeAt(int i);
    
    void setTransform(const TransformRef &tf);
    const TransformRef getTransform() const;
    
    bool contains(const NodeRef &obj) const override;
    bool hasChildren() const override;
    
    void applyTransform() override;
    void applyTransform(const StaticTransform &tf) override;
    void applyTransformToLeaf(const StaticTransform &tf) override;
    void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) override;
};

class SurfaceObjectNode : public Node {
protected:
    bool m_ready;
    SingleSurfaceObject** m_refinedObjs;
    size_t m_numRefinedObjs;
public:
    SurfaceObjectNode() : m_ready(false) { };
    virtual ~SurfaceObjectNode();
    
    virtual void createSurfaceObjects() = 0;
    void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) final;
};

class ReferenceNode : public Node {
    NodeRef m_node;
    bool m_ready;
    RenderingData m_subData;
    SurfaceObject* m_surfObj;
public:
    ReferenceNode(const NodeRef node) : m_node(node), m_ready(false) { };
    bool isInstanced() const override { return true; };
    
    void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) final;
};

class CameraNode : public Node {
protected:
    bool m_ready;
    Camera* m_camera;
public:
    CameraNode() : m_ready(false) { };
    virtual ~CameraNode();
    
    virtual void createCamera() = 0;
    void getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) final;
};

#endif
