//
//  nodes.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "nodes.h"
#include "SurfaceObject.h"
#include "../Memory/ArenaAllocator.h"

bool InternalNode::addChildNode(const NodeRef &node) {
    // Create a shared_ptr for passing it to contains() with a No-Op deleter so that "this" will not be deleted.
    NodeRef thisRef(this, [](void* ptr){});
    if (node->contains(thisRef) || this == node.get()) {
        printf("This causes a circular reference.\n");
        return false;
    }
    if (node->isInstanced()) {
        auto it = std::find(m_childNodes.begin(), m_childNodes.end(), node);
        if (it != m_childNodes.end()) {
            printf("Another instanced node already exists in this node.\n");
            return false;
        }
    }
    else {
        if (contains(node)) {
            printf("This node already has the given node.\n");
            return false;
        }
    }
    m_childNodes.push_back(node);
    node->setParent(this);
    return true;
}

const NodeRef &InternalNode::childNodeAt(int i) const {
    return m_childNodes[i];
}

NodeRef &InternalNode::childNodeAt(int i) {
    return m_childNodes[i];
}


void InternalNode::setTransform(const TransformRef &tf) {
    m_localToWorld = tf;
}

const TransformRef InternalNode::getTransform() const {
    return m_localToWorld;
}


bool InternalNode::contains(const NodeRef &obj) const {
    for (int i = 0; i < m_childNodes.size(); ++i) {
        if (m_childNodes[i] == obj || m_childNodes[i]->contains(obj))
            return true;
    }
    return false;
}

bool InternalNode::hasChildren() const {
    return m_childNodes.size() > 0;
}


void InternalNode::applyTransform() {
    TransformRef tf = getTransform();
    if (tf->isStatic()) {
        for (int i = 0; i < m_childNodes.size(); ++i)
            m_childNodes[i]->applyTransform(*(StaticTransform*)tf.get());
        m_localToWorld = createShared<StaticTransform>(Matrix4x4::Identity);
    }
    else {
        SLRAssert(false, "AnimatedTransform case is not implemented.");
    }
}

void InternalNode::applyTransform(const StaticTransform &tf) {
    SLRAssert(false, "Not implemented.");
//    m_localToWorld = tf * m_localToWorld;
}

void InternalNode::applyTransformToLeaf(const StaticTransform &tf) {
    SLRAssert(false, "Not implemented.");
//    StaticTransform cumT = tf * m_localToWorld;
//    for (int i = 0; i < m_childNodes.size(); ++i) {
//        NodeRef &child = m_childNodes[i];
//        if (child->hasChildren())
//            child->applyTransformToLeaf(cumT);
//        else
//            child->applyTransform(cumT);
//    }
//    m_localToWorld = Matrix4x4::Identity;
}

void InternalNode::getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) {
    const Transform* reduced;
    if (subTF)
        reduced = ChainedTransform(subTF, m_localToWorld.get()).reduce(mem);// &m_localToWorld or m_localToWorld.copy(tfMem)?
    else
        reduced = m_localToWorld.get();
    
    if (m_localToWorld->isStatic()) {
        for (int i = 0; i < m_childNodes.size(); ++i)
            m_childNodes[i]->getRenderingData(mem, reduced, data);
    }
    else {
        RenderingData subData;
        for (int i = 0; i < m_childNodes.size(); ++i)
            m_childNodes[i]->getRenderingData(mem, nullptr, &subData);
        if (subData.surfObjs.size() > 1) {
            SurfaceObjectAggregate* aggr = mem.create<SurfaceObjectAggregate>(subData.surfObjs);
            TransformedSurfaceObject* tfobj = mem.create<TransformedSurfaceObject>(aggr, reduced);
            data->surfObjs.push_back(tfobj);
        }
        else if (subData.surfObjs.size() == 1) {
            TransformedSurfaceObject* tfobj = mem.create<TransformedSurfaceObject>(subData.surfObjs[0], reduced);
            data->surfObjs.push_back(tfobj);
        }
        
        if (subData.camera) {
            data->camera = subData.camera;
            data->camParent = subData.camParent;
        }
    }
}


SurfaceObjectNode::~SurfaceObjectNode() {
    if (!m_ready)
        return;
    for (int i = 0; i < m_numRefinedObjs; ++i)
        delete m_refinedObjs[i];
    delete[] m_refinedObjs;
}

void SurfaceObjectNode::getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) {
    if (!m_ready) {
        createSurfaceObjects();
        m_ready = true;
    }
    size_t curSize = data->surfObjs.size();
    data->surfObjs.resize(curSize + m_numRefinedObjs);
    
    StaticTransform transform;
    if (subTF) {
        SLRAssert(subTF->isStatic(), "Transformation given to SurfaceObjectNode must be static.");
        subTF->sample(0.0f, &transform);
    }
    if (transform.isIdentity()) {
        for (int i = 0; i < m_numRefinedObjs; ++i) {
            SingleSurfaceObject* obj = m_refinedObjs[i];
            data->surfObjs[curSize + i] = obj;
        }
    }
    else {
        for (int i = 0; i < m_numRefinedObjs; ++i) {
            SingleSurfaceObject* obj = m_refinedObjs[i];
            TransformedSurfaceObject* tfObj = mem.create<TransformedSurfaceObject>(obj, subTF);
            data->surfObjs[curSize + i] = tfObj;
        }
    }
}


void ReferenceNode::getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) {
    if (!m_ready) {
        m_node->getRenderingData(mem, nullptr, &m_subData);
        if (m_subData.surfObjs.size() > 1)
            m_surfObj = mem.create<SurfaceObjectAggregate>(m_subData.surfObjs);
        else
            m_surfObj = m_subData.surfObjs[0];
        m_ready = true;
    }
    data->surfObjs.push_back(mem.create<TransformedSurfaceObject>(m_surfObj, subTF));
}


CameraNode::~CameraNode() {
    if (!m_ready)
        return;
    delete m_camera;
}

void CameraNode::getRenderingData(ArenaAllocator &mem, const Transform* subTF, RenderingData *data) {
    if (!m_ready) {
        createCamera();
        m_ready = true;
    }
    data->camera = m_camera;
    data->camParent = m_parent;
}
