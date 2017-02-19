//
//  node_constructor.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph_node_constructor__
#define __SLRSceneGraph_node_constructor__

#include <libSLR/defines.h>
#include <assimp/scene.h>
#include "declarations.h"

namespace SLRSceneGraph {
    struct SurfaceAttributeTuple {
        SurfaceMaterialRef material;
        NormalTextureRef normalMap;
        FloatTextureRef alphaMap;
        SurfaceAttributeTuple(const SurfaceMaterialRef &mat, const NormalTextureRef &norm, const FloatTextureRef &alpha) :
        material(mat), normalMap(norm), alphaMap(alpha) {}
    };
    
    typedef std::function<SurfaceAttributeTuple(const aiMaterial*, const std::string &, SLR::Allocator*)> CreateMaterialFunction;
    SurfaceAttributeTuple createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem);
    SurfaceAttributeTuple createMaterialFunction(const Function &userMatProc, ExecuteContext &context, ErrorMessage* err,
                                                 const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem);
    
    typedef std::function<bool(const std::string &, const TriangleMeshNodeRef &, const SLR::Point3D &, const SLR::Point3D &)> MeshCallback;
    bool meshCallbackDefaultFunction(const std::string &name, const TriangleMeshNodeRef &mesh, const SLR::Point3D &minP, const SLR::Point3D &maxP);
    bool meshCallbackFunction(const Function &meshProc, ExecuteContext &context, ErrorMessage* err,
                              const std::string &name, const TriangleMeshNodeRef &mesh, const SLR::Point3D &minP, const SLR::Point3D &maxP);
    
    SLR_SCENEGRAPH_API void construct(const std::string &filePath, InternalNodeRef &nodeOut,
                                      const CreateMaterialFunction &materialFunc = createMaterialDefaultFunction,
                                      const MeshCallback &meshCallback = meshCallbackDefaultFunction);
}

#endif /* __SLRSceneGraph_node_constructor__ */
