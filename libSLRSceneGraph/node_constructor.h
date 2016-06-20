//
//  node_constructor.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__node_constructor__
#define __SLRSceneGraph__node_constructor__

#include <libSLR/defines.h>
#include "references.h"
#include <assimp/scene.h>

namespace SLRSceneGraph {
    struct SurfaceAttributeTuple {
        SurfaceMaterialRef material;
        Normal3DTextureRef normalMap;
        FloatTextureRef alphaMap;
        SurfaceAttributeTuple(const SurfaceMaterialRef &mat, const Normal3DTextureRef &norm, const FloatTextureRef &alpha) :
        material(mat), normalMap(norm), alphaMap(alpha) {}
    };
    
    typedef std::function<SurfaceAttributeTuple(const aiMaterial*, const std::string &, SLR::Allocator*)> CreateMaterialFunction;
    SLR_SCENEGRAPH_API SurfaceAttributeTuple createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem);
    SLR_SCENEGRAPH_API void construct(const std::string &filePath, InternalNodeRef &nodeOut,
                                      const CreateMaterialFunction &materialFunc = createMaterialDefaultFunction);
}

#endif
