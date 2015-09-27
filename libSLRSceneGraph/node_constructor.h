//
//  node_constructor.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLRSceneGraph__node_constructor__
#define __SLRSceneGraph__node_constructor__

#include <libSLR/defines.h>
#include <libSLR/references.h>
#include "references.h"
#include <assimp/scene.h>

namespace SLRSceneGraph {
    typedef SLR::SurfaceMaterialRef (*createMaterialFunction)(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem);
    SLR::SurfaceMaterialRef createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem);
    void construct(const aiScene &objSrc, const std::string &pathPrefix, InternalNodeRef &nodeOut, const createMaterialFunction materialFunc = createMaterialDefaultFunction);
}

#endif
