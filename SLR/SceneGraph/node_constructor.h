//
//  node_constructor.h
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__node_constructor__
#define __SLR__node_constructor__

#include "../defines.h"
#include "../references.h"
#include <assimp/scene.h>

namespace Assimp {
    typedef SurfaceMaterialRef (*createMaterialFunction)(const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem);
    SurfaceMaterialRef createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem);
    void construct(const aiScene &objSrc, const std::string &pathPrefix, InternalNodeRef &nodeOut, const createMaterialFunction materialFunc = createMaterialDefaultFunction);
}

#endif
