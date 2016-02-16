//
//  node_constructor.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "node_constructor.h"
#include "TriangleMeshNode.h"
#include "API.hpp"
#include "surface_materials.hpp"
#include "textures.hpp"
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <libSLR/Memory/Allocator.h>
#include <libSLR/Core/Transform.h>

template <typename RealType>
inline void makeTangent(RealType nx, RealType ny, RealType nz, RealType* s) {
    if (std::fabs(nx) > std::fabs(ny)) {
        RealType invLen = 1.0f / std::sqrt(nx * nx + nz * nz);
        s[0] = -nz * invLen;
        s[1] = 0.0f;
        s[2] = nx * invLen;
    }
    else {
        RealType invLen = 1.0f / std::sqrt(ny * ny + nz * nz);
        s[0] = 0.0f;
        s[1] = nz * invLen;
        s[2] = -ny * invLen;
    }
}

namespace SLRSceneGraph {
    static void recursiveConstruct(const aiScene* objSrc, const aiNode* nodeSrc, const std::vector<SurfaceMaterialRef> &materials, InternalNodeRef &nodeOut) {
        if (nodeSrc->mNumMeshes == 0 && nodeSrc->mNumChildren == 0) {
            nodeOut = nullptr;
            return;
        }
        nodeOut = createShared<InternalNode>();
        nodeOut->setName(nodeSrc->mName.C_Str());
        
        const aiMatrix4x4 &tf = nodeSrc->mTransformation;
        float tfElems[] = {
            tf.a1, tf.a2, tf.a3, tf.a4,
            tf.b1, tf.b2, tf.b3, tf.b4,
            tf.c1, tf.c2, tf.c3, tf.c4,
            tf.d1, tf.d2, tf.d3, tf.d4,
        };
        nodeOut->setTransform(createShared<SLR::StaticTransform>(SLR::Matrix4x4(tfElems)));
        
        for (int m = 0; m < nodeSrc->mNumMeshes; ++m) {
            const aiMesh* mesh = objSrc->mMeshes[nodeSrc->mMeshes[m]];
            if (mesh->mPrimitiveTypes != aiPrimitiveType_TRIANGLE) {
                printf("ignored non triangle mesh.\n");
                continue;
            }
            
            TriangleMeshNodeRef surfMesh = createShared<TriangleMeshNode>();
            const SurfaceMaterialRef &surfMat = materials[mesh->mMaterialIndex];
            
            for (int v = 0; v < mesh->mNumVertices; ++v) {
                const aiVector3D &p = mesh->mVertices[v];
                const aiVector3D &n = mesh->mNormals[v];
                float tangent[3];
                if (mesh->mTangents == nullptr)
                    makeTangent(n.x, n.y, n.z, tangent);
                const aiVector3D &t = mesh->mTangents ? mesh->mTangents[v] : aiVector3D(tangent[0], tangent[1], tangent[2]);
                const aiVector3D &uv = mesh->mNumUVComponents[0] > 0 ? mesh->mTextureCoords[0][v] : aiVector3D(0, 0, 0);
                
                SLR::Vertex outVtx{SLR::Point3D(p.x, p.y, p.z), SLR::Normal3D(n.x, n.y, n.z), SLR::Tangent3D(t.x, t.y, t.z), SLR::TexCoord2D(uv.x, uv.y)};
                float dotNT = dot(outVtx.normal, outVtx.tangent);
                if (std::fabs(dotNT) >= 0.01f)
                    outVtx.tangent = SLR::normalize(outVtx.tangent - dotNT * outVtx.normal);
                //SLRAssert(absDot(outVtx.normal, outVtx.tangent) < 0.01f, "shading normal and tangent must be orthogonal: %g", absDot(outVtx.normal, outVtx.tangent));
                surfMesh->addVertex(outVtx);
            }
            
            for (int f = 0; f < mesh->mNumFaces; ++f) {
                const aiFace &face = mesh->mFaces[f];
                surfMesh->addTriangle(face.mIndices[0], face.mIndices[1], face.mIndices[2], surfMat);
            }
            
            surfMesh->setName(mesh->mName.C_Str());
            nodeOut->addChildNode(surfMesh);
        }
        
        if (nodeSrc->mNumChildren) {
            for (int c = 0; c < nodeSrc->mNumChildren; ++c) {
                InternalNodeRef subNode;
                recursiveConstruct(objSrc, nodeSrc->mChildren[c], materials, subNode);
                if (subNode != nullptr)
                    nodeOut->addChildNode(subNode);
            }
        }
    }
    
    SLR_SCENEGRAPH_API SurfaceMaterialRef createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem) {
        using namespace SLR;
        aiReturn ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        SpectrumTextureRef diffuseTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = Image::createTiledImage((pathPrefix + strValue.C_Str()).c_str(), mem, SpectrumType::Reflectance);
            diffuseTex = createShared<ImageSpectrumTexture>(image);
        }
        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
            InputSpectrumRef sp = Spectrum::create(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, color[0], color[1], color[2]);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        else {
            InputSpectrumRef sp = Spectrum::create(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        
        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
    };
    
    SLR_SCENEGRAPH_API void construct(const std::string &filePath, InternalNodeRef &nodeOut, const CreateMaterialFunction &materialFunc) {
        using namespace SLR;
        DefaultAllocator &defMem = DefaultAllocator::instance();
        
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(filePath, 0);
        if (!scene) {
            printf("Failed to load %s.\n", filePath.c_str());
            return;
        }
        printf("Reading: %s done.\n", filePath.c_str());
        
        std::string pathPrefix = filePath.substr(0, filePath.find_last_of("/") + 1);
        
        // マテリアルの生成。
        std::vector<SurfaceMaterialRef> materials;
        std::vector<Normal3DTextureRef> normalMaps;
        for (int m = 0; m < scene->mNumMaterials; ++m) {
            const aiMaterial* aiMat = scene->mMaterials[m];
            SurfaceMaterialRef surfMat = materialFunc(aiMat, pathPrefix, &defMem);
            materials.push_back(surfMat);
            
            aiString strValue;
            Normal3DTextureRef normalTex;
            if (aiMat->Get(AI_MATKEY_TEXTURE_DISPLACEMENT(0), strValue) == aiReturn_SUCCESS) {
                TiledImage2DRef image = Image::createTiledImage((pathPrefix + strValue.C_Str()).c_str(), &defMem, SpectrumType::Illuminant);
                normalTex = createShared<ImageNormal3DTexture>(image);
            }
            normalMaps.push_back(normalTex);
        }
        
        recursiveConstruct(scene, scene->mRootNode, materials, nodeOut);
        
        nodeOut->setName(filePath);
        
        printf("Constructing: %s done.\n", filePath.c_str());
    }
}
