//
//  node_constructor.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "node_constructor.h"
#include "../Memory/Allocator.h"
#include "../Core/Image.h"
#include "../Core/textures.h"
#include "../Core/surface_material.h"
#include "../Textures/constant_textures.h"
#include "../Textures/image_textures.h"
#include "TriangleMeshNode.h"

inline void makeTangent(float nx, float ny, float nz, float* s) {
    if (fabsf(nx) > fabsf(ny)) {
        float invLen = 1.0f / sqrtf(nx * nx + nz * nz);
        s[0] = -nz * invLen;
        s[1] = 0.0f;
        s[2] = nx * invLen;
    }
    else {
        float invLen = 1.0f / sqrtf(ny * ny + nz * nz);
        s[0] = 0.0f;
        s[1] = nz * invLen;
        s[2] = -ny * invLen;
    }
}

namespace Assimp {
    void recursiveConstruct(const aiScene &objSrc, const aiNode* nodeSrc, const std::vector<SurfaceMaterialRef> &materials, InternalNodeRef &nodeOut) {
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
        nodeOut->setTransform(createShared<StaticTransform>(Matrix4x4(tfElems)));
        
        for (int m = 0; m < nodeSrc->mNumMeshes; ++m) {
            const aiMesh* mesh = objSrc.mMeshes[nodeSrc->mMeshes[m]];
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
                
                Vertex outVtx{Point3D(p.x, p.y, p.z), Normal3D(n.x, n.y, n.z), Tangent3D(t.x, t.y, t.z), TexCoord2D(uv.x, uv.y)};
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
    
    SurfaceMaterialRef createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem) {
        aiReturn ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        SpectrumTextureRef diffuseTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), mem);
            diffuseTex = createShared<ImageSpectrumTexture>(image);
        }
        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
            InputSpectrumRef sp = createShared<UpsampledContinuousSpectrum>(ColorSpace::sRGB_NonLinear, color[0], color[1], color[2]);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        else {
            InputSpectrumRef sp = createShared<UpsampledContinuousSpectrum>(ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        
        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
    };
    
    void construct(const aiScene &objSrc, const std::string &pathPrefix, InternalNodeRef &nodeOut, const createMaterialFunction materialFunc) {
        DefaultAllocator &defMem = DefaultAllocator::instance();
        
        // マテリアルの生成。
        std::vector<SurfaceMaterialRef> materials;
        std::vector<Normal3DTextureRef> normalMaps;
        for (int m = 0; m < objSrc.mNumMaterials; ++m) {
            const aiMaterial* aiMat = objSrc.mMaterials[m];
            SurfaceMaterialRef surfMat = materialFunc(aiMat, pathPrefix, &defMem);
            materials.push_back(surfMat);
            
            aiString strValue;
            Normal3DTextureRef normalTex;
            if (aiMat->Get(AI_MATKEY_TEXTURE_DISPLACEMENT(0), strValue) == aiReturn_SUCCESS) {
                TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), &defMem);
                normalTex = createShared<ImageNormal3DTexture>(image);
            }
            normalMaps.push_back(normalTex);
        }
        
        recursiveConstruct(objSrc, objSrc.mRootNode, materials, nodeOut);
    }
}

//namespace OBJ2 {
//    void construct(const File &file, SceneV2 &scene, ObjectRef &objects) {
//        std::vector<IColorTextureRef> colorTextures;
//        std::vector<FloatTextureRef> floatTextures;
//        std::vector<SurfaceMaterialRef> materials;
//        for (int mIdx = 0; mIdx < file.materials.size(); ++mIdx) {
//            const Material &matRef = file.materials[mIdx];
//            
//            if (matRef.Ni == 1.0f) {
//                IColorTextureRef diffuseTex;
//                if (matRef.map_Kd != "")
//                    diffuseTex = createShared<ImageTexture<Color>>(scene.getTextureMem(), matRef.map_Kd);
//                else
//                    diffuseTex = createShared<ConstantTexture<Color>>(scene.getTextureMem(), Color(matRef.Kd.R, matRef.Kd.G, matRef.Kd.B));
//                colorTextures.push_back(diffuseTex);
//                
//                FloatTextureRef sigmaTex = createShared<ConstantFloatTexture>(scene.getTextureMem(), 0.0f);
//                floatTextures.push_back(sigmaTex);
//                
//                materials.emplace_back(SurfaceMaterial::createMatte(scene.getMaterialMem(), nullptr, diffuseTex, sigmaTex));
//            }
//            else {
//                IColorTextureRef coeffRTex = createShared<ConstantTexture<Color>>(scene.getTextureMem(), Color(0.95f, 0.95f, 0.95f));
//                IColorTextureRef coeffTTex = createShared<ConstantTexture<Color>>(scene.getTextureMem(), Color(0.95f, 0.95f, 0.95f));
//                materials.emplace_back(SurfaceMaterial::createGlass(scene.getMaterialMem(), nullptr, coeffRTex, coeffTTex, matRef.Ni, 1.0f));
//            }
//        }
//        
//        ObjectGroupRef rootGroup = createShared<ObjectGroup>();
//        objects = rootGroup;
//        
//        for (int oIdx = 0; oIdx < file.objects.size(); ++oIdx) {
//            const Object &objSrc = file.objects[oIdx];
//            TriangleMeshRef obj = createShared<TriangleMesh>(scene);
//            rootGroup->addChild(obj);
//            
//            for (int p = 0; p < objSrc.positions.size(); ++p) {
//                const Position &pos = objSrc.positions[p];
//                obj->addPoint(::Point(pos.x, pos.y, pos.z));
//            }
//            for (int n = 0; n < objSrc.normals.size(); ++n) {
//                const Normal &norm = objSrc.normals[n];
//                obj->addNormal(::Normal(norm.x, norm.y, norm.z));
//            }
//            for (int t = 0; t < objSrc.texCoords.size(); ++t) {
//                const TexCoord &uv = objSrc.texCoords[t];
//                obj->addTexCoord(::TexCoord(uv.u, uv.v));
//            }
//            
//            for (int gIdx = 0; gIdx < objSrc.groups.size(); ++gIdx) {
//                const Group &groupSrc = objSrc.groups[gIdx];
//                
//                SurfaceMaterialRef mat = materials[groupSrc.matIdx];
//                obj->addMaterial(mat);
//                
//                SurfaceMaterial* matRef = mat.get();
//                
//                for (int fIdx = 0; fIdx < groupSrc.faces.size(); ++fIdx) {
//                    const Face &face = groupSrc.faces[fIdx];
//                    Triangle triangle{
//                        face.p0, face.p1, face.p2,
//                        face.n0, face.n1, face.n2,
//                        UINT32_MAX, UINT32_MAX, UINT32_MAX,
//                        face.t0, face.t1, face.t2,
//                        matRef
//                    };
//                    obj->addTriangle(triangle);
//                }
//            }
//        }
//    }
//}