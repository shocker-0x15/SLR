//
//  node_constructor.cpp
//
//  Created by 渡部 心 on 2015/03/24.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "node_constructor.h"

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <libSLR/MemoryAllocators/Allocator.h>
#include <libSLR/Core/transform.h>
#include <libSLR/SurfaceShape/TriangleSurfaceShape.h>
#include "textures.h"
#include "surface_materials.h"
#include "Scene/TriangleMeshNode.h"
#include "API.h"

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
    static void recursiveConstruct(const aiScene* objSrc, const aiNode* nodeSrc,
                                   const std::vector<SurfaceMaterialRef> &materials, const std::vector<NormalTextureRef> &normalMaps, const std::vector<FloatTextureRef> &alphaMaps,
                                   const MeshCallback &meshCallback, InternalNodeRef &nodeOut) {
        if (nodeSrc->mNumMeshes == 0 && nodeSrc->mNumChildren == 0) {
            nodeOut = nullptr;
            return;
        }
        
        const aiMatrix4x4 &tf = nodeSrc->mTransformation;
        float tfElems[] = {
            tf.a1, tf.a2, tf.a3, tf.a4,
            tf.b1, tf.b2, tf.b3, tf.b4,
            tf.c1, tf.c2, tf.c3, tf.c4,
            tf.d1, tf.d2, tf.d3, tf.d4,
        };
        
        nodeOut = createShared<InternalNode>(createShared<SLR::StaticTransform>(SLR::Matrix4x4(tfElems)));
        nodeOut->setName(nodeSrc->mName.C_Str());
        
        std::vector<Triangle> meshIndices;
        for (int m = 0; m < nodeSrc->mNumMeshes; ++m) {
            const aiMesh* mesh = objSrc->mMeshes[nodeSrc->mMeshes[m]];
            if (mesh->mPrimitiveTypes != aiPrimitiveType_TRIANGLE) {
                printf("ignored non triangle mesh.\n");
                continue;
            }
            
            TriangleMeshNodeRef surfMesh = createShared<TriangleMeshNode>();
            const SurfaceMaterialRef &surfMat = materials[mesh->mMaterialIndex];
            const NormalTextureRef &normalMap = normalMaps[mesh->mMaterialIndex];
            const FloatTextureRef &alphaMap = alphaMaps[mesh->mMaterialIndex];
            
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
            
            SLR::BoundingBox3D bbox;
            meshIndices.clear();
            for (int f = 0; f < mesh->mNumFaces; ++f) {
                const aiFace &face = mesh->mFaces[f];
                meshIndices.emplace_back(face.mIndices[0], face.mIndices[1], face.mIndices[2]);
                
                const aiVector3D &p = mesh->mVertices[face.mIndices[0]];
                bbox.unify(SLR::Point3D(p.x, p.y, p.z));
            }
            surfMesh->addMaterialGroup(surfMat, normalMap, alphaMap, std::move(meshIndices));
            
            MeshAttributeTuple meshAttr = meshCallback(surfMesh->getName(), surfMesh, bbox.minP, bbox.maxP);
            
            surfMesh->setName(mesh->mName.C_Str());

            surfMesh->useOnlyForBoundary(!meshAttr.render);
            surfMesh->setAxisForRadialTangent(meshAttr.axisForRadialTangent);
            
            nodeOut->addChildNode(surfMesh);
        }
        
        if (nodeSrc->mNumChildren) {
            for (int c = 0; c < nodeSrc->mNumChildren; ++c) {
                InternalNodeRef subNode;
                recursiveConstruct(objSrc, nodeSrc->mChildren[c], materials, normalMaps, alphaMaps, meshCallback, subNode);
                if (subNode != nullptr)
                    nodeOut->addChildNode(subNode);
            }
        }
    }
    
    SurfaceAttributeTuple createMaterialDefaultFunction(const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem) {
        using namespace SLR;
        aiReturn ret;
        (void)ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        const Texture2DMappingRef &mapping = Texture2DMapping::sharedInstanceRef();
        
        SpectrumTextureRef diffuseTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = Image::createTiledImage((pathPrefix + strValue.C_Str()).c_str(), mem, ImageStoreMode::AsIs, SpectrumType::Reflectance);
            diffuseTex = createShared<ImageSpectrumTexture>(mapping, image);
        }
        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
            AssetSpectrumRef sp = Spectrum::create(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, color[0], color[1], color[2]);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        else {
            AssetSpectrumRef sp = Spectrum::create(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        
        SurfaceMaterialRef mat = SurfaceMaterial::createMatte(diffuseTex, nullptr);
        
        NormalTextureRef normalTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DISPLACEMENT(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = Image::createTiledImage((pathPrefix + strValue.C_Str()).c_str(), mem, ImageStoreMode::NormalTexture, SpectrumType::Reflectance);
            normalTex = createShared<ImageNormalTexture>(mapping, image);
        }
        
        FloatTextureRef alphaTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_OPACITY(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = Image::createTiledImage((pathPrefix + strValue.C_Str()).c_str(), mem, ImageStoreMode::AlphaTexture, SpectrumType::Reflectance);
            alphaTex = createShared<ImageFloatTexture>(mapping, image);
        }
        
        return SurfaceAttributeTuple(mat, normalTex, alphaTex);
    }
    
#ifdef SLR_Platform_Windows_MSVC
    class FuncGetPathElement {
        std::string pathPrefix;
    public:
        FuncGetPathElement(const std::string &prefix) : pathPrefix(prefix) {}
        Element operator()(const aiString &str) const {
            return Element(pathPrefix + std::string(str.C_Str()));
        }
    };
#endif
    
    SurfaceAttributeTuple createMaterialFunction(const Function &userMatProc, ExecuteContext &context, ErrorMessage* err,
                                                 const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem) {
        using namespace SLR;
        aiString aiStr;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, aiStr);
        Element matName = Element::create<TypeMap::String>(std::string(aiStr.C_Str()));
        Element matAttrs = Element::create<TypeMap::Tuple>();
        
        auto &attrs = matAttrs.raw<TypeMap::Tuple>();
#ifdef SLR_Platform_Windows_MSVC
        FuncGetPathElement getPathElement{ pathPrefix };
#else
        auto getPathElement = [&pathPrefix](const aiString &str) {
            return Element::create<TypeMap::String>(pathPrefix + std::string(str.C_Str()));
        };
#endif
        auto getRGBElement = [](const float* RGB) {
            ParameterListRef values = createShared<ParameterList>();
            values->add("", Element(RGB[0]));
            values->add("", Element(RGB[1]));
            values->add("", Element(RGB[2]));
            return Element::createFromReference<TypeMap::Tuple>(values);
        };
        
        Element diffuseTexturesElement = Element::create<TypeMap::Tuple>();
        auto &diffuseTextures = diffuseTexturesElement.raw<TypeMap::Tuple>();
        for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_DIFFUSE); ++i) {
            if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(i), aiStr) == aiReturn_SUCCESS)
                diffuseTextures.add("", getPathElement(aiStr));
        }
        attrs.add("diffuse textures", diffuseTexturesElement);
        
        Element specularTexturesElement = Element::create<TypeMap::Tuple>();
        auto &specularTextures = specularTexturesElement.raw<TypeMap::Tuple>();
        for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_SPECULAR); ++i) {
            if (aiMat->Get(AI_MATKEY_TEXTURE_SPECULAR(i), aiStr) == aiReturn_SUCCESS)
                specularTextures.add("", getPathElement(aiStr));
        }
        attrs.add("specular textures", specularTexturesElement);
        
        Element emissiveTexturesElement = Element::create<TypeMap::Tuple>();
        auto &emissiveTextures = emissiveTexturesElement.raw<TypeMap::Tuple>();
        for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_EMISSIVE); ++i) {
            if (aiMat->Get(AI_MATKEY_TEXTURE_EMISSIVE(i), aiStr) == aiReturn_SUCCESS)
                emissiveTextures.add("", getPathElement(aiStr));
        }
        attrs.add("emissive textures", emissiveTexturesElement);
        
        Element heightTexturesElement = Element::create<TypeMap::Tuple>();
        auto &heightTextures = heightTexturesElement.raw<TypeMap::Tuple>();
        for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_HEIGHT); ++i) {
            if (aiMat->Get(AI_MATKEY_TEXTURE_HEIGHT(i), aiStr) == aiReturn_SUCCESS)
                heightTextures.add("", getPathElement(aiStr));
        }
        attrs.add("height textures", heightTexturesElement);
        
        Element normalTexturesElement = Element::create<TypeMap::Tuple>();
        auto &normalTextures = normalTexturesElement.raw<TypeMap::Tuple>();
        for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_NORMALS); ++i) {
            if (aiMat->Get(AI_MATKEY_TEXTURE_NORMALS(i), aiStr) == aiReturn_SUCCESS)
                normalTextures.add("", getPathElement(aiStr));
        }
        attrs.add("normal textures", normalTexturesElement);
        
        if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS)
            attrs.add("diffuse color", getRGBElement(color));
        
        if (aiMat->Get(AI_MATKEY_COLOR_SPECULAR, color, nullptr) == aiReturn_SUCCESS)
            attrs.add("specular color", getRGBElement(color));
        
        if (aiMat->Get(AI_MATKEY_COLOR_EMISSIVE, color, nullptr) == aiReturn_SUCCESS)
            attrs.add("emissive color", getRGBElement(color));
        
        ParameterList params;
        params.add("", matName);
        params.add("", matAttrs);
        
        // perform the user-defined material create function.
        Element result = userMatProc(params, context, err);
        if (err->error) {
            
        }
        
        if (result.type == Type::Tuple) {
            const ParameterList &tuple = result.raw<TypeMap::Tuple>();
            const Element &eMat = tuple(0);
            const Element &eNormalMap = tuple(1);
            const Element &eAlphaMap = tuple(2);
            if (eMat.type == Type::SurfaceMaterial) {
                SurfaceMaterialRef mat = eMat.rawRef<TypeMap::SurfaceMaterial>();
                NormalTextureRef normalMap;
                if (eNormalMap.type == Type::NormalTexture)
                    normalMap = eNormalMap.rawRef<TypeMap::NormalTexture>();
                FloatTextureRef alphaMap;
                if (eAlphaMap.type == Type::FloatTexture)
                    alphaMap = eAlphaMap.rawRef<TypeMap::FloatTexture>();
                
                return SurfaceAttributeTuple(mat, normalMap, alphaMap);
            }
        }
        else if (result.type == Type::SurfaceMaterial) {
            return SurfaceAttributeTuple(result.rawRef<TypeMap::SurfaceMaterial>(), nullptr, nullptr);
        }
        
        printf("User defined material function is invalid, fall back to the default function.\n");
        return createMaterialDefaultFunction(aiMat, pathPrefix, mem);
    }
    
    MeshAttributeTuple meshCallbackDefaultFunction(const std::string &name, const TriangleMeshNodeRef &mesh, const SLR::Point3D &minP, const SLR::Point3D &maxP) {
        return MeshAttributeTuple(true, -1);
    }
    
    MeshAttributeTuple meshCallbackFunction(const Function &meshProc, ExecuteContext &context, ErrorMessage* err,
                                            const std::string &name, const TriangleMeshNodeRef &mesh, const SLR::Point3D &minP, const SLR::Point3D &maxP) {
        Element elName = Element::create<TypeMap::String>(name);
        Element elMesh = Element::createFromReference<TypeMap::SurfaceNode>(mesh);
        Element elMinP = Element::create<TypeMap::Point>(minP);
        Element elMaxP = Element::create<TypeMap::Point>(maxP);
        
        ParameterList params;
        params.add("", elName);
        params.add("", elMesh);
        params.add("", elMinP);
        params.add("", elMaxP);
        
        // perform the user-defined mesh callback function.
        Element result = meshProc(params, context, err);
        if (err->error) {
            
        }
        
        if (result.type == Type::Bool) {
            return MeshAttributeTuple(result.asRaw<TypeMap::Bool>(), -1);
        }
        else if (result.type == Type::Tuple) {
            const ParameterList &returnedValues = result.asRaw<TypeMap::Tuple>();
            const Element &elRender = returnedValues(0);
            const Element &elAxisForRadialTangent = returnedValues(1);
            if (elRender.type == Type::Bool) {
                bool render = elRender.raw<TypeMap::Bool>();
                int32_t axisForRadialTangent = -1;
                if (elAxisForRadialTangent.type == Type::Integer) {
                    axisForRadialTangent = elAxisForRadialTangent.raw<TypeMap::Integer>();
                    if (axisForRadialTangent < -1 || axisForRadialTangent > 2)
                        axisForRadialTangent = -1;
                }
                
                return MeshAttributeTuple(render, axisForRadialTangent);
            }
        }
        
        printf("User defined mesh callback function is invalid, fall back to the default function.\n");
        return meshCallbackDefaultFunction(name, mesh, minP, maxP);
    }
    
    SLR_SCENEGRAPH_API void construct(const std::string &filePath, InternalNodeRef &nodeOut,
                                      const CreateMaterialFunction &materialFunc, const MeshCallback &meshCallback) {
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
        
        // create materials
        std::vector<SurfaceMaterialRef> materials;
        std::vector<NormalTextureRef> normalMaps;
        std::vector<FloatTextureRef> alphaMaps;
        for (int m = 0; m < scene->mNumMaterials; ++m) {
            const aiMaterial* aiMat = scene->mMaterials[m];
            
            SurfaceAttributeTuple surfAttr = materialFunc(aiMat, pathPrefix, &defMem);
            materials.push_back(surfAttr.material);
            normalMaps.push_back(surfAttr.normalMap);
            alphaMaps.push_back(surfAttr.alphaMap);
        }
        
        recursiveConstruct(scene, scene->mRootNode, materials, normalMaps, alphaMaps, meshCallback, nodeOut);
        
        nodeOut->setName(filePath);
        
        printf("Constructing: %s done.\n", filePath.c_str());
    }
}
