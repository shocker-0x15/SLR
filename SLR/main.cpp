//
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "defines.h"
#include "references.h"
#include <fstream>

#include "Helper/StopWatch.h"

#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include "SceneGraph/node_constructor.h"
#include "SceneGraph/Scene.h"
#include "SceneGraph/nodes.h"
#include "SceneGraph/camera_nodes.h"
#include "SceneGraph/TriangleMeshNode.h"
#include "SceneGraph/InfiniteSphereNode.h"

#include "Core/geometry.h"
#include "Core/SurfaceObject.h"
#include "Core/surface_material.h"
#include "Core/XORShift.h"
#include "Core/Image.h"
#include "Core/RenderSettings.h"
#include "Memory/ArenaAllocator.h"
#include "Textures/constant_textures.h"
#include "Textures/image_textures.h"
#include "Textures/checker_board_textures.h"
#include "Renderers/PathTracingRenderer.h"
#include "Renderers/AMCMCPPMRenderer.h"

Scene scene;

// RTC3 提出シーン
void buildScene(StopWatch &sw) {
    const auto removeModelFileName = [](const std::string &filePath) {
        return filePath.substr(0, filePath.find_last_of("/") + 1);
    };
    const auto addChild = [&removeModelFileName](InternalNodeRef &parent, const std::string &filePath, InternalNodeRef* readObj = nullptr,
                                                 Assimp::createMaterialFunction matFunc = Assimp::createMaterialDefaultFunction) {
        Assimp::Importer importer;
        const aiScene* sceneFile = importer.ReadFile(filePath, 0);
        printf("Reading: %s done.\n", filePath.c_str());
        
        InternalNodeRef objNode;
        Assimp::construct(*sceneFile, removeModelFileName(filePath), objNode, matFunc);
        objNode->setName(filePath);
        parent->addChildNode(objNode);
        
        if (readObj)
            *readObj = objNode;
    };
    
    InternalNodeRef &root = scene.rootNode();
    
    TiledImage2DRef IBLImage = TiledImage2D::create("images/Malibu_Overlook_3k.exr", &DefaultAllocator::instance(), SpectrumType::Illuminant);
    //    TiledImage2DRef IBLImage = TiledImage2D::create("images/Playa_Sunrise.exr", &DefaultAllocator::instance());
    //    TiledImage2DRef IBLImage = TiledImage2D::create("images/Etnies_Park_Center_3k.exr", &DefaultAllocator::instance());
    SpectrumTextureRef IBLTex = createShared<ImageSpectrumTexture>(IBLImage);
    InfiniteSphereNodeRef envNode = createShared<InfiniteSphereNode>(&scene, IBLTex, 10.0f);
    scene.setEnvNode(envNode);
    
    addChild(root, "models/plain_proto.assbin", nullptr);
    
    {
        ArenaAllocator mem;
        InternalNodeRef tempRoot = createShared<InternalNode>();
        tempRoot->setTransform(createShared<StaticTransform>(StaticTransform()));
        addChild(tempRoot, "models/plain_proto.assbin", nullptr);
        RenderingData renderData;
        tempRoot->getRenderingData(mem, nullptr, &renderData);
        auto aggr = createUnique<SurfaceObjectAggregate>(renderData.surfObjs);
        
        BoundingBox3D plainBound = aggr->bounds();
        
        InternalNodeRef grassNode = createShared<InternalNode>();
        grassNode->setTransform(createShared<StaticTransform>());
        InternalNodeRef dummy = createShared<InternalNode>();
        auto grassMaterial = [](const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem) {
            aiReturn ret;
            aiString strValue;
            float color[3];
            
            aiMat->Get(AI_MATKEY_NAME, strValue);
            
            SpectrumTextureRef diffuseTex;
            if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
                TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), mem, SpectrumType::Reflectance);
                diffuseTex = createShared<ImageSpectrumTexture>(image);
            }
            else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
                InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 2.0f * color[0], 2.0f * color[1], 2.0f * color[2]);
                diffuseTex = createShared<ConstantSpectrumTexture>(sp);
            }
            else {
                InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
                diffuseTex = createShared<ConstantSpectrumTexture>(sp);
            }
            
            return SurfaceMaterial::createMatte(diffuseTex, nullptr);
        };
        addChild(dummy, "models/grass.assbin", &grassNode, grassMaterial);
        ReferenceNodeRef grassRefNode = createShared<ReferenceNode>(grassNode);
        const float grassWidth = 1.162f;
        const float grassDepth = 0.752f;
        const float grassHeight = 0.180f;
        const uint32_t numX = 60;
        const uint32_t numY = 60;
        XORShift rng(41957139);
        for (int i = 0; i < numX; ++i) {
            for (int j = 0; j < numY; ++j) {
                InternalNodeRef translateNode = createShared<InternalNode>();
                Ray ray(Point3D(plainBound.minP.x + (plainBound.maxP.x - plainBound.minP.x) * (i + 0.5f) / numX,
                                100.0f,
                                plainBound.minP.z + (plainBound.maxP.z - plainBound.minP.z) * (j + 0.5f) / numY),
                        Vector3D(0, -1, 0), 0.0f);
                Intersection isect;
                aggr->intersect(ray, &isect);
                SurfacePoint surfPt;
                isect.getSurfacePoint(&surfPt);
                
                Matrix4x4 t = translate(Vector3D(surfPt.p));
                Vector3D axis = cross(Vector3D(0, 1, 0), surfPt.shadingFrame.z);
                float angle = std::acos(std::clamp(dot(Vector3D(0, 1, 0), surfPt.shadingFrame.z), -1.0f, 1.0f));
                if (angle < 0.0001f)
                    axis = Vector3D(1, 0, 0);
                Matrix4x4 r = rotate(angle, axis);
                Matrix4x4 s = scale(0.25f);
                Matrix4x4 ry = rotateY(float(2 * M_PI * rng.getFloat0cTo1o()));
                translateNode->setTransform(createShared<StaticTransform>(t * r * s * ry));
                translateNode->addChildNode(grassRefNode);
                root->addChildNode(translateNode);
            }
        }
    }
    
    InternalNodeRef lowpolyTreeNode;
    Assimp::createMaterialFunction treeMaterial = [](const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem) {
        aiReturn ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        SpectrumTextureRef diffuseTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), mem, SpectrumType::Reflectance);
            diffuseTex = createShared<ImageSpectrumTexture>(image);
        }
        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
            InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 2.0f * color[0], 2.0f * color[1], 2.0f * color[2]);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        else {
            InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        
        //        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
        InputSpectrumRef RsCol = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB, 0.025f, 0.025f, 0.025f);
        SpectrumTextureRef RsTex = createShared<ConstantSpectrumTexture>(RsCol);
        FloatTextureRef shininessTex = createShared<ConstantFloatTexture>(100.0f);
        return SurfaceMaterial::createAshikhminShirley(diffuseTex, RsTex, shininessTex, shininessTex);
    };
    addChild(root, "models/lowpoly_tree.assbin", &lowpolyTreeNode, treeMaterial);
    lowpolyTreeNode->setTransform(createShared<StaticTransform>(translate(0.3727f, 0.8596f, 0.8495f) * scale(0.25f)));
    
    InternalNodeRef toycarNode;
    Assimp::createMaterialFunction toycarMaterial = [](const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem) {
        aiReturn ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        SpectrumTextureRef diffuseTex;
        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
            TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), mem, SpectrumType::Reflectance);
            diffuseTex = createShared<ImageSpectrumTexture>(image);
        }
        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
            InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, color[0], color[1], color[2]);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        else {
            InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 1.0f, 0.0f, 1.0f);
            diffuseTex = createShared<ConstantSpectrumTexture>(sp);
        }
        
        //        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
        InputSpectrumRef RsCol = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB, 0.025f, 0.025f, 0.025f);
        SpectrumTextureRef RsTex = createShared<ConstantSpectrumTexture>(RsCol);
        FloatTextureRef shininessTex = createShared<ConstantFloatTexture>(100.0f);
        return SurfaceMaterial::createAshikhminShirley(diffuseTex, RsTex, shininessTex, shininessTex);
    };
    addChild(root, "models/toycar.assbin", &toycarNode, toycarMaterial);
    auto toycarTFStart = createUnique<StaticTransform>(translate(0.0f, 0.72593f, 1.4f) *
                                                       rotateZ(0.0f) * rotateY(float(105.319f * M_PI / 180.0f)) * rotateX(-float(16.663f * M_PI / 180.0f)) *
                                                       scale(0.1f));
    auto toycarTFEnd = createUnique<StaticTransform>(toycarTFStart->getMatrix4x4() * translate(0.0f, 0.0f, 0.25f));
    toycarNode->setTransform(createShared<StaticTransform>(toycarTFStart->getMatrix4x4()));
    //    toycarNode->setTransform(createShared<AnimatedTransform>(toycarTFStart.get(), toycarTFEnd.get(), 0.0f, 1.0f));
    
    InternalNodeRef cameraNode = createShared<InternalNode>();
    cameraNode->setName("Camera");
    CameraNodeRef camera = createShared<PerspectiveCameraNode>(16.0f / 9.0f, 30 * (M_PI / 180.0f), 0.0025f, 1.0f, 2.4f);
    cameraNode->addChildNode(camera);
    cameraNode->setTransform(createShared<StaticTransform>(translate(2.158f, 0.825f, 2.38f) * rotateY(-float((180.0f - 49.9f) * M_PI / 180)) * rotateX(float(-5.44f * M_PI / 180))));
    //    cameraNode->setTransform(createShared<StaticTransform>(
    //                             invert(lookAt(Point3D(0.0f, 12.0f, 12.0f), Point3D(0.0f, 0.0f, 0.0f), Vector3D(0, 1, 0))) * rotateY(float(M_PI))
    //                             ));
    root->addChildNode(cameraNode);
    
    root->setTransform(createShared<StaticTransform>(rotateY(-float(M_PI / 2))));
    
    scene.build();
}




int main(int argc, const char * argv[]) {
    using namespace std::chrono;
    initSpectrum();
    
    StopWatch stopwatch;
    StopWatchHiRes stopwatchHiRes;
    
    std::time_t ctimeLaunch = system_clock::to_time_t(stopwatch.start());
    printf("%s\n", std::ctime(&ctimeLaunch));
    
    stopwatch.start();
    buildScene(stopwatch);
    printf("build: %g [s]\n", stopwatch.stop() * 1e-3f);
    
    RenderSettings settings;
    settings.addItem(RenderSettingItem::ImageWidth, 1920);
    settings.addItem(RenderSettingItem::ImageHeight, 1080);
    settings.addItem(RenderSettingItem::TimeStart, 0.0f);
    settings.addItem(RenderSettingItem::TimeEnd, 1.0f);
    settings.addItem(RenderSettingItem::NumSamples, 65536);
    settings.addItem(RenderSettingItem::RNGSeed, 1509761209);
//    settings.addItem(RenderSettingItem::SensorResponse, 3.0f);
    settings.addItem(RenderSettingItem::SensorResponse, float(1.0f / (M_PI * 0.0025f * 0.0025f)));
    PathTracingRenderer renderer;
    renderer.render(scene, settings);
//    AMCMCPPMRenderer renderer;
//    renderer.render(scene, settings);
    
    return 0;
}
