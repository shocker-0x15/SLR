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

// Cornell Box (box-box)
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
    
    InternalNodeRef CBNode;
    addChild(root, "models/Cornell_box_RB.assbin", &CBNode);
    {
        InternalNodeRef lightNode = createShared<InternalNode>();
        lightNode->setName("AreaLight");
        CBNode->addChildNode(lightNode);
        
        lightNode->setTransform(createShared<StaticTransform>(translate(0.0f, 0.999f, 0.0f)));
        {
            TriangleMeshNodeRef lightMesh = createShared<TriangleMeshNode>();
            
            InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, 0.9f, 0.9f, 0.9f);
            SpectrumTextureRef diffuseTex = createShared<ConstantSpectrumTexture>(sp);
            SurfaceMaterialRef scatterMat = SurfaceMaterial::createMatte(diffuseTex, nullptr);
            
            InputSpectrumRef lightColor = createInputSpectrum(SpectrumType::Illuminant, StandardIlluminant::MinWavelength, StandardIlluminant::MaxWavelength,
                                                              StandardIlluminant::D65, StandardIlluminant::NumSamples);
//            InputSpectrumRef lightColor = createInputSpectrum(ColorSpace::sRGB, 500.0f, 500.0f, 500.0f);
            SpectrumTextureRef difLightTex = createShared<ConstantSpectrumTexture>(lightColor);
            EmitterSurfacePropertyRef emitMat = EmitterSurfaceProperty::createDiffuseEmitter(difLightTex);
            
            SurfaceMaterialRef surfMat = SurfaceMaterial::createEmitterSurfaceMaterial(scatterMat, emitMat);
            
            lightMesh->addVertex(Vertex(Point3D(-0.25f, 0, -0.25f), Normal3D(0, -1, 0), Tangent3D(1, 0, 0), TexCoord2D(0, 0)));
            lightMesh->addVertex(Vertex(Point3D( 0.25f, 0, -0.25f), Normal3D(0, -1, 0), Tangent3D(1, 0, 0), TexCoord2D(1, 0)));
            lightMesh->addVertex(Vertex(Point3D( 0.25f, 0,  0.25f), Normal3D(0, -1, 0), Tangent3D(1, 0, 0), TexCoord2D(1, 1)));
            lightMesh->addVertex(Vertex(Point3D(-0.25f, 0,  0.25f), Normal3D(0, -1, 0), Tangent3D(1, 0, 0), TexCoord2D(0, 1)));
            
            lightMesh->addTriangle(0, 1, 2, surfMat);
            lightMesh->addTriangle(0, 2, 3, surfMat);
            
            lightNode->addChildNode(lightMesh);
        }
        
        {
            TriangleMeshNodeRef mesh = createShared<TriangleMeshNode>();
            
            std::vector<SurfaceMaterialRef> mats;
            for (int i = 0; i < 24; ++i) {
                InputSpectrumRef sp = createInputSpectrum(SpectrumType::Reflectance, ColorChecker::MinWavelength, ColorChecker::MaxWavelength,
                                                          ColorChecker::Spectra[i], ColorChecker::NumSamples);
                SpectrumTextureRef diffuseTex = createShared<ConstantSpectrumTexture>(sp);
                SurfaceMaterialRef scatterMat = SurfaceMaterial::createMatte(diffuseTex, nullptr);
                mats.push_back(scatterMat);
            }
            
            for (int i = 0; i < 24; ++i) {
                Matrix4x4 transform = rotateY(float(M_PI / 2)) * translate(0.0f, 0.0f, -0.999f) * scale(0.9f / 3.0f) * translate(-3.0f + float(i % 6), 1.0f - float(i / 6), 0.0f);
                mesh->addVertex(Vertex(transform * Point3D(0.0f, 0.0f, 0.0f), transform * Normal3D(0, -1, 0), transform * Tangent3D(1, 0, 0), TexCoord2D(0, 0)));
                mesh->addVertex(Vertex(transform * Point3D(1.0f, 0.0f, 0.0f), transform * Normal3D(0, -1, 0), transform * Tangent3D(1, 0, 0), TexCoord2D(1, 0)));
                mesh->addVertex(Vertex(transform * Point3D(1.0f, 1.0f, 0.0f), transform * Normal3D(0, -1, 0), transform * Tangent3D(1, 0, 0), TexCoord2D(1, 1)));
                mesh->addVertex(Vertex(transform * Point3D(0.0f, 1.0f, 0.0f), transform * Normal3D(0, -1, 0), transform * Tangent3D(1, 0, 0), TexCoord2D(0, 1)));
                
                uint32_t idxBase = 4 * i;
                mesh->addTriangle(idxBase + 0, idxBase + 1, idxBase + 2, mats[i]);
                mesh->addTriangle(idxBase + 0, idxBase + 2, idxBase + 3, mats[i]);
            }
            
            CBNode->addChildNode(mesh);
        }
    }
    CBNode->setTransform(createShared<StaticTransform>(rotateY((float)(-M_PI / 2))));
    
    auto boxMaterial = [](const aiMaterial* aiMat, const std::string &pathPrefix, Allocator* mem) {
        aiReturn ret;
        aiString strValue;
        float color[3];
        
        aiMat->Get(AI_MATKEY_NAME, strValue);
        
        //        SpectrumTextureRef diffuseTex;
        //        if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(0), strValue) == aiReturn_SUCCESS) {
        //            TiledImage2DRef image = TiledImage2D::create((pathPrefix + strValue.C_Str()).c_str(), mem);
        //            diffuseTex = createShared<ImageSpectrumTexture>(image);
        //        }
        //        else if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color, nullptr) == aiReturn_SUCCESS) {
        //            diffuseTex = createShared<ConstantSpectrumTexture>(inverseGammaCorrection(2.0f * Spectrum(color[0], color[1], color[2])));
        //        }
        //        else {
        //            diffuseTex = createShared<ConstantSpectrumTexture>(Spectrum(1.0f, 0.0f, 1.0f));
        //        }
        //
        //        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
        
//        InputSpectrumRef sp = createShared<UpsampledContinuousSpectrum>(ColorSpace::sRGB_NonLinear, 0.5f, 0.5f, 0.5f);
//        SpectrumTextureRef diffuseTex = createShared<ConstantSpectrumTexture>(sp);
//        return SurfaceMaterial::createMatte(diffuseTex, nullptr);
        
        float coeffRSrc[] = {1.0f, 1.0f};
        InputSpectrumRef coeffR = createInputSpectrum(SpectrumType::Reflectance, 360.0f, 830.0f, coeffRSrc, 2);
        SpectrumTextureRef texCoeffR = createShared<ConstantSpectrumTexture>(coeffR);
        SpectrumTextureRef texCoeffT = createShared<ConstantSpectrumTexture>(coeffR);
        using namespace SpectrumLibrary;
        float lambdas[] = {360.0f, 830.0f};
        float etas[] = {1.0f, 2.4f};
        InputSpectrumRef etaExt = createInputSpectrum(SpectrumType::IndexOfRefraction, Air.minLambdas, Air.maxLambdas, Air.etas, Air.numSamples);
        InputSpectrumRef etaInt = createInputSpectrum(SpectrumType::IndexOfRefraction, lambdas, etas, 2);
        SpectrumTextureRef texEtaExt = createShared<ConstantSpectrumTexture>(etaExt);
        SpectrumTextureRef texEtaInt = createShared<ConstantSpectrumTexture>(etaInt);
        SpatialFresnelDielectricRef fresnel = createShared<SpatialFresnelDielectric>(texEtaExt, texEtaInt);
        return SurfaceMaterial::createGlass(texCoeffR, texCoeffT, fresnel);
    };
    
//    InternalNodeRef leftBoxNode;
//    addChild(CBNode, "models/box.assbin", &leftBoxNode, boxMaterial);
//    leftBoxNode->setTransform(createShared<StaticTransform>(rotateY(float(M_PI / 2)) *
//                                                            translate(-0.4f, -0.999f, -0.4f) * rotateY(float(M_PI / 12)) * scale(0.25f, 0.25f, 0.25f) *
//                                                            translate(0.0f, 1.0f, 0.0f)
//                                                            ));
//    InternalNodeRef rightBoxNode;
//    addChild(CBNode, "models/box.assbin", &rightBoxNode, boxMaterial);
//    rightBoxNode->setTransform(createShared<StaticTransform>(rotateY(float(M_PI / 2)) *
//                                                             translate(0.4f, -0.999f, 0.4f) * rotateY(-float(M_PI / 12)) * scale(0.25f) *
//                                                             translate(0.0f, 1.0f, 0.0f)
//                                                             ));
    InternalNodeRef pikachuNode;
    addChild(CBNode, "models/Pikachu_corrected.assbin", &pikachuNode, boxMaterial);
    pikachuNode->setTransform(createShared<StaticTransform>(translate(0.0f, -1.0f, 0.0f) * rotateY(float(135 * M_PI / 180.0f)) * scale(0.2f)
                                                            ));
    
    InternalNodeRef cameraNode = createShared<InternalNode>();
    cameraNode->setName("Camera");
    CameraNodeRef camera = createShared<PerspectiveCameraNode>(1.0f, 30 * (M_PI / 180.0f), 0.0f, 1.0f, 6.0f);
    cameraNode->addChildNode(camera);
    cameraNode->setTransform(createShared<StaticTransform>(translate(0.0f, 0.0f, 5.0f) * rotateY(-float(M_PI))));
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
    settings.addItem(RenderSettingItem::ImageWidth, 1024);
    settings.addItem(RenderSettingItem::ImageHeight, 1024);
    settings.addItem(RenderSettingItem::TimeStart, 0.0f);
    settings.addItem(RenderSettingItem::TimeEnd, 1.0f);
    settings.addItem(RenderSettingItem::NumSamples, 65536);
    settings.addItem(RenderSettingItem::RNGSeed, 1509761209);
//    settings.addItem(RenderSettingItem::SensorResponse, 3.0f);
    settings.addItem(RenderSettingItem::SensorResponse, 10.0f);
    PathTracingRenderer renderer;
    renderer.render(scene, settings);
//    AMCMCPPMRenderer renderer;
//    renderer.render(scene, settings);
    
    return 0;
}
