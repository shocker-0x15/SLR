
//
//  API.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#include "API.hpp"

#include "Parser/SceneParsingDriver.h"
#include <libSLR/Core/Transform.h>

#include "TriangleMeshNode.h"
#include "camera_nodes.h"

#include "nodes.h"
#include "node_constructor.h"

#include <libSLR/Renderers/PathTracingRenderer.h>

#include <libSLR/BasicTypes/Spectrum.h>
#include <libSLR/Core/Image.h>
#include "image_loader.h"
#include "textures.hpp"
#include "surface_materials.hpp"

namespace SLRSceneGraph {
    bool readScene(const std::string &filePath, Scene* scene, RenderingContext* context) {
        TypeInfo::init();
        SceneParsingDriver parser;
//        parser.traceParsing = true;
        StatementsRef statements = parser.parse(filePath);
        SLRAssert(statements, "Failed to parse scene file: %s", filePath.c_str());
        
        ExecuteContext executeContext;
        ErrorMessage errMsg;
        executeContext.scene = scene;
        executeContext.renderingContext = context;
        executeContext.stackVariables["root"] = Element(TypeMap::Node(), scene->rootNode());
        for (int i = 0; i < statements->size(); ++i) {
            StatementRef statement = statements->at(i);
            if (!statement->perform(executeContext, &errMsg)) {
                printf("%s\n", errMsg.message.c_str());
                return false;
            }
        }
        return true;
    }
    
    class Function {
        struct ArgInfo {
            std::string name;
            Type expectedType;
            Element defaultValue;
        };
        typedef std::function<Element(const std::map<std::string, Element> &, ErrorMessage*)> A;
        typedef std::function<Element(const std::map<std::string, Element> &, RenderingContext*, ErrorMessage*)> B;
        
        const std::vector<ArgInfo> signature;
        const A procA;
        const B procB;
        
        static Element NoOpA(const std::map<std::string, Element> &args, ErrorMessage* err) { SLRAssert(false, "Not implemented."); return Element(); };
        static Element NoOpB(const std::map<std::string, Element> &args, RenderingContext* context, ErrorMessage* err) { SLRAssert(false, "Not implemented."); return Element(); };
    public:
        Function(const std::vector<ArgInfo> &sig) : signature(sig), procA(NoOpA), procB(NoOpB) { };
        Function(const std::vector<ArgInfo> &sig, const A &proc) : signature(sig), procA(proc), procB(NoOpB) { };
        Function(const std::vector<ArgInfo> &sig, const B &proc) : signature(sig), procA(NoOpA), procB(proc) { };
        
        bool mapParamsToArgs(const ParameterList &params, std::map<std::string, Element>* args) const {
            size_t numArgs = signature.size();
            std::vector<bool> assigned(numArgs, false);
            for (auto namedParam : params.named) {
                const std::string key = namedParam.first;
                const Element value = namedParam.second;
                size_t idx = std::distance(std::begin(signature),
                                           std::find_if(std::begin(signature), std::end(signature),
                                                        [&key, &value](const ArgInfo &arg) {
                                                            return arg.name == key && value.isConvertibleTo(arg.expectedType);
                                                        }));
                if (idx == numArgs) {
                    args->clear();
                    return false;
                }
                const ArgInfo &argInfo = signature[idx];
                (*args)[argInfo.name] = value.convertTo(argInfo.expectedType);
                assigned[idx] = true;
            }
            for (auto it = std::begin(params.unnamed); it != std::end(params.unnamed); ++it) {
                const Element &value = *it;
                size_t idx = std::distance(std::begin(signature),
                                           std::find_if(std::begin(signature), std::end(signature),
                                                        [this, &value, &assigned](const ArgInfo &arg) {
                                                            size_t curIdx = std::distance(&signature[0], &arg);
                                                            return assigned[curIdx] == false && value.isConvertibleTo(arg.expectedType);
                                                        }));
                if (idx == numArgs) {
                    args->clear();
                    return false;
                }
                const ArgInfo &argInfo = signature[idx];
                (*args)[argInfo.name] = value.convertTo(argInfo.expectedType);
                assigned[idx] = true;
            }
            for (int i = 0; i < numArgs; ++i) {
                if (assigned[i])
                    continue;
                const ArgInfo &argInfo = signature[i];
                if (argInfo.defaultValue.type == Type::Void) {
                    args->clear();
                    return false;
                }
                (*args)[argInfo.name] = argInfo.defaultValue;
            }
            
            return true;
        };
        
        Element operator()(const ParameterList &params, ErrorMessage* err) const {
            std::map<std::string, Element> args;
            if (mapParamsToArgs(params, &args)) {
                *err = ErrorMessage();
                return procA(args, err);
            }
            *err = ErrorMessage("Parameters are invalid.");
            return Element();
        };
        
        Element operator()(const ParameterList &params, RenderingContext* context, ErrorMessage* err) const {
            std::map<std::string, Element> args;
            if (mapParamsToArgs(params, &args)) {
                *err = ErrorMessage();
                return procB(args, context, err);
            }
            *err = ErrorMessage("Parameters are invalid.");
            return Element();
        };
        
        template <typename T>
        T perform(const std::function<T(const std::map<std::string, Element> &)> &proc, const ParameterList &params, const T &failValue = T(), ErrorMessage* err = nullptr) const {
            std::map<std::string, Element> args;
            if (mapParamsToArgs(params, &args)) {
                if (err)
                    *err = ErrorMessage();
                return proc(args);
            }
            if (err)
                *err = ErrorMessage("Parameters are invalid.");
            return failValue;
        };
    };
    
    // tx, ty, tz
    Element Translate(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float tx = args.at("x").as<double>();
                float ty = args.at("y").as<double>();
                float tz = args.at("z").as<double>();
                
                return Element(TypeMap::Matrix(), SLR::translate(tx, ty, tz));
            }
        };
        return proc(params, err);
    }
    // angle
    Element RotateX(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"angle", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float angle = args.at("angle").as<double>();
                return Element(TypeMap::Matrix(), SLR::rotateX(angle));
            }
        };
        return proc(params, err);
    }
    // angle
    Element RotateY(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"angle", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float angle = args.at("angle").as<double>();
                return Element(TypeMap::Matrix(), SLR::rotateY(angle));
            }
        };
        return proc(params, err);
    }
    // angle
    Element RotateZ(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"angle", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float angle = args.at("angle").as<double>();
                return Element(TypeMap::Matrix(), SLR::rotateZ(angle));
            }
        };
        return proc(params, err);
    }
    // sx, sy, sz
    // s
    Element Scale(const ParameterList &params, ErrorMessage* err) {
        static const Function proc0{
            {{"s", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float s = args.at("s").as<double>();
                return Element(TypeMap::Matrix(), SLR::scale(s));
            }
        };
        static const Function proc1{
            {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float sx = args.at("x").as<double>();
                float sy = args.at("y").as<double>();
                float sz = args.at("z").as<double>();
                return Element(TypeMap::Matrix(), SLR::scale(sx, sy, sz));
            }
        };
        static const auto procs = std::make_array<Function>(proc0, proc1);
        for (int i = 0; i < procs.size(); ++i) {
            Element elem = procs[i](params, err);
            if (!err->error)
                return elem;
        }
        return Element();
    }
    
    static bool strToSpectrumType(const std::string &str, SLR::SpectrumType* type) {
        if (str == "Reflectance")
            *type = SLR::SpectrumType::Reflectance;
        else if (str == "Illuminant")
            *type = SLR::SpectrumType::Illuminant;
        else if (str == "IndexOfRefraction")
            *type = SLR::SpectrumType::IndexOfRefraction;
        else
            return false;
        return true;
    }
    
    static bool strToColorSpace(const std::string &str, SLR::ColorSpace* space) {
        if (str == "Rec709")
            *space = SLR::ColorSpace::sRGB;
        else if (str == "sRGB")
            *space = SLR::ColorSpace::sRGB_NonLinear;
        else if (str == "xyY")
            *space = SLR::ColorSpace::xyY;
        else if (str == "XYZ")
            *space = SLR::ColorSpace::XYZ;
        else
            return false;
        return true;
    }
    
    // type = Reflectance, e
    // type = Reflectance, space = sRGB, e0, e1, e2
    // minWL, maxWL, values
    // wls, values
    // ID
    Element CreateSpectrum(const ParameterList &params, ErrorMessage* err) {
        static const Function proc0{
            {
                {"type", Type::String},
                {"value", Type::RealNumber}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                std::string typeStr = args.at("type").as<std::string>();
                SLR::SpectrumType type;
                if (!strToSpectrumType(typeStr, &type)) {
                    *err = ErrorMessage("Specified spectrum type is invalid.");
                    return Element();
                }
                float value = args.at("value").as<double>();
                
                return Element(TypeMap::Spectrum(), Spectrum::create(type, SLR::ColorSpace::sRGB, value, value, value));
            }
        };
        static const Function proc1{
            {
                {"type", Type::String, Element(TypeMap::String(), "Reflectance")},
                {"space", Type::String, Element(TypeMap::String(), "sRGB")},
                {"e0", Type::RealNumber},
                {"e1", Type::RealNumber},
                {"e2", Type::RealNumber}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                std::string typeStr = args.at("type").as<std::string>();
                SLR::SpectrumType type;
                if (!strToSpectrumType(typeStr, &type)) {
                    *err = ErrorMessage("Specified spectrum type is invalid.");
                    return Element();
                }
                std::string spaceStr = args.at("space").as<std::string>();
                SLR::ColorSpace space;
                if (!strToColorSpace(spaceStr, &space)) {
                    *err = ErrorMessage("Specified color space is invalid.");
                    return Element();
                }
                float e0 = args.at("e0").as<double>();
                float e1 = args.at("e1").as<double>();
                float e2 = args.at("e2").as<double>();
                
                return Element(TypeMap::Spectrum(), Spectrum::create(type, space, e0, e1, e2));
            }
        };
        static const Function proc2{
            {
                {"minWL", Type::RealNumber},
                {"maxWL", Type::RealNumber},
                {"values", Type::RealNumber}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                return Element();
            }
        };
        static const Function proc3{
            {
                {"wls", Type::RealNumber},
                {"values", Type::RealNumber}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err)  {
                return Element();
            }
        };
        static const Function proc4{
            {
                {"ID", Type::String}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err)  {
                using namespace SLR;
                InputSpectrumRef spectrum;
                if (args.at("ID").as<std::string>() == "D65")
                    spectrum = Spectrum::create(SpectrumType::Illuminant, StandardIlluminant::MinWavelength, StandardIlluminant::MaxWavelength,
                                                StandardIlluminant::D65, StandardIlluminant::NumSamples);
                else
                    *err = ErrorMessage("unrecognized spectrum ID.");
                return Element(TypeMap::Spectrum(), spectrum);
            }
        };
        static const auto procs = std::make_array<Function>(proc0, proc1, proc2, proc3, proc4);
        for (int i = 0; i < procs.size(); ++i) {
            Element elem = procs[i](params, err);
            if (!err->error)
                return elem;
        }
        return Element();
    }
    
    Element CreateSpectrumTexture(const ParameterList &params, ErrorMessage* err) {
        static const Function proc0{
            {{"spectrum", Type::Spectrum}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                InputSpectrumRef spectrum = args.at("spectrum").asRef<SLR::InputSpectrum>();
                return Element(TypeMap::SpectrumTexture(), createShared<ConstantSpectrumTexture>(spectrum));
            }
        };
        static const auto procs = std::make_array<Function>(proc0);
        for (int i = 0; i < procs.size(); ++i) {
            Element elem = procs[i](params, err);
            if (!err->error)
                return elem;
        }
        return Element();
    }
    
    Element CreateMatte(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"reflectance", Type::SpectrumTexture}, {"sigma", Type::FloatTexture, Element(TypeMap::FloatTexture(), nullptr)}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                SpectrumTextureRef reflectance = args.at("reflectance").asRef<SpectrumTexture>();
                FloatTextureRef sigma = args.at("sigma").asRef<FloatTexture>();
                return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createMatte(reflectance, sigma));
            }
        };
        return proc(params, err);
    }
    
    Element CreateDiffuseEmitter(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"emittance", Type::SpectrumTexture}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                SpectrumTextureRef emittance = args.at("emittance").asRef<SpectrumTexture>();
                return Element(TypeMap::EmitterSurfaceProperty(), SurfaceMaterial::createDiffuseEmitter(emittance));
            }
        };
        return proc(params, err);
    }
    
    Element CreateEmitterSurfaceMaterial(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"mat", Type::SurfaceMaterial}, {"emit", Type::EmitterSurfaceProperty}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                SurfaceMaterialRef mat = args.at("mat").asRef<SurfaceMaterial>();
                EmitterSurfacePropertyRef emit = args.at("emit").asRef<EmitterSurfaceProperty>();
                return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createEmitterSurfaceMaterial(mat, emit));
            }
        };
        return proc(params, err);
    }
    
    // vertices, faces
    Element CreateMesh(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"vertices", Type::Tuple}, {"faces", Type::Tuple}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                const std::vector<Element> &vertices = args.at("vertices").as<ParameterList>().unnamed;
                const std::vector<Element> &faces = args.at("faces").as<ParameterList>().unnamed;
                
                TriangleMeshNodeRef lightMesh = createShared<TriangleMeshNode>();
                for (int i = 0; i < vertices.size(); ++i) {
                    static const Function sigVertex{
                        {{"position", Type::Tuple}, {"normal", Type::Tuple}, {"tangent", Type::Tuple}, {"texCoord", Type::Tuple}}
                    };
                    static const auto procVertex = [](const std::map<std::string, Element> &arg) {
                        static const Function sigPosition = Function({{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                        static const Function sigNormal = Function({{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                        static const Function sigTangent = Function({{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                        static const Function sigTexCoord = Function({{"u", Type::RealNumber}, {"v", Type::RealNumber}});
                        static const auto procPosition = [](const std::map<std::string, Element> &arg) {
                            return SLR::Point3D(arg.at("x").as<double>(), arg.at("y").as<double>(), arg.at("z").as<double>());
                        };
                        static const auto procNormal = [](const std::map<std::string, Element> &arg) {
                            return SLR::Normal3D(arg.at("x").as<double>(), arg.at("y").as<double>(), arg.at("z").as<double>());
                        };
                        static const auto procTangent = [](const std::map<std::string, Element> &arg) {
                            return SLR::Tangent3D(arg.at("x").as<double>(), arg.at("y").as<double>(), arg.at("z").as<double>());
                        };
                        static const auto procTexCoord = [](const std::map<std::string, Element> &arg) {
                            return SLR::TexCoord2D(arg.at("u").as<double>(), arg.at("v").as<double>());
                        };
                        return SLR::Vertex(sigPosition.perform<SLR::Point3D>(procPosition, arg.at("position").as<ParameterList>()),
                                           sigNormal.perform<SLR::Normal3D>(procNormal, arg.at("normal").as<ParameterList>()),
                                           sigTangent.perform<SLR::Tangent3D>(procTangent, arg.at("tangent").as<ParameterList>()),
                                           sigTexCoord.perform<SLR::TexCoord2D>(procTexCoord, arg.at("texCoord").as<ParameterList>()));
                    };
                    lightMesh->addVertex(sigVertex.perform<SLR::Vertex>(procVertex, vertices[i].as<ParameterList>()));
                }
                for (int i = 0; i < faces.size(); ++i) {
                    struct Face {
                        uint32_t v0, v1, v2;
                        SurfaceMaterialRef mat;
                    };
                    static const Function sigFace{
                        {{"indices", Type::Tuple}, {"material", Type::SurfaceMaterial}}
                    };
                    static const auto procFace = [](const std::map<std::string, Element> &arg) {
                        static const Function sigIndices = Function({{"v0", Type::Integer}, {"v1", Type::Integer}, {"v2", Type::Integer}});
                        std::map<std::string, Element> indices;
                        sigIndices.mapParamsToArgs(arg.at("indices").as<ParameterList>(), &indices);
                        uint32_t v0 = (uint32_t)indices.at("v0").as<int32_t>();
                        uint32_t v1 = (uint32_t)indices.at("v1").as<int32_t>();
                        uint32_t v2 = (uint32_t)indices.at("v2").as<int32_t>();
                        return Face{v0, v1, v2, arg.at("material").asRef<SurfaceMaterial>()};
                    };
                    Face f = sigFace.perform<Face>(procFace, faces[i].as<ParameterList>());
                    lightMesh->addTriangle(f.v0, f.v1, f.v2, f.mat);
                }
                
                return Element(TypeMap::Mesh(), lightMesh);
            }
        };
        return proc(params, err);
    }
    
    Element CreateNode(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                return Element(TypeMap::Node(), InternalNode());
            }
        };
        return proc(params, err);
    }
    
    Element SetTransform(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"node", Type::Node}, {"transform", Type::Transform}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                InternalNodeRef node = args.at("node").asRef<InternalNode>();
                TransformRef tf = args.at("transform").asRef<SLR::Transform>();
                node->setTransform(tf);
                return Element();
            }
        };
        return proc(params, err);
    }
    
    Element AddChild(const ParameterList &params, ErrorMessage* err) {
        static const Function proc0{
            {{"parent", Type::Node}, {"child", Type::Node}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                InternalNodeRef parent = args.at("parent").asRef<InternalNode>();
                InternalNodeRef child = args.at("child").asRef<InternalNode>();
                parent->addChildNode(child);
                return Element();
            }
        };
        static const Function proc1{
            {{"parent", Type::Node}, {"child", Type::Mesh}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                InternalNodeRef parent = args.at("parent").asRef<InternalNode>();
                TriangleMeshNodeRef child = args.at("child").asRef<TriangleMeshNode>();
                parent->addChildNode(child);
                return Element();
            }
        };
        static const Function proc2{
            {{"parent", Type::Node}, {"child", Type::Camera}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                InternalNodeRef parent = args.at("parent").asRef<InternalNode>();
                CameraNodeRef child = args.at("child").asRef<CameraNode>();
                parent->addChildNode(child);
                return Element();
            }
        };
        static const auto procs = std::make_array<Function>(proc0, proc1, proc2);
        for (int i = 0; i < procs.size(); ++i) {
            Element elem = procs[i](params, err);
            if (!err->error)
                return elem;
        }
        return Element();
    }
    
    Element Load3DModel(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {{"path", Type::String}},
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                std::string path = args.at("path").as<std::string>();
                
                InternalNodeRef modelNode;
                construct(path, modelNode);
                modelNode->setName(path);
                
                return Element(TypeMap::Node(), modelNode);
            }
        };
        return proc(params, err);
    }
    
    // width = 1024, height = 1024, aspect = 1.0, fovY = 0.5235987756, sensitivity = 0,
    // radius = 0.0, imgDist = 0.02, objDist = 5
    Element CreatePerspectiveCamera(const ParameterList &params, ErrorMessage* err) {
        static const Function proc{
            {
                {"sensitivity", Type::RealNumber, Element(0.0)},
                {"aspect", Type::RealNumber, Element(1.0)}, {"fovY", Type::RealNumber, Element(0.5235987756)},
                {"radius", Type::RealNumber, Element(0.0)}, {"imgDist", Type::RealNumber, Element(0.02)}, {"objDist", Type::RealNumber, Element(5.0)}
            },
            [](const std::map<std::string, Element> &args, ErrorMessage* err) {
                float sensitivity = args.at("sensitivity").as<double>();
                float aspect = args.at("aspect").as<double>();
                float fovY = args.at("fovY").as<double>();
                float radius = args.at("radius").as<double>();
                float imgDist = args.at("imgDist").as<double>();
                float objDist = args.at("objDist").as<double>();
                
                return Element(TypeMap::Camera(), createShared<PerspectiveCameraNode>(sensitivity, aspect, fovY, radius, imgDist, objDist));
            }
        };
        return proc(params, err);
    }
    
    Element SetRenderer(const ParameterList &params, RenderingContext* context, ErrorMessage* err) {
        Element paramMethod = params("method", Type::String);
        if (paramMethod.type == Type::Void) {
            *err = ErrorMessage("Rendering method is not specified.");
            return Element();
        }
        std::string method = paramMethod.as<std::string>();
        if (method == "PT") {
            static const Function proc{
                {{"method", Type::String}, {"samples", Type::Integer}},
                [](const std::map<std::string, Element> &args, RenderingContext* context, ErrorMessage* err) {
                    uint32_t spp = args.at("samples").as<uint32_t>();
                    context->renderer = createUnique<SLR::PathTracingRenderer>(spp);
                    return Element();
                }
            };
            return proc(params, context, err);
        }
        else {
            *err = ErrorMessage("Unknown method is specified.");
        }
        return Element();
    }
    
    Element SetRenderSettings(const ParameterList &params, RenderingContext* context, ErrorMessage* err) {
        static const Function proc{
            {
                {"width", Type::Integer, Element(1024)}, {"height", Type::Integer, Element(1024)},
                {"timeStart", Type::RealNumber, Element(0.0)}, {"timeEnd", Type::RealNumber, Element(0.0)},
                {"brightness", Type::RealNumber, Element(1.0f)}, 
                {"rngSeed", Type::Integer, Element(1509761209)},
            },
            [](const std::map<std::string, Element> &args, RenderingContext* context, ErrorMessage* err) {
                context->width = args.at("width").as<int32_t>();
                context->height = args.at("height").as<int32_t>();
                context->timeStart = args.at("timeStart").as<double>();
                context->timeEnd = args.at("timeEnd").as<double>();
                context->brightness = args.at("brightness").as<double>();
                context->rngSeed = args.at("rngSeed").as<int32_t>();
                
                return Element();
            }
        };
        return proc(params, context, err);
    }
    
    namespace Spectrum {
        using namespace SLR;
        
#ifdef Use_Spectral_Representation
        InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
            return createShared<UpsampledContinuousSpectrum>(spType, space, e0, e1, e2);
        }
        InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<RegularContinuousSpectrum>(minLambda, maxLambda, values, numSamples);
        }
        InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<IrregularContinuousSpectrum>(lambdas, values, numSamples);
        }
#else
        static void spectrum_to_XYZ(SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples, SpectrumFloat XYZ[3]) {
            const SpectrumFloat CMFBinWidth = (WavelengthHighBound - WavelengthLowBound) / (NumCMFSamples - 1);
            const SpectrumFloat binWidth = (maxLambda - minLambda) / (numSamples - 1);
            uint32_t curCMFIdx = 0;
            uint32_t baseIdx = 0;
            SpectrumFloat curWL = WavelengthLowBound;
            SpectrumFloat prev_xbarVal = 0, prev_ybarVal = 0, prev_zbarVal = 0;
            SpectrumFloat prevValue = 0;
            SpectrumFloat halfWidth = 0;
            CompensatedSum<SpectrumFloat> X(0), Y(0), Z(0);
            while (true) {
                SpectrumFloat xbarValue, ybarValue, zbarValue;
                if (curWL == WavelengthLowBound + curCMFIdx * CMFBinWidth) {
                    xbarValue = xbar_2deg[curCMFIdx];
                    ybarValue = ybar_2deg[curCMFIdx];
                    zbarValue = zbar_2deg[curCMFIdx];
                    ++curCMFIdx;
                }
                else {
                    uint32_t idx = std::min(uint32_t((curWL - WavelengthLowBound) / CMFBinWidth), NumCMFSamples - 1);
                    SpectrumFloat CMFBaseWL = WavelengthLowBound + idx * CMFBinWidth;
                    SpectrumFloat t = (curWL - CMFBaseWL) / CMFBinWidth;
                    xbarValue = (1 - t) * xbar_2deg[idx] + t * xbar_2deg[idx + 1];
                    ybarValue = (1 - t) * ybar_2deg[idx] + t * ybar_2deg[idx + 1];
                    zbarValue = (1 - t) * zbar_2deg[idx] + t * zbar_2deg[idx + 1];
                }
                
                SpectrumFloat value;
                if (curWL < minLambda) {
                    value = values[0];
                }
                else if (curWL > maxLambda) {
                    value = values[numSamples - 1];
                }
                else if (curWL == minLambda + baseIdx * binWidth) {
                    value = values[baseIdx];
                    ++baseIdx;
                }
                else {
                    uint32_t idx = std::min(uint32_t((curWL - minLambda) / binWidth), numSamples - 1);
                    SpectrumFloat baseWL = minLambda + idx * binWidth;
                    SpectrumFloat t = (curWL - baseWL) / binWidth;
                    value = (1 - t) * values[idx] + t * values[idx + 1];
                }
                
                SpectrumFloat avgValue = (prevValue + value) * 0.5f;
                X += avgValue * (prev_xbarVal + xbarValue) * halfWidth;
                Y += avgValue * (prev_ybarVal + ybarValue) * halfWidth;
                Z += avgValue * (prev_zbarVal + zbarValue) * halfWidth;
                
                prev_xbarVal = xbarValue;
                prev_ybarVal = ybarValue;
                prev_zbarVal = zbarValue;
                prevValue = value;
                SpectrumFloat prevWL = curWL;
                curWL = std::min(WavelengthLowBound + curCMFIdx * CMFBinWidth,
                                 baseIdx < numSamples ? (minLambda + baseIdx * binWidth) : INFINITY);
                halfWidth = (curWL - prevWL) * 0.5f;
                
                if (curCMFIdx == NumCMFSamples)
                    break;
            }
            XYZ[0] = X / integralCMF;
            XYZ[1] = Y / integralCMF;
            XYZ[2] = Z / integralCMF;
        }
        
        static void spectrum_to_XYZ(const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples, SpectrumFloat XYZ[3]) {
            const SpectrumFloat CMFBinWidth = (WavelengthHighBound - WavelengthLowBound) / (NumCMFSamples - 1);
            uint32_t curCMFIdx = 0;
            uint32_t baseIdx = 0;
            SpectrumFloat curWL = WavelengthLowBound;
            SpectrumFloat prev_xbarVal = 0, prev_ybarVal = 0, prev_zbarVal = 0;
            SpectrumFloat prevValue = 0;
            SpectrumFloat halfWidth = 0;
            CompensatedSum<SpectrumFloat> X(0), Y(0), Z(0);
            while (true) {
                SpectrumFloat xbarValue, ybarValue, zbarValue;
                if (curWL == WavelengthLowBound + curCMFIdx * CMFBinWidth) {
                    xbarValue = xbar_2deg[curCMFIdx];
                    ybarValue = ybar_2deg[curCMFIdx];
                    zbarValue = zbar_2deg[curCMFIdx];
                    ++curCMFIdx;
                }
                else {
                    uint32_t idx = std::min(uint32_t((curWL - WavelengthLowBound) / CMFBinWidth), NumCMFSamples - 1);
                    SpectrumFloat CMFBaseWL = WavelengthLowBound + idx * CMFBinWidth;
                    SpectrumFloat t = (curWL - CMFBaseWL) / CMFBinWidth;
                    xbarValue = (1 - t) * xbar_2deg[idx] + t * xbar_2deg[idx + 1];
                    ybarValue = (1 - t) * ybar_2deg[idx] + t * ybar_2deg[idx + 1];
                    zbarValue = (1 - t) * zbar_2deg[idx] + t * zbar_2deg[idx + 1];
                }
                
                SpectrumFloat value;
                if (curWL < lambdas[0]) {
                    value = values[0];
                }
                else if (curWL > lambdas[1]) {
                    value = values[numSamples - 1];
                }
                else if (curWL == lambdas[baseIdx]) {
                    value = values[baseIdx];
                    ++baseIdx;
                }
                else {
                    const SpectrumFloat* lb = std::lower_bound(lambdas + std::max(baseIdx - 1, 0u), lambdas + numSamples, curWL);
                    uint32_t idx = std::max(int32_t(std::distance(lambdas, lb)) - 1, 0);
                    SpectrumFloat t = (curWL - lambdas[idx]) / (lambdas[idx + 1] - lambdas[idx]);
                    value = (1 - t) * values[idx] + t * values[idx + 1];
                }
                
                SpectrumFloat avgValue = (prevValue + value) * 0.5f;
                X += avgValue * (prev_xbarVal + xbarValue) * halfWidth;
                Y += avgValue * (prev_ybarVal + ybarValue) * halfWidth;
                Z += avgValue * (prev_zbarVal + zbarValue) * halfWidth;
                
                prev_xbarVal = xbarValue;
                prev_ybarVal = ybarValue;
                prev_zbarVal = zbarValue;
                prevValue = value;
                SpectrumFloat prevWL = curWL;
                curWL = std::min(WavelengthLowBound + curCMFIdx * CMFBinWidth, baseIdx < numSamples ? lambdas[baseIdx] : INFINITY);
                halfWidth = (curWL - prevWL) * 0.5f;
                
                if (curCMFIdx == NumCMFSamples)
                    break;
            }
            XYZ[0] = X / integralCMF;
            XYZ[1] = Y / integralCMF;
            XYZ[2] = Z / integralCMF;
        }
        
        InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
            SLRAssert(e0 >= 0.0 && e1 >= 0.0 && e2 >= 0.0, "Values should not be minus.");
            switch (space) {
                case ColorSpace::sRGB:
                    return createShared<RGBInputSpectrum>(e0, e1, e2);
                case ColorSpace::sRGB_NonLinear: {
                    e0 = sRGB_degamma(e0);
                    e1 = sRGB_degamma(e1);
                    e2 = sRGB_degamma(e2);
                    return createShared<RGBInputSpectrum>(e0, e1, e2);
                }
                case ColorSpace::xyY: {
                    SpectrumFloat xyY[3] = {e0, e1, e2};
                    SpectrumFloat XYZ[3];
                    xyY_to_XYZ(xyY, XYZ);
                    e0 = XYZ[0];
                    e1 = XYZ[1];
                    e2 = XYZ[2];
                }
                case ColorSpace::XYZ: {
                    SpectrumFloat XYZ[3] = {e0, e1, e2};
                    SpectrumFloat RGB[3];
                    switch (spType) {
                        case SpectrumType::Reflectance:
                            XYZ_to_sRGB_E(XYZ, RGB);
                            break;
                        case SpectrumType::Illuminant:
                            XYZ_to_sRGB(XYZ, RGB);
                            break;
                        case SpectrumType::IndexOfRefraction:
                            XYZ_to_sRGB_E(XYZ, RGB);
                            break;
                        default:
                            break;
                    }
                    RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
                    RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
                    RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
                    return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
                }
                default:
                    SLRAssert(false, "Invalid color space.");
                    return createShared<RGBInputSpectrum>();
            }
        }
        InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
            SpectrumFloat XYZ[3];
            spectrum_to_XYZ(minLambda, maxLambda, values, numSamples, XYZ);
            SpectrumFloat RGB[3];
            switch (spType) {
                case SpectrumType::Reflectance:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                case SpectrumType::Illuminant:
                    XYZ_to_sRGB(XYZ, RGB);
                    break;
                case SpectrumType::IndexOfRefraction:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                default:
                    break;
            }
            RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
            RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
            RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
            return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
        }
        InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
            SpectrumFloat XYZ[3];
            spectrum_to_XYZ(lambdas, values, numSamples, XYZ);
            SpectrumFloat RGB[3];
            switch (spType) {
                case SpectrumType::Reflectance:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                case SpectrumType::Illuminant:
                    XYZ_to_sRGB(XYZ, RGB);
                    break;
                case SpectrumType::IndexOfRefraction:
                    XYZ_to_sRGB_E(XYZ, RGB);
                    break;
                default:
                    break;
            }
            RGB[0] = RGB[0] < 0.0f ? 0.0f : RGB[0];
            RGB[1] = RGB[1] < 0.0f ? 0.0f : RGB[1];
            RGB[2] = RGB[2] < 0.0f ? 0.0f : RGB[2];
            return createShared<RGBInputSpectrum>(RGB[0], RGB[1], RGB[2]);
        }
#endif
    } // namespace Spectrum
    
    namespace Image {
        using namespace SLR;
        std::map<std::string, Image2DRef> s_imageDB;
        
        std::shared_ptr<SLR::TiledImage2D> createTiledImage(const std::string &filepath, SLR::Allocator *mem, SLR::SpectrumType spType, bool gammaCorrection) {
            if (s_imageDB.count(filepath) > 0) {
                return std::static_pointer_cast<SLR::TiledImage2D>(s_imageDB[filepath]);
            }
            else {
                uint64_t requiredSize;
                bool imgSuccess;
                uint32_t width, height;
                ::ColorFormat colorFormat;
                imgSuccess = getImageInfo(filepath, &width, &height, &requiredSize, &colorFormat);
                SLRAssert(imgSuccess, "Error occured during getting image information.\n%s", filepath.c_str());
                
                void* linearData = malloc(requiredSize);
                imgSuccess = loadImage(filepath, (uint8_t*)linearData, gammaCorrection);
                SLRAssert(imgSuccess, "failed to load the image\n%s", filepath.c_str());
                
                SLR::ColorFormat internalFormat = (SLR::ColorFormat)colorFormat;
                TiledImage2D* texData = new SLR::TiledImage2D(linearData, width, height, internalFormat, mem, spType);
                free(linearData);
                
                std::shared_ptr<TiledImage2D> ret = std::shared_ptr<SLR::TiledImage2D>(texData);
                s_imageDB[filepath] = ret;
                return ret;
            }
        };

    } // namespace Image
}
