//
//  API.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "API.h"

#include <thread>

#include <libSLR/MemoryAllocators/ArenaAllocator.h>
#include <libSLR/BasicTypes/spectrum_library.h>
#include <libSLR/Core/transform.h>
#include <libSLR/Core/image_2d.h>
#include <libSLR/RNG/XORShiftRNG.h>
#include <libSLR/SurfaceShape/TriangleSurfaceShape.h>
#include <libSLR/Scene/Scene.h>
#include <libSLR/Renderer/DebugRenderer.h>
#include <libSLR/Renderer/PTRenderer.h>
#include <libSLR/Renderer/BPTRenderer.h>
#include <libSLR/Renderer/VolumetricPTRenderer.h>
#include <libSLR/Renderer/VolumetricBPTRenderer.h>

#include "images.h"
#include "textures.h"
#include "surface_materials.h"
#include "medium_materials.h"
#include "Scene/node.h"
#include "Scene/camera_nodes.h"
#include "Scene/TriangleMeshNode.h"
#include "Scene/medium_nodes.h"
#include "node_constructor.h"

#include "Parser/SceneParsingDriver.h"
#include "Parser/BuiltinFunctions/builtin_math.h"
#include "Parser/BuiltinFunctions/builtin_transform.h"
#include "Parser/BuiltinFunctions/builtin_texture.h"

namespace SLRSceneGraph {
    static bool strToImageStoreMode(const std::string &str, SLR::ImageStoreMode* mode) {
        if (str == "AsIs")
            *mode = SLR::ImageStoreMode::AsIs;
        else if (str == "Normal")
            *mode = SLR::ImageStoreMode::NormalTexture;
        else if (str == "Alpha")
            *mode = SLR::ImageStoreMode::AlphaTexture;
        else
            return false;
        return true;
    }
    
    static bool strToSpectrumType(const std::string &str, SLR::SpectrumType* type) {
        if (str == "Reflectance")
            *type = SLR::SpectrumType::Reflectance;
        else if (str == "Illuminant")
            *type = SLR::SpectrumType::Illuminant;
        else if (str == "IndexOfRefraction" || str == "IoR" || str == "IOR")
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
    
    SLR_SCENEGRAPH_API bool readScene(const std::string &filePath, const SceneRef &scene, RenderingContext* context) {
        TypeInfo::init();
        ExecuteContext executeContext;
        ErrorMessage errMsg;
        char curDir[256];
        SLR_getcwd(sizeof(curDir), curDir);
        std::string strCurDir = curDir;
        std::replace(strCurDir.begin(), strCurDir.end(), '\\', '/');
        std::string mFilePath = filePath;
        std::replace(mFilePath.begin(), mFilePath.end(), '\\', '/');
        std::string pathPrefix = mFilePath.substr(0, mFilePath.find_last_of("/") + 1);
        executeContext.absFileDirPath = strCurDir +"/" + pathPrefix;
        executeContext.scene = scene;
        executeContext.renderingContext = context;
        {
            LocalVariables &stack = executeContext.stackVariables.current();
            stack["root"] = Element::createFromReference<TypeMap::InternalNode>(scene->rootNode());
            
            stack["print"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"value", Type::Any}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::cout << args.at("value") << std::endl;
                                                   return Element();
                                               }
                                               );
            stack["addItem"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"tuple", Type::Tuple}, {"key", Type::String, Element::create<TypeMap::String>("")}, {"item", Type::Any}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   auto tuple = args.at("tuple").rawRef<TypeMap::Tuple>();
                                                   auto key = args.at("key").raw<TypeMap::String>();
                                                   tuple->add(key, args.at("item"));
                                                   
                                                   return Element::createFromReference<TypeMap::Tuple>(tuple);
                                               }
                                               );
            stack["numElements"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"tuple", Type::Tuple}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   auto tuple = args.at("tuple").rawRef<TypeMap::Tuple>();
                                                   return Element((TypeMap::Integer::InternalType)tuple->numParams());
                                               }
                                               );
            stack["Point"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   auto x = args.at("x").raw<TypeMap::RealNumber>();
                                                   auto y = args.at("y").raw<TypeMap::RealNumber>();
                                                   auto z = args.at("z").raw<TypeMap::RealNumber>();
                                                   return Element(SLR::Point3D(x, y, z));
                                               }
                                               );
            stack["Vector"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   auto x = args.at("x").raw<TypeMap::RealNumber>();
                                                   auto y = args.at("y").raw<TypeMap::RealNumber>();
                                                   auto z = args.at("z").raw<TypeMap::RealNumber>();
                                                   return Element(SLR::Vector3D(x, y, z));
                                               }
                                               );
            stack["getX"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {{"point", Type::Point}},
                                                   {{"vector", Type::Vector}},
                                                   {{"normal", Type::Normal}}
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &point = args.at("point").raw<TypeMap::Point>();
                                                       return Element(point.x);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &vector = args.at("vector").raw<TypeMap::Vector>();
                                                       return Element(vector.x);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &normal = args.at("normal").raw<TypeMap::Normal>();
                                                       return Element(normal.x);
                                                   }
                                               }
                                               );
            stack["getY"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {{"point", Type::Point}},
                                                   {{"vector", Type::Vector}},
                                                   {{"normal", Type::Normal}}
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &point = args.at("point").raw<TypeMap::Point>();
                                                       return Element(point.y);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &vector = args.at("vector").raw<TypeMap::Vector>();
                                                       return Element(vector.y);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &normal = args.at("normal").raw<TypeMap::Normal>();
                                                       return Element(normal.y);
                                                   }
                                               }
                                               );
            stack["getZ"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {{"point", Type::Point}},
                                                   {{"vector", Type::Vector}},
                                                   {{"normal", Type::Normal}}
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &point = args.at("point").raw<TypeMap::Point>();
                                                       return Element(point.z);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &vector = args.at("vector").raw<TypeMap::Vector>();
                                                       return Element(vector.z);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &normal = args.at("normal").raw<TypeMap::Normal>();
                                                       return Element(normal.z);
                                                   }
                                               }
                                               );
            
            stack["HSVtoRGB"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"H", Type::RealNumber}, {"S", Type::RealNumber}, {"V", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   auto H = args.at("H").raw<TypeMap::RealNumber>();
                                                   auto S = args.at("S").raw<TypeMap::RealNumber>();
                                                   auto V = args.at("V").raw<TypeMap::RealNumber>();
                                                   
                                                   float R = V;
                                                   float G = V;
                                                   float B = V;
                                                   if (S > 0.0f) {
                                                       H *= 6.0f;
                                                       int i = (int)H;
                                                       float f = H - (float)i;
                                                       switch (i) {
                                                           case 0:
                                                               G *= 1 - S * (1 - f);
                                                               B *= 1 - S;
                                                               break;
                                                           case 1:
                                                               R *= 1 - S * f;
                                                               B *= 1 - S;
                                                               break;
                                                           case 2:
                                                               R *= 1 - S;
                                                               B *= 1 - S * (1 - f);
                                                               break;
                                                           case 3:
                                                               R *= 1 - S;
                                                               G *= 1 - S * f;
                                                               break;
                                                           case 4:
                                                               R *= 1 - S * (1 - f);
                                                               G *= 1 - S;
                                                               break;
                                                           case 5:
                                                               G *= 1 - S;
                                                               B *= 1 - S * f;
                                                               break;
                                                           default:
                                                               break;
                                                       }
                                                   }
                                                   
                                                   ParameterListRef rgbTuple = createShared<ParameterList>();
                                                   rgbTuple->add("", Element(R));
                                                   rgbTuple->add("", Element(G));
                                                   rgbTuple->add("", Element(B));
                                                   
                                                   return Element::createFromReference<TypeMap::Tuple>(rgbTuple);
                                               }
                                               );
            
            stack["random"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   static SLR::XORShiftRNG rng{2112984105};
                                                   return Element(rng.getFloat0cTo1o());
                                               }
                                               );
            
            stack["abs"] = BuiltinFunctions::Math::abs;
            stack["min"] = BuiltinFunctions::Math::min;
            stack["max"] = BuiltinFunctions::Math::max;
            stack["clamp"] = BuiltinFunctions::Math::clamp;
            stack["sqrt"] = BuiltinFunctions::Math::sqrt;
            stack["pow"] = BuiltinFunctions::Math::pow;
            stack["exp"] = BuiltinFunctions::Math::exp;
            stack["ln"] = BuiltinFunctions::Math::ln;
            stack["log2"] = BuiltinFunctions::Math::log2;
            stack["log10"] = BuiltinFunctions::Math::log10;
            stack["sin"] = BuiltinFunctions::Math::sin;
            stack["cos"] = BuiltinFunctions::Math::cos;
            stack["tan"] = BuiltinFunctions::Math::tan;
            stack["asin"] = BuiltinFunctions::Math::asin;
            stack["acos"] = BuiltinFunctions::Math::acos;
            stack["atan"] = BuiltinFunctions::Math::atan;
            stack["dot"] = BuiltinFunctions::Math::dot;
            stack["cross"] = BuiltinFunctions::Math::cross;
            stack["distance"] = BuiltinFunctions::Math::distance;
            
            stack["translate"] = BuiltinFunctions::Transform::translate;
            stack["rotate"] = BuiltinFunctions::Transform::rotate;
            stack["rotateX"] = BuiltinFunctions::Transform::rotateX;
            stack["rotateY"] = BuiltinFunctions::Transform::rotateY;
            stack["rotateZ"] = BuiltinFunctions::Transform::rotateZ;
            stack["scale"] = BuiltinFunctions::Transform::scale;
            stack["lookAt"] = BuiltinFunctions::Transform::lookAt;
            stack["AnimatedTransform"] = BuiltinFunctions::Transform::AnimatedTransform;
            
            stack["Texture2DMapping"] = BuiltinFunctions::Texture::Texture2DMapping;
            stack["Texture3DMapping"] = BuiltinFunctions::Texture::Texture3DMapping;
            stack["SpectrumTexture"] = BuiltinFunctions::Texture::SpectrumTexture;
            stack["NormalTexture"] = BuiltinFunctions::Texture::NormalTexture;
            stack["FloatTexture"] = BuiltinFunctions::Texture::FloatTexture;
            
            stack["createVertex"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"position", Type::Tuple}, {"normal", Type::Tuple}, {"tangent", Type::Tuple}, {"texCoord", Type::Tuple}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   static const Function sigPosition = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                   static const Function sigNormal = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                   static const Function sigTangent = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                   static const Function sigTexCoord = Function(1, {{"u", Type::RealNumber}, {"v", Type::RealNumber}});
                                                   static const auto procPosition = [](const std::map<std::string, Element> &arg) {
                                                       return SLR::Point3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                   };
                                                   static const auto procNormal = [](const std::map<std::string, Element> &arg) {
                                                       return SLR::Normal3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                   };
                                                   static const auto procTangent = [](const std::map<std::string, Element> &arg) {
                                                       return SLR::Tangent3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                   };
                                                   static const auto procTexCoord = [](const std::map<std::string, Element> &arg) {
                                                       return SLR::TexCoord2D(arg.at("u").raw<TypeMap::RealNumber>(), arg.at("v").raw<TypeMap::RealNumber>());
                                                   };
                                                   return Element::create<TypeMap::Vertex>(sigPosition.perform<SLR::Point3D>(procPosition, args.at("position").raw<TypeMap::Tuple>()),
                                                                                           sigNormal.perform<SLR::Normal3D>(procNormal, args.at("normal").raw<TypeMap::Tuple>()),
                                                                                           sigTangent.perform<SLR::Tangent3D>(procTangent, args.at("tangent").raw<TypeMap::Tuple>()),
                                                                                           sigTexCoord.perform<SLR::TexCoord2D>(procTexCoord, args.at("texCoord").raw<TypeMap::Tuple>()));
                                               }
                                               );
            stack["Spectrum"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {
                                                       {"type", Type::String},
                                                       {"value", Type::RealNumber}
                                                   },
                                                   {
                                                       {"type", Type::String, Element::create<TypeMap::String>("Reflectance")},
                                                       {"space", Type::String, Element::create<TypeMap::String>("sRGB")},
                                                       {"e0", Type::RealNumber},
                                                       {"e1", Type::RealNumber},
                                                       {"e2", Type::RealNumber}
                                                   },
                                                   {
                                                       {"type", Type::String, Element::create<TypeMap::String>("Reflectance")},
                                                       {"minWL", Type::RealNumber},
                                                       {"maxWL", Type::RealNumber},
                                                       {"values", Type::Tuple}
                                                   },
                                                   {
                                                       {"type", Type::String, Element::create<TypeMap::String>("Reflectance")},
                                                       {"wls", Type::Tuple},
                                                       {"values", Type::Tuple}
                                                   },
                                                   {
                                                       {"type", Type::String},
                                                       {"name", Type::String},
                                                       {"idx", Type::Integer, Element(0)}
                                                   }
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       std::string typeStr = args.at("type").raw<TypeMap::String>();
                                                       SLR::SpectrumType type;
                                                       if (!strToSpectrumType(typeStr, &type)) {
                                                           *err = ErrorMessage("Specified spectrum type is invalid.");
                                                           return Element();
                                                       }
                                                       float value = args.at("value").raw<TypeMap::RealNumber>();
                                                       
                                                       return Element::createFromReference<TypeMap::Spectrum>(Spectrum::create(type, SLR::ColorSpace::sRGB, value, value, value));
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       std::string typeStr = args.at("type").raw<TypeMap::String>();
                                                       SLR::SpectrumType type;
                                                       if (!strToSpectrumType(typeStr, &type)) {
                                                           *err = ErrorMessage("Specified spectrum type is invalid.");
                                                           return Element();
                                                       }
                                                       std::string spaceStr = args.at("space").raw<TypeMap::String>();
                                                       SLR::ColorSpace space;
                                                       if (!strToColorSpace(spaceStr, &space)) {
                                                           *err = ErrorMessage("Specified color space is invalid.");
                                                           return Element();
                                                       }
                                                       float e0 = args.at("e0").raw<TypeMap::RealNumber>();
                                                       float e1 = args.at("e1").raw<TypeMap::RealNumber>();
                                                       float e2 = args.at("e2").raw<TypeMap::RealNumber>();
                                                       
                                                       return Element::createFromReference<TypeMap::Spectrum>(Spectrum::create(type, space, e0, e1, e2));
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       std::string typeStr = args.at("type").raw<TypeMap::String>();
                                                       SLR::SpectrumType type;
                                                       if (!strToSpectrumType(typeStr, &type)) {
                                                           *err = ErrorMessage("Specified spectrum type is invalid.");
                                                           return Element();
                                                       }
                                                       float minWL = args.at("minWL").raw<TypeMap::RealNumber>();
                                                       float maxWL = args.at("maxWL").raw<TypeMap::RealNumber>();
                                                       const ParameterList &valueList = args.at("values").raw<TypeMap::Tuple>();
                                                       size_t numSamples = valueList.numUnnamed();
                                                       std::vector<float> values;
                                                       values.resize(numSamples);
                                                       for (int i = 0; i < numSamples; ++i) {
                                                           const Element &el = valueList(i);
                                                           values[i] = el.asRaw<TypeMap::RealNumber>();
                                                       }
                                                       
                                                       return Element::createFromReference<TypeMap::Spectrum>(Spectrum::create(type, minWL, maxWL, values.data(), (uint32_t)numSamples));
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       std::string typeStr = args.at("type").raw<TypeMap::String>();
                                                       SLR::SpectrumType type;
                                                       if (!strToSpectrumType(typeStr, &type)) {
                                                           *err = ErrorMessage("Specified spectrum type is invalid.");
                                                           return Element();
                                                       }
                                                       const ParameterList &wavelengthList = args.at("wls").raw<TypeMap::Tuple>();
                                                       const ParameterList &valueList = args.at("values").raw<TypeMap::Tuple>();
                                                       size_t numSamples = wavelengthList.numUnnamed();
                                                       if (numSamples != valueList.numUnnamed()) {
                                                           *err = ErrorMessage("The sizes of the wavelengths and the values are different.");
                                                           return Element();
                                                       }
                                                       std::vector<float> wls;
                                                       std::vector<float> values;
                                                       wls.resize(numSamples);
                                                       values.resize(numSamples);
                                                       for (int i = 0; i < numSamples; ++i) {
                                                           const Element &elWavelength = wavelengthList(i);
                                                           const Element &elValue = valueList(i);
                                                           wls[i] = elWavelength.asRaw<TypeMap::RealNumber>();
                                                           values[i] = elValue.asRaw<TypeMap::RealNumber>();
                                                       }
                                                       
                                                       return Element::createFromReference<TypeMap::Spectrum>(Spectrum::create(type, wls.data(), values.data(), (uint32_t)numSamples));
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err)  {
                                                       using namespace SLR;
                                                       AssetSpectrumRef spectrum;
                                                       std::string type = args.at("type").raw<TypeMap::String>();
                                                       std::string name = args.at("name").raw<TypeMap::String>();
                                                       auto idx = args.at("idx").raw<TypeMap::Integer>();
                                                       SpectrumType spType = SpectrumType::Illuminant;
                                                       SpectrumLibrary::Data data;
                                                       bool success = true;
                                                       if (type == "Illuminant") {
                                                           success = SpectrumLibrary::queryIlluminantSpectrum(name, idx, &data);
                                                           spType = SpectrumType::Illuminant;
                                                       }
                                                       else if (type == "Reflectance") {
                                                           success = SpectrumLibrary::queryReflectanceSpectrum(name, idx, &data);
                                                           spType = SpectrumType::Reflectance;
                                                       }
                                                       else if (type == "IoR") {
                                                           success = SpectrumLibrary::queryIoRSpectrum(name, idx, &data);
                                                           spType = SpectrumType::IndexOfRefraction;
                                                       }
                                                       if (success) {
                                                           if (data.dType == SpectrumLibrary::DistributionType::Regular)
                                                               spectrum = Spectrum::create(spType, data.minLambdas, data.maxLambdas, data.values, data.numSamples);
                                                           else
                                                               spectrum = Spectrum::create(spType, data.lambdas, data.values, data.numSamples);
                                                       }
                                                       else {
                                                           *err = ErrorMessage("Spectrum %s not found.", name.c_str());
                                                           return Element();
                                                       }
                                                       return Element::createFromReference<TypeMap::Spectrum>(spectrum);
                                                   }
                                               }
                                               );
            stack["scaleAndOffset"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"spectrum", Type::Spectrum},
                                                   {"scale", Type::RealNumber},
                                                   {"offset", Type::RealNumber}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   AssetSpectrumRef spectrum = args.at("spectrum").rawRef<TypeMap::Spectrum>();
                                                   float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                                   float offset = args.at("offset").raw<TypeMap::RealNumber>();
                                                   return Element::createFromReference<TypeMap::Spectrum>(AssetSpectrumRef(spectrum->createScaledAndOffset(scale, offset)));
                                               }
                                               );
            stack["Image2D"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"path", Type::String},
                                                   {"mode", Type::String, Element::create<TypeMap::String>("AsIs")},
                                                   {"type", Type::String, Element::create<TypeMap::String>("Reflectance")}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string path = args.at("path").raw<TypeMap::String>();
                                                   std::string modeStr = args.at("mode").raw<TypeMap::String>();
                                                   std::string typeStr = args.at("type").raw<TypeMap::String>();
                                                   SLR::ImageStoreMode mode;
                                                   if (!strToImageStoreMode(modeStr, &mode)) {
                                                       *err = ErrorMessage("Specified image store mode is invalid.");
                                                       return Element();
                                                   }
                                                   SLR::SpectrumType spType;
                                                   if (!strToSpectrumType(typeStr, &spType)) {
                                                       *err = ErrorMessage("Specified spectrum type is invalid.");
                                                       return Element();
                                                   }
                                                   
                                                   return Element::createFromReference<TypeMap::Image2D>(createImage2D(path, mode, spType, false));
                                               }
                                               );
            
            stack["createSurfaceMaterial"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"type", Type::String}, {"params", Type::Tuple}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string type = args.at("type").raw<TypeMap::String>();
                                                   const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                   if (type == "matte") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"reflectance", Type::SpectrumTexture},
                                                               {"sigma", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(nullptr)}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef reflectance = args.at("reflectance").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef sigma = args.at("sigma").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createMatte(reflectance, sigma));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "metal") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"coeffR", Type::SpectrumTexture},
                                                               {"eta", Type::SpectrumTexture}, {"k", Type::SpectrumTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef coeffR = args.at("coeffR").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef eta = args.at("eta").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef k = args.at("k").rawRef<TypeMap::SpectrumTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createMetal(coeffR, eta, k));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "glass") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"coeff", Type::SpectrumTexture},
                                                               {"etaExt", Type::SpectrumTexture}, {"etaInt", Type::SpectrumTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef coeff = args.at("coeff").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef etaExt = args.at("etaExt").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef etaInt = args.at("etaInt").rawRef<TypeMap::SpectrumTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createGlass(coeff, etaExt, etaInt));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "Ward") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"R", Type::SpectrumTexture},
                                                               {"anisoX", Type::FloatTexture}, {"anisoY", Type::FloatTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef R = args.at("R").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef anisoX = args.at("anisoX").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef anisoY = args.at("anisoY").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createModifiedWardDur(R, anisoX, anisoY));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "Ashikhmin") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"Rd", Type::SpectrumTexture}, {"Rs", Type::SpectrumTexture},
                                                               {"nx", Type::FloatTexture}, {"ny", Type::FloatTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef Rd = args.at("Rd").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef Rs = args.at("Rs").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef nx = args.at("nx").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef ny = args.at("ny").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createAshikhminShirley(Rd, Rs, nx, ny));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "microfacet metal") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"eta", Type::SpectrumTexture}, {"k", Type::SpectrumTexture},
                                                               {"alpha_gx", Type::FloatTexture}, {"alpha_gy", Type::FloatTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef eta = args.at("eta").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef k = args.at("k").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef alpha_gx = args.at("alpha_gx").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef alpha_gy = args.at("alpha_gy").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createMicrofacetMetal(eta, k, alpha_gx, alpha_gy));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "microfacet glass") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"etaExt", Type::SpectrumTexture}, {"etaInt", Type::SpectrumTexture},
                                                               {"alpha_gx", Type::FloatTexture}, {"alpha_gy", Type::FloatTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef etaExt = args.at("etaExt").rawRef<TypeMap::SpectrumTexture>();
                                                               SpectrumTextureRef etaInt = args.at("etaInt").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef alpha_gx = args.at("alpha_gx").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef alpha_gy = args.at("alpha_gy").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createMicrofacetGlass(etaExt, etaInt, alpha_gx, alpha_gy));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "disney reflection") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"baseColor", Type::SpectrumTexture},
                                                               {"subsurface", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"metallic", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"specular", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.5f))}, 
                                                               {"specularTint", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"roughness", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.5f))}, 
                                                               {"anisotropic", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"sheen", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"sheenTint", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.5f))}, 
                                                               {"clearcoat", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(0.0f))}, 
                                                               {"clearcoatGloss", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(createShared<ConstantFloatTexture>(1.0f))}, 
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef baseColor = args.at("baseColor").rawRef<TypeMap::SpectrumTexture>();
                                                               FloatTextureRef subsurface = args.at("subsurface").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef metallic = args.at("metallic").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef specular = args.at("specular").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef specularTint = args.at("specularTint").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef roughness = args.at("roughness").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef anisotropic = args.at("anisotropic").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef sheen = args.at("sheen").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef sheenTint = args.at("sheenTint").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef clearcoat = args.at("clearcoat").rawRef<TypeMap::FloatTexture>();
                                                               FloatTextureRef clearcoatGloss = args.at("clearcoatGloss").rawRef<TypeMap::FloatTexture>();
                                                               
                                                               SurfaceMaterialRef ret = SurfaceMaterial::createDisneyReflection(
                                                                                            baseColor, 
                                                                                            subsurface, metallic, specular, specularTint, 
                                                                                            roughness, anisotropic, sheen, sheenTint, 
                                                                                            clearcoat, clearcoatGloss);
                                                               
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(ret);
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "flipped") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"base", Type::SurfaceMaterial}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SurfaceMaterialRef base = args.at("base").rawRef<TypeMap::SurfaceMaterial>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createFlippedMaterial(base));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "emitter") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"scatter", Type::SurfaceMaterial},
                                                               {"emitter", Type::EmitterSurfaceProperty}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SurfaceMaterialRef scatter = args.at("scatter").rawRef<TypeMap::SurfaceMaterial>();
                                                               EmitterSurfacePropertyRef emitter = args.at("emitter").rawRef<TypeMap::EmitterSurfaceProperty>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createEmitterSurfaceMaterial(scatter, emitter));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "mix") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"mat0", Type::SurfaceMaterial}, {"mat1", Type::SurfaceMaterial},
                                                               {"factor", Type::FloatTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SurfaceMaterialRef mat0 = args.at("mat0").rawRef<TypeMap::SurfaceMaterial>();
                                                               SurfaceMaterialRef mat1 = args.at("mat1").rawRef<TypeMap::SurfaceMaterial>();
                                                               FloatTextureRef factor = args.at("factor").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createMixedMaterial(mat0, mat1, factor));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "sum") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"mat0", Type::SurfaceMaterial}, {"mat1", Type::SurfaceMaterial}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SurfaceMaterialRef mat0 = args.at("mat0").rawRef<TypeMap::SurfaceMaterial>();
                                                               SurfaceMaterialRef mat1 = args.at("mat1").rawRef<TypeMap::SurfaceMaterial>();
                                                               return Element::createFromReference<TypeMap::SurfaceMaterial>(SurfaceMaterial::createSummedMaterial(mat0, mat1));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   *err = ErrorMessage("Specified material type is invalid.");
                                                   return Element();
                                               }
                                               );
            stack["createEmitterSurfaceProperty"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"type", Type::String}, {"params", Type::Tuple}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string type = args.at("type").raw<TypeMap::String>();
                                                   const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                   if (type == "diffuse") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"emittance", Type::SpectrumTexture}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef emittance = args.at("emittance").rawRef<TypeMap::SpectrumTexture>();
                                                               return Element::createFromReference<TypeMap::EmitterSurfaceProperty>(SurfaceMaterial::createDiffuseEmitter(emittance));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "ideal directional") {
                                                       const static Function configFunc{
                                                           0, {
                                                               {"emittance", Type::SpectrumTexture},
                                                               {"direction", Type::Vector, SLR::Vector3D(0, 0, 1)}
                                                           },
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               SpectrumTextureRef emittance = args.at("emittance").rawRef<TypeMap::SpectrumTexture>();
                                                               SLR::Vector3D direction = SLR::normalize(args.at("direction").raw<TypeMap::Vector>());
                                                               return Element::createFromReference<TypeMap::EmitterSurfaceProperty>(SurfaceMaterial::createIdealDirectionalEmitter(emittance, direction));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   *err = ErrorMessage("Specified material type is invalid.");
                                                   return Element();
                                               }
                                               );
            stack["createMediumMaterial"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"type", Type::String}, {"params", Type::Tuple}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string type = args.at("type").raw<TypeMap::String>();
                                                   const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                   if (type == "isotropic") {
                                                       const static Function configFunc{
                                                           0, {},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               return Element::createFromReference<TypeMap::MediumMaterial>(MediumMaterial::createIsotropic());
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   else if (type == "Henyey-Greenstein") {
                                                       const static Function configFunc{
                                                           0, {{"g", Type::FloatTexture}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               FloatTextureRef g = args.at("g").rawRef<TypeMap::FloatTexture>();
                                                               return Element::createFromReference<TypeMap::MediumMaterial>(MediumMaterial::createHenyeyGreenstein(g));
                                                           }
                                                       };
                                                       return configFunc(params, context, err);
                                                   }
                                                   *err = ErrorMessage("Specified material type is invalid.");
                                                   return Element();
                                               }
                                               );
            stack["createMesh"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"vertices", Type::Tuple},
                                                   {"matGroups", Type::Tuple},
                                                   {"axisForRadialTangent", Type::Integer, -1}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const std::vector<Element> &vertices = args.at("vertices").raw<TypeMap::Tuple>().unnamed;
                                                   const std::vector<Element> &matGroups = args.at("matGroups").raw<TypeMap::Tuple>().unnamed;
                                                   int32_t axisForRadialTangent = args.at("axisForRadialTangent").raw<TypeMap::Integer>();
                                                   
                                                   TriangleMeshNode::MaterialGroup resultMatGroup;
                                                   
                                                   static const Function sigMatGroup{
                                                       1, {
                                                           {"mat", Type::SurfaceMaterial},
                                                           {"normal", Type::NormalTexture, Element::createFromReference<TypeMap::NormalTexture>(nullptr)},
                                                           {"alpha", Type::FloatTexture, Element::createFromReference<TypeMap::FloatTexture>(nullptr)},
                                                           {"triangles", Type::Tuple}
                                                       }
                                                   };
                                                   static const auto procMatGroup = [&resultMatGroup, &err](const std::map<std::string, Element> &args) {
                                                       resultMatGroup.material = args.at("mat").rawRef<TypeMap::SurfaceMaterial>();
                                                       resultMatGroup.normalMap = args.at("normal").rawRef<TypeMap::NormalTexture>();
                                                       resultMatGroup.alphaMap = args.at("alpha").rawRef<TypeMap::FloatTexture>();
                                                       
                                                       static const Function sigTriangle{
                                                           1, {{"v0", Type::Integer}, {"v1", Type::Integer}, {"v2", Type::Integer}}
                                                       };
                                                       static const auto procTriangle = [](const std::map<std::string, Element> &args) {
                                                           return Triangle(args.at("v0").raw<TypeMap::Integer>(),
                                                                           args.at("v1").raw<TypeMap::Integer>(),
                                                                           args.at("v2").raw<TypeMap::Integer>());
                                                       };
                                                       
                                                       const std::vector<Element> &triangles = args.at("triangles").raw<TypeMap::Tuple>().unnamed;
                                                       for (int i = 0; i < triangles.size(); ++i) {
                                                           resultMatGroup.triangles.push_back(sigTriangle.perform<Triangle>(procTriangle, triangles[i].raw<TypeMap::Tuple>(), Triangle(), err));
                                                       }
                                                       
                                                       return 0;
                                                   };
                                                   
                                                   TriangleMeshNodeRef mesh = createShared<TriangleMeshNode>();
                                                   for (int i = 0; i < vertices.size(); ++i) {
                                                       if (vertices[i].type == Type::Vertex) {
                                                           mesh->addVertex(vertices[i].raw<TypeMap::Vertex>());
                                                       }
                                                       else {
                                                           const Function &CreateVertex = context.stackVariables["createVertex"].raw<TypeMap::Function>();
                                                           Element vtx = CreateVertex(vertices[i].raw<TypeMap::Tuple>(), context, err);
                                                           if (err->error)
                                                               return Element();
                                                           mesh->addVertex(vtx.raw<TypeMap::Vertex>());
                                                       }
                                                   }
                                                   for (int i = 0; i < matGroups.size(); ++i) {
                                                       resultMatGroup.material = nullptr;
                                                       resultMatGroup.normalMap = nullptr;
                                                       resultMatGroup.alphaMap = nullptr;
                                                       resultMatGroup.triangles.clear();
                                                       sigMatGroup.perform<uint32_t>(procMatGroup, matGroups[i].raw<TypeMap::Tuple>(), 0, err);
                                                       if (err->error)
                                                           return Element();
                                                       mesh->addMaterialGroup(resultMatGroup.material, resultMatGroup.normalMap, resultMatGroup.alphaMap, std::move(resultMatGroup.triangles));
                                                   }
                                                   
                                                   if (axisForRadialTangent < -1 || axisForRadialTangent > 2)
                                                       axisForRadialTangent = -1;
                                                   mesh->setAxisForRadialTangent(axisForRadialTangent);
                                                   
                                                   return Element::createFromReference<TypeMap::SurfaceNode>(mesh);
                                               }
                                               );
            stack["createInfinitesimalPoint"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"position", Type::Point}, {"direction", Type::Vector}, {"material", Type::SurfaceMaterial}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const auto &position = args.at("position").raw<TypeMap::Point>();
                                                   const auto &direction = args.at("direction").raw<TypeMap::Point>();
                                                   const SurfaceMaterialRef material = args.at("material").rawRef<TypeMap::SurfaceMaterial>();
                                                   
                                                   SurfaceNodeRef node = createShared<InfinitesimalPointNode>(position, direction, material);
                                                   
                                                   return Element::createFromReference<TypeMap::SurfaceNode>(node);
                                               }
                                               );
            stack["createVacuum"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"min", Type::Point}, {"max", Type::Point}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const auto &minP = args.at("min").raw<TypeMap::Point>();
                                                   const auto &maxP = args.at("max").raw<TypeMap::Point>();
                                                   MediumNodeRef mediumNode = createShared<VacuumMediumNode>(SLR::BoundingBox3D(minP, maxP));
                                                   return Element::createFromReference<TypeMap::MediumNode>(mediumNode);
                                               }
                                               );
            stack["createHomogeneousMedium"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"min", Type::Point}, {"max", Type::Point}, {"sigma_s", Type::Spectrum}, {"sigma_e", Type::Spectrum}, {"mat", Type::MediumMaterial}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const auto &minP = args.at("min").raw<TypeMap::Point>();
                                                   const auto &maxP = args.at("max").raw<TypeMap::Point>();
                                                   AssetSpectrumRef sigma_s = args.at("sigma_s").rawRef<TypeMap::Spectrum>();
                                                   AssetSpectrumRef sigma_e = args.at("sigma_e").rawRef<TypeMap::Spectrum>();
                                                   MediumMaterialRef mat = args.at("mat").rawRef<TypeMap::MediumMaterial>();
                                                   MediumNodeRef mediumNode = createShared<HomogeneousMediumNode>(SLR::BoundingBox3D(minP, maxP), sigma_s, sigma_e, mat);
                                                   return Element::createFromReference<TypeMap::MediumNode>(mediumNode);
                                               }
                                               );
            stack["createGridMedium"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {
                                                       {"min", Type::Point}, {"max", Type::Point},
                                                       {"base_sigma_s", Type::Spectrum}, {"base_sigma_e", Type::Spectrum},
                                                       {"density_grid", Type::Tuple}, {"numX", Type::Integer}, {"numY", Type::Integer}, {"numZ", Type::Integer},
                                                       {"mat", Type::MediumMaterial}
                                                   },
                                                   {
                                                       {"min", Type::Point}, {"max", Type::Point},
                                                       {"sigma_s", Type::Tuple}, {"sigma_e", Type::Tuple}, {"mat", Type::MediumMaterial}
                                                   }
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       const auto &minP = args.at("min").raw<TypeMap::Point>();
                                                       const auto &maxP = args.at("max").raw<TypeMap::Point>();
                                                       AssetSpectrumRef base_sigma_s = args.at("base_sigma_s").rawRef<TypeMap::Spectrum>();
                                                       AssetSpectrumRef base_sigma_e = args.at("base_sigma_e").rawRef<TypeMap::Spectrum>();
                                                       const ParameterList &density_grid = args.at("density_grid").raw<TypeMap::Tuple>();
                                                       uint32_t numX = args.at("numX").raw<TypeMap::Integer>();
                                                       uint32_t numY = args.at("numY").raw<TypeMap::Integer>();
                                                       uint32_t numZ = args.at("numZ").raw<TypeMap::Integer>();
                                                       SLRAssert(numZ == density_grid.numUnnamed(), "The number of z-slices of density_grid and specified grid z dimension do not match.");
                                                       std::vector<std::vector<float>> densityArray;
                                                       densityArray.resize(numZ);
                                                       for (int z = 0; z < numZ; ++z) {
                                                           const ParameterList &zSlice = density_grid(z).raw<TypeMap::Tuple>();
                                                           SLRAssert(numY * numX == zSlice.numUnnamed(), "The number of elements of a z-slice of density_grid and specified grid x dimension times y dimension do not match.");
                                                           std::vector<float> zSliceArray;
                                                           zSliceArray.resize(numY * numX);
                                                           for (int y = 0; y < numY; ++y) {
                                                               for (int x = 0; x < numX; ++x) {
                                                                   float density = zSlice(numX * y + x).asRaw<TypeMap::RealNumber>();
                                                                   zSliceArray[numX * y + x] = density;
                                                               }
                                                           }
                                                           densityArray[z] = std::move(zSliceArray);
                                                       }
                                                       MediumMaterialRef mat = args.at("mat").rawRef<TypeMap::MediumMaterial>();
                                                       MediumNodeRef mediumNode = createShared<DensityGridMediumNode>(SLR::BoundingBox3D(minP, maxP), base_sigma_s, base_sigma_e, 
                                                                                                                      std::move(densityArray), numX, numY, numZ, mat);
                                                       return Element::createFromReference<TypeMap::MediumNode>(mediumNode);
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       SLRAssert_NotImplemented();
                                                       return Element();
                                                   }
                                               }
                                               );
            stack["createNode"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   return Element::create<TypeMap::InternalNode>(createShared<SLR::StaticTransform>());
                                               }
                                               );
            stack["copyNode"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"src", Type::Node}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   NodeRef node = args.at("src").rawRef<TypeMap::Node>();
                                                   NodeRef copied = std::dynamic_pointer_cast<Node>(node->copy());
                                                   return Element::createFromReference<TypeMap::Node>(copied);
                                               }
                                               );
            stack["createReferenceNode"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"node", Type::Node}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   NodeRef node = args.at("node").rawRef<TypeMap::Node>();
                                                   NodeRef refNode = createShared<ReferenceNode>(node);
                                                   return Element::createFromReference<TypeMap::Node>(refNode);
                                               }
                                               );
            stack["setTransform"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"node", Type::InternalNode}, {"transform", Type::Transform}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   InternalNodeRef node = args.at("node").rawRef<TypeMap::InternalNode>();
                                                   TransformRef tf = args.at("transform").rawRef<TypeMap::Transform>();
                                                   node->setTransform(tf);
                                                   return Element();
                                               }
                                               );
            stack["addChild"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"parent", Type::InternalNode}, {"child", Type::Node}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   InternalNodeRef parent = args.at("parent").rawRef<TypeMap::InternalNode>();
                                                   NodeRef child = args.at("child").rawRef<TypeMap::Node>();
                                                   parent->addChildNode(child);
                                                   return Element();
                                               }
                                               );
            stack["setInternalMedium"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"surface", Type::SurfaceNode}, {"medium", Type::MediumNode}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   SurfaceNodeRef surface = args.at("surface").rawRef<TypeMap::SurfaceNode>();
                                                   MediumNodeRef medium = args.at("medium").rawRef<TypeMap::MediumNode>();
                                                   surface->setInternalMedium(medium);
                                                   return Element();
                                               }
                                               );
            stack["load3DModel"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"path", Type::String},
                                                   {"matProc", Type::Function, Element::createFromReference<TypeMap::Function>(nullptr)},
                                                   {"meshProc", Type::Function, Element::createFromReference<TypeMap::Function>(nullptr)}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string path = context.absFileDirPath + args.at("path").raw<TypeMap::String>();
                                                   auto userMatProcRef = args.at("matProc").rawRef<TypeMap::Function>();
                                                   auto meshProcRef = args.at("meshProc").rawRef<TypeMap::Function>();
                                                   
                                                   CreateMaterialFunction matProc = createMaterialDefaultFunction;
                                                   if (userMatProcRef) {
                                                       matProc = [&userMatProcRef, &context, &err](const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem) {
                                                           return createMaterialFunction(*userMatProcRef.get(), context, err, aiMat, pathPrefix, mem);
                                                       };
                                                   }
                                                   
                                                   MeshCallback meshCallback = meshCallbackDefaultFunction;
                                                   if (meshProcRef) {
                                                       meshCallback = [&meshProcRef, &context, &err](const std::string &name, const TriangleMeshNodeRef &mesh, const SLR::Point3D &minP, const SLR::Point3D &maxP) {
                                                           return meshCallbackFunction(*meshProcRef.get(), context, err, name, mesh, minP, maxP);
                                                       };
                                                   }
                                                   
                                                   InternalNodeRef modelNode;
                                                   construct(path, modelNode, matProc, meshCallback);
                                                   if (!modelNode) {
                                                       *err = ErrorMessage("Some errors occur during loading a 3D model.");
                                                       return Element();
                                                   }
                                                   modelNode->setName(path);
                                                   
                                                   return Element::createFromReference<TypeMap::Node>(modelNode);
                                               }
                                               );
            stack["scanXZFromYPlus"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"node", Type::Node},
                                                   {"numX", Type::Integer},
                                                   {"numY", Type::Integer},
                                                   {"randomness", Type::RealNumber, 0.0f},
                                                   {"callback", Type::Function}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   using namespace SLR;
                                                   NodeRef node = args.at("node").rawRef<TypeMap::Node>();
                                                   uint32_t numX = args.at("numX").raw<TypeMap::Integer>();
                                                   uint32_t numY = args.at("numY").raw<TypeMap::Integer>();
                                                   float randomness = args.at("randomness").raw<TypeMap::RealNumber>();
                                                   const Function &callback = args.at("callback").raw<TypeMap::Function>();
                                                   
                                                   SLR::XORShiftRNG rng(50287412);
                                                   
                                                   node->prepareForRendering();
                                                   SLR::Node &rawNode = *node->getRaw();
                                                   SLR::RenderingData renderingData(nullptr);
                                                   SLR::ArenaAllocator mem;
                                                   rawNode.createRenderingData(&mem, nullptr, &renderingData);
                                                   auto aggregate = createUnique<SurfaceObjectAggregate>(renderingData.surfObjs);
                                                   
                                                   SLR::BoundingBox3D bounds = aggregate->bounds();
                                                   for (int i = 0; i < numY; ++i) {
                                                       for (int j = 0; j < numX; ++j) {
                                                           float purturbX = randomness * (rng.getFloat0cTo1o() - 0.5f);
                                                           float purturbZ = randomness * (rng.getFloat0cTo1o() - 0.5f);
                                                           Ray ray(SLR::Point3D(bounds.minP.x + (bounds.maxP.x - bounds.minP.x) * (j + 0.5f + purturbX) / numX,
                                                                           bounds.maxP.y * 1.5f,
                                                                           bounds.minP.z + (bounds.maxP.z - bounds.minP.z) * (i + 0.5f + purturbZ) / numY),
                                                                   SLR::Vector3D(0, -1, 0), 0.0f);
                                                           SurfaceInteraction si;
                                                           if (!aggregate->intersect(ray, RaySegment(), &si))
                                                               continue;
                                                           SurfacePoint surfPt;
                                                           si.calculateSurfacePoint(&surfPt);
                                                           
                                                           ParameterList params;
                                                           
                                                           SLR::ReferenceFrame shadingFrame = surfPt.getShadingFrame();
                                                           Element elPosition = Element(surfPt.getPosition());
                                                           Element elNormal = Element(Normal3D(shadingFrame.z));
                                                           Element elTangent = Element(shadingFrame.x);
                                                           Element elBitangent = Element(shadingFrame.y);
                                                           params.add("", elPosition);
                                                           params.add("", elTangent);
                                                           params.add("", elBitangent);
                                                           params.add("", elNormal);
                                                           
                                                           callback(params, context, err);
                                                       }
                                                   }
                                                   
                                                   rawNode.destroyRenderingData(&mem);
                                                   
                                                   return Element();
                                               }
                                               );
            stack["createPerspectiveCamera"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"sensitivity", Type::RealNumber, Element(0.0)},
                                                   {"aspect", Type::RealNumber, Element(1.0)},
                                                   {"fovY", Type::RealNumber, Element(0.5235987756)},
                                                   {"radius", Type::RealNumber, Element(0.0)},
                                                   {"imgDist", Type::RealNumber, Element(0.02)},
                                                   {"objDist", Type::RealNumber, Element(5.0)}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float sensitivity = args.at("sensitivity").raw<TypeMap::RealNumber>();
                                                   float aspect = args.at("aspect").raw<TypeMap::RealNumber>();
                                                   float fovY = args.at("fovY").raw<TypeMap::RealNumber>();
                                                   float radius = args.at("radius").raw<TypeMap::RealNumber>();
                                                   float imgDist = args.at("imgDist").raw<TypeMap::RealNumber>();
                                                   float objDist = args.at("objDist").raw<TypeMap::RealNumber>();
                                                   
                                                   NodeRef rawRef = createShared<PerspectiveCameraNode>(sensitivity, aspect, fovY, radius, imgDist, objDist);
                                                   return Element::createFromReference<TypeMap::Node>(rawRef);
                                               }
                                               );
            stack["createEquirectangularCamera"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"sensitivity", Type::RealNumber, Element(0.0)},
                                                   {"phi", Type::RealNumber},
                                                   {"theta", Type::RealNumber}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float sensitivity = args.at("sensitivity").raw<TypeMap::RealNumber>();
                                                   float phi = args.at("phi").raw<TypeMap::RealNumber>();
                                                   float theta = args.at("theta").raw<TypeMap::RealNumber>();
                                                   
                                                   NodeRef rawRef = createShared<EquirectangularCameraNode>(sensitivity, phi, theta);
                                                   return Element::createFromReference<TypeMap::Node>(rawRef);
                                               }
                                               );
            stack["setRenderer"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"method", Type::String}, {"config", Type::Tuple, Element::create<TypeMap::Tuple>()}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   std::string method = args.at("method").raw<TypeMap::String>();
                                                   const ParameterList &config = args.at("config").raw<TypeMap::Tuple>();
                                                   if (method == "PT") {
                                                       const static Function configPT{
                                                           0, {{"samples", Type::Integer, Element(8)}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               uint32_t spp = args.at("samples").raw<TypeMap::Integer>();
                                                               context.renderingContext->renderer = createUnique<SLR::PTRenderer>(spp);
                                                               return Element();
                                                           }
                                                       };
                                                       return configPT(config, context, err);
                                                   }
                                                   else if (method == "BPT") {
                                                       const static Function configBPT{
                                                           0, {{"samples", Type::Integer, Element(8)}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               uint32_t spp = args.at("samples").raw<TypeMap::Integer>();
                                                               context.renderingContext->renderer = createUnique<SLR::BPTRenderer>(spp);
                                                               return Element();
                                                           }
                                                       };
                                                       return configBPT(config, context, err);
                                                   }
                                                   else if (method == "Volumetric PT") {
                                                       const static Function configVolumetricPT{
                                                           0, {{"samples", Type::Integer, Element(8)}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               uint32_t spp = args.at("samples").raw<TypeMap::Integer>();
                                                               context.renderingContext->renderer = createUnique<SLR::VolumetricPTRenderer>(spp);
                                                               return Element();
                                                           }
                                                       };
                                                       return configVolumetricPT(config, context, err);
                                                   }
                                                   else if (method == "Volumetric BPT") {
                                                       const static Function configVolumetricBPT{
                                                           0, {{"samples", Type::Integer, Element(8)}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               uint32_t spp = args.at("samples").raw<TypeMap::Integer>();
                                                               context.renderingContext->renderer = createUnique<SLR::VolumetricBPTRenderer>(spp);
                                                               return Element();
                                                           }
                                                       };
                                                       return configVolumetricBPT(config, context, err);
                                                   }
                                                   else if (method == "debug") {
                                                       const static Function configDebug{
                                                           0, {{"outputs", Type::Tuple}},
                                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                               const ParameterList &outputs = args.at("outputs").raw<TypeMap::Tuple>();
                                                               bool chFlags[(int)SLR::ExtraChannel::NumChannels];
                                                               for (int i = 0; i < (int)SLR::ExtraChannel::NumChannels; ++i)
                                                                   chFlags[i] = false;
                                                               for (int i = 0; i < outputs.numUnnamed(); ++i) {
                                                                   Element elem = outputs(i);
                                                                   if (elem.type != SLRSceneGraph::Type::String)
                                                                       continue;
                                                                   std::string chName = elem.raw<TypeMap::String>();
                                                                   if (chName == "geometric normal")
                                                                       chFlags[(int)SLR::ExtraChannel::GeometricNormal] = true;
                                                                   else if (chName == "shading normal")
                                                                       chFlags[(int)SLR::ExtraChannel::ShadingNormal] = true;
                                                                   else if (chName == "shading tangent")
                                                                       chFlags[(int)SLR::ExtraChannel::ShadingTangent] = true;
                                                                   else if (chName == "distance")
                                                                       chFlags[(int)SLR::ExtraChannel::Distance] = true;
                                                               }
                                                               context.renderingContext->renderer = createUnique<SLR::DebugRenderer>(chFlags);
                                                               return Element();
                                                           }
                                                       };
                                                       return configDebug(config, context, err);
                                                   }
                                                   else {
                                                       *err = ErrorMessage("Unknown method is specified.");
                                                   }
                                                   return Element();
                                               }
                                               );
            stack["setRenderSettings"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{
                                                   {"numThreads", Type::Integer, Element((int32_t)std::thread::hardware_concurrency())},
                                                   {"width", Type::Integer, Element(1024)},
                                                   {"height", Type::Integer, Element(1024)},
                                                   {"timeStart", Type::RealNumber, Element(0.0)},
                                                   {"timeEnd", Type::RealNumber, Element(0.0)},
                                                   {"brightness", Type::RealNumber, Element(1.0f)},
                                                   {"rngSeed", Type::Integer, Element(1509761209)}
                                               },
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   RenderingContext* renderCtx = context.renderingContext;
                                                   renderCtx->width = args.at("width").raw<TypeMap::Integer>();
                                                   renderCtx->height = args.at("height").raw<TypeMap::Integer>();
                                                   renderCtx->timeStart = args.at("timeStart").raw<TypeMap::RealNumber>();
                                                   renderCtx->timeEnd = args.at("timeEnd").raw<TypeMap::RealNumber>();
                                                   renderCtx->brightness = args.at("brightness").raw<TypeMap::RealNumber>();
                                                   renderCtx->rngSeed = args.at("rngSeed").raw<TypeMap::Integer>();
                                                   
                                                   return Element();
                                               }
                                               );
            stack["setEnvironment"] =
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {
                                                       {"path", Type::String}, 
                                                       {"scale", Type::RealNumber, Element(1.0)}
                                                   }, 
                                                   {
                                                       {"image", Type::SpectrumTexture}, 
                                                       {"scale", Type::RealNumber, Element(1.0)}
                                                   }
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       std::string path = context.absFileDirPath + args.at("path").raw<TypeMap::String>();
                                                       float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                                       
                                                       Image2DRef img = createImage2D(path, SLR::ImageStoreMode::AsIs, SLR::SpectrumType::Illuminant, false);
                                                       const Texture2DMappingRef &mapping = Texture2DMapping::sharedInstanceRef();
                                                       SpectrumTextureRef IBLTex = createShared<ImageSpectrumTexture>(mapping, img);
                                                       std::weak_ptr<Scene> sceneWRef = context.scene;
                                                       InfiniteSphereNodeRef infSphere = createShared<InfiniteSphereNode>(sceneWRef, IBLTex, scale);
                                                       
                                                       context.scene->setEnvironmentNode(infSphere);
                                                       
                                                       return Element();
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       SpectrumTextureRef image = args.at("image").rawRef<TypeMap::SpectrumTexture>();
                                                       float scale = args.at("scale").raw<TypeMap::RealNumber>();

                                                       std::weak_ptr<Scene> sceneWRef = context.scene;
                                                       InfiniteSphereNodeRef infSphere = createShared<InfiniteSphereNode>(sceneWRef, image, scale);
                                                       
                                                       context.scene->setEnvironmentNode(infSphere);
                                                       
                                                       return Element();
                                                   }
                                               }
                                               );
        }
        
        SceneParsingDriver parser;
//        parser.traceParsing = true;
        StatementsRef statements = parser.parse(filePath);
        if (!statements) {
            printf("Failed to parse scene file: %s\n", filePath.c_str());
            return false;
        }
        
        for (int i = 0; i < statements->size(); ++i) {
            StatementRef statement = statements->at(i);
            if (!statement->perform(executeContext, &errMsg)) {
                printf("%s\n", errMsg.message.c_str());
                return false;
            }
        }
        return true;
    }

    namespace Spectrum {
        using namespace SLR;
        
#ifdef SLR_Use_Spectral_Representation
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
            return createShared<UpsampledContinuousSpectrum>(spType, space, e0, e1, e2);
        }
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<RegularContinuousSpectrum>(minLambda, maxLambda, values, numSamples);
        }
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<IrregularContinuousSpectrum>(lambdas, values, numSamples);
        }
#else
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
            SLRAssert(e0 >= 0.0 && e1 >= 0.0 && e2 >= 0.0, "Values should not be minus.");
            switch (space) {
                case ColorSpace::sRGB:
                    return createShared<RGBSpectrum>(e0, e1, e2);
                case ColorSpace::sRGB_NonLinear: {
                    e0 = sRGB_degamma(e0);
                    e1 = sRGB_degamma(e1);
                    e2 = sRGB_degamma(e2);
                    return createShared<RGBSpectrum>(e0, e1, e2);
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
                    return createShared<RGBSpectrum>(RGB[0], RGB[1], RGB[2]);
                }
                default:
                    SLRAssert(false, "Invalid color space.");
                    return createShared<RGBSpectrum>();
            }
        }
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<RGBSpectrum>(RegularContinuousSpectrum(minLambda, maxLambda, values, numSamples).convertToRGB(spType));
        }
        SLR_SCENEGRAPH_API AssetSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<RGBSpectrum>(IrregularContinuousSpectrum(lambdas, values, numSamples).convertToRGB(spType));
        }
#endif
    } // namespace Spectrum
}
