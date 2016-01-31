
//
//  API.cpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright c 2015年 渡部 心. All rights reserved.
//

#include "API.hpp"

#include "Parser/SceneParsingDriver.h"

#include <libSLR/Core/Image.h>
#include <libSLR/Core/Transform.h>
#include <libSLR/Core/SurfaceObject.h>
#include <libSLR/Core/XORShift.h>
#include <libSLR/BasicTypes/Spectrum.h>
#include <libSLR/Memory/ArenaAllocator.h>
#include <libSLR/Renderers/PathTracingRenderer.h>

#include "nodes.h"
#include "node_constructor.h"
#include "TriangleMeshNode.h"
#include "camera_nodes.h"
#include "InfiniteSphereNode.h"

#include "Helper/image_loader.h"
#include "textures.hpp"
#include "surface_materials.hpp"

namespace SLRSceneGraph {
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

#ifdef SLR_Defs_Windows
    class FuncGetPathElement {
        std::string pathPrefix;
    public:
        FuncGetPathElement(const std::string &prefix) : pathPrefix(prefix) {}
        Element operator()(const aiString &str) const {
            return Element(pathPrefix + std::string(str.C_Str()));
        }
    };
#endif
    
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
            stack["root"] = Element(TypeMap::Node(), scene->rootNode());
            
            stack["print"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"value", Type::Any}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::cout << args.at("value") << std::endl;
                                 return Element();
                             })
                    );
            stack["addItem"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"tuple", Type::Tuple}, {"key", Type::String, Element(TypeMap::String(), "")}, {"item", Type::Any}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 auto tuple = args.at("tuple").rawRef<TypeMap::Tuple>();
                                 auto key = args.at("key").raw<TypeMap::String>();
                                 tuple->add(key, args.at("item"));
                                 
                                 return Element(TypeMap::Tuple(), tuple);
                             })
                    );
            stack["numElements"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"tuple", Type::Tuple}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 auto tuple = args.at("tuple").rawRef<TypeMap::Tuple>();
                                 return Element((TypeMap::Integer::InternalType)tuple->numParams());
                             })
                    );
            stack["Vector"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 auto x = args.at("x").raw<TypeMap::RealNumber>();
                                 auto y = args.at("y").raw<TypeMap::RealNumber>();
                                 auto z = args.at("z").raw<TypeMap::RealNumber>();
                                 return Element(SLR::Vector3D(x, y, z));
                             })
                    );
            stack["getX"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"point", Type::Point}},
                                 {{"vector", Type::Vector}},
                                 {{"normal", Type::Normal}}
                             },
                             {
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
                                 }})
                    );
            stack["getY"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"point", Type::Point}},
                                 {{"vector", Type::Vector}},
                                 {{"normal", Type::Normal}}
                             },
                             {
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
                                 }})
                    );
            stack["getZ"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"point", Type::Point}},
                                 {{"vector", Type::Vector}},
                                 {{"normal", Type::Normal}}
                             },
                             {
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
                                 }})
                    );
            stack["min"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                 float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                 return Element(std::min(x0, x1));
                             })
                    );
            stack["max"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                 float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                 return Element(std::max(x0, x1));
                             })
                    );
            stack["clamp"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}, {"min", Type::RealNumber}, {"max", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 float min = args.at("min").raw<TypeMap::RealNumber>();
                                 float max = args.at("max").raw<TypeMap::RealNumber>();
                                 return Element(std::clamp(x, min, max));
                             })
                    );
            stack["sqrt"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::sqrt(x));
                             })
                    );
            stack["pow"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}, {"e", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 float e = args.at("e").raw<TypeMap::RealNumber>();
                                 return Element(std::pow(x, e));
                             })
                    );
            stack["sin"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::sin(x));
                             })
                    );
            stack["cos"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::cos(x));
                             })
                    );
            stack["tan"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::tan(x));
                             })
                    );
            stack["asin"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::asin(x));
                             })
                    );
            stack["acos"] =
            Element(TypeMap::Function(),
                    Function(1, {{"x", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float x = args.at("x").raw<TypeMap::RealNumber>();
                                 return Element(std::acos(x));
                             })
                    );
            stack["atan"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"x", Type::RealNumber}},
                                 {{"y", Type::RealNumber}, {"x", Type::RealNumber}}
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     float x = args.at("x").raw<TypeMap::RealNumber>();
                                     return Element(std::atan(x));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     float y = args.at("y").raw<TypeMap::RealNumber>();
                                     float x = args.at("x").raw<TypeMap::RealNumber>();
                                     return Element(std::atan2(y, x));
                                 }
                             })
                    );
            stack["dot"] =
            Element(TypeMap::Function(),
                    Function(1, {{"v0", Type::Vector}, {"v1", Type::Vector}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 const SLR::Vector3D v0 = args.at("v0").raw<TypeMap::Vector>();
                                 const SLR::Vector3D v1 = args.at("v1").raw<TypeMap::Vector>();
                                 return Element(SLR::dot(v0, v1));
                             })
                    );
            stack["cross"] =
            Element(TypeMap::Function(),
                    Function(1, {{"v0", Type::Vector}, {"v1", Type::Vector}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 const SLR::Vector3D v0 = args.at("v0").raw<TypeMap::Vector>();
                                 const SLR::Vector3D v1 = args.at("v1").raw<TypeMap::Vector>();
                                 return Element(SLR::cross(v0, v1));
                             })
                    );
            stack["random"] =
            Element(TypeMap::Function(),
                    Function(1, {},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 static SLR::XORShift rng{2112984105};
                                 return Element(rng.getFloat0cTo1o());
                             })
                    );
            stack["translate"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float tx = args.at("x").raw<TypeMap::RealNumber>();
                                 float ty = args.at("y").raw<TypeMap::RealNumber>();
                                 float tz = args.at("z").raw<TypeMap::RealNumber>();
                                 
                                 return Element(TypeMap::Matrix(), SLR::translate(tx, ty, tz));
                             })
                    );
            stack["rotate"] =
            Element(TypeMap::Function(),
                    Function(1, {{"angle", Type::RealNumber}, {"axis", Type::Vector}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                 SLR::Vector3D axis = args.at("axis").raw<TypeMap::Vector>();
                                 return Element(TypeMap::Matrix(), SLR::rotate(angle, axis));
                             })
                    );
            stack["rotateX"] =
            Element(TypeMap::Function(),
                    Function(1, {{"angle", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                 return Element(TypeMap::Matrix(), SLR::rotateX(angle));
                             })
                    );
            stack["rotateY"] =
            Element(TypeMap::Function(),
                    Function(1, {{"angle", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                 return Element(TypeMap::Matrix(), SLR::rotateY(angle));
                             })
                    );
            stack["rotateZ"] =
            Element(TypeMap::Function(),
                    Function(1, {{"angle", Type::RealNumber}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                 return Element(TypeMap::Matrix(), SLR::rotateZ(angle));
                             })
                    );
            stack["scale"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"s", Type::RealNumber}},
                                 {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}}
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     float s = args.at("s").raw<TypeMap::RealNumber>();
                                     return Element(TypeMap::Matrix(), SLR::scale(s));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     float sx = args.at("x").raw<TypeMap::RealNumber>();
                                     float sy = args.at("y").raw<TypeMap::RealNumber>();
                                     float sz = args.at("z").raw<TypeMap::RealNumber>();
                                     return Element(TypeMap::Matrix(), SLR::scale(sx, sy, sz));
                                 }
                             })
                    );
            stack["createVertex"] =
            Element(TypeMap::Function(),
                    Function(1, {{"position", Type::Tuple}, {"normal", Type::Tuple}, {"tangent", Type::Tuple}, {"texCoord", Type::Tuple}},
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
                                 SLR::Vertex vtx = SLR::Vertex(sigPosition.perform<SLR::Point3D>(procPosition, args.at("position").raw<TypeMap::Tuple>()),
                                                               sigNormal.perform<SLR::Normal3D>(procNormal, args.at("normal").raw<TypeMap::Tuple>()),
                                                               sigTangent.perform<SLR::Tangent3D>(procTangent, args.at("tangent").raw<TypeMap::Tuple>()),
                                                               sigTexCoord.perform<SLR::TexCoord2D>(procTexCoord, args.at("texCoord").raw<TypeMap::Tuple>()));
                                 return Element(TypeMap::Vertex(), vtx);
                             })
                    );
            stack["Spectrum"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {
                                     {"type", Type::String},
                                     {"value", Type::RealNumber}
                                 },
                                 {
                                     {"type", Type::String, Element(TypeMap::String(), "Reflectance")},
                                     {"space", Type::String, Element(TypeMap::String(), "sRGB")},
                                     {"e0", Type::RealNumber},
                                     {"e1", Type::RealNumber},
                                     {"e2", Type::RealNumber}
                                 },
                                 {
                                     {"type", Type::String, Element(TypeMap::String(), "Reflectance")},
                                     {"minWL", Type::RealNumber},
                                     {"maxWL", Type::RealNumber},
                                     {"values", Type::Tuple}
                                 },
                                 {
                                     {"type", Type::String, Element(TypeMap::String(), "Reflectance")},
                                     {"wls", Type::Tuple},
                                     {"values", Type::Tuple}
                                 },
                                 {
                                     {"ID", Type::String},
                                     {"idx", Type::Integer, Element(0)}
                                 }
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     std::string typeStr = args.at("type").raw<TypeMap::String>();
                                     SLR::SpectrumType type;
                                     if (!strToSpectrumType(typeStr, &type)) {
                                         *err = ErrorMessage("Specified spectrum type is invalid.");
                                         return Element();
                                     }
                                     float value = args.at("value").raw<TypeMap::RealNumber>();
                                     
                                     return Element(TypeMap::Spectrum(), Spectrum::create(type, SLR::ColorSpace::sRGB, value, value, value));
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
                                     
                                     return Element(TypeMap::Spectrum(), Spectrum::create(type, space, e0, e1, e2));
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
                                         values.push_back(el.raw<TypeMap::RealNumber>());
                                     }
                                     
                                     return Element(TypeMap::Spectrum(), Spectrum::create(type, minWL, maxWL, values.data(), (uint32_t)numSamples));
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
                                         wls.push_back(elWavelength.raw<TypeMap::RealNumber>());
                                         values.push_back(elValue.raw<TypeMap::RealNumber>());
                                     }
                                     
                                     return Element(TypeMap::Spectrum(), Spectrum::create(type, wls.data(), values.data(), (uint32_t)numSamples));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err)  {
                                     using namespace SLR;
                                     InputSpectrumRef spectrum;
                                     std::string ID = args.at("ID").raw<TypeMap::String>();
                                     auto idx = args.at("idx").raw<TypeMap::Integer>();
                                     if (ID == "D65") {
                                         if (idx != 0) {
                                             *err = ErrorMessage("Specified index is out of range for this spectrum.");
                                             return Element();
                                         }
                                         spectrum = Spectrum::create(SpectrumType::Illuminant, StandardIlluminant::MinWavelength, StandardIlluminant::MaxWavelength,
                                                                     StandardIlluminant::D65, StandardIlluminant::NumSamples);
                                     }
                                     else if (ID == "ColorChecker") {
                                         if (idx < 0 || idx >= 24) {
                                             *err = ErrorMessage("Specified index is out of range for this spectrum.");
                                             return Element();
                                         }
                                         spectrum = Spectrum::create(SpectrumType::Reflectance, ColorChecker::MinWavelength, ColorChecker::MaxWavelength,
                                                                     ColorChecker::Spectra[idx], ColorChecker::NumSamples);
                                     }
                                     else {
                                         if (SpectrumLibrary::IORs.count(ID) > 0) {
                                             if (idx < 0 || idx >= 2) {
                                                 *err = ErrorMessage("Specified index is out of range.");
                                                 return Element();
                                             }
                                             const SpectrumLibrary::IndexOfRefraction &IOR = SpectrumLibrary::IORs.at(ID);
                                             const float* values = idx == 0 ? IOR.etas : IOR.ks;
                                             if (!values) {
                                                 *err = ErrorMessage("This IOR doesn't have the spectrum corresponding to the index specified.");
                                                 return Element();
                                             }
                                             if (IOR.dType == SpectrumLibrary::DistributionType::Regular)
                                                 spectrum = Spectrum::create(SpectrumType::IndexOfRefraction, IOR.minLambdas, IOR.maxLambdas, values, IOR.numSamples);
                                             else
                                                 spectrum = Spectrum::create(SpectrumType::IndexOfRefraction, IOR.lambdas, values, IOR.numSamples);
                                         }
                                         else {
                                             *err = ErrorMessage("unrecognized spectrum ID.");
                                             return Element();
                                         }
                                     }
                                     return Element(TypeMap::Spectrum(), spectrum);
                                 }
                             })
                    );
            stack["Image2D"] =
            Element(TypeMap::Function(),
                    Function(1, {{"path", Type::String}, {"type", Type::String, Element(TypeMap::String(), "Reflectance")}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::string path = args.at("path").raw<TypeMap::String>();
                                 std::string typeStr = args.at("type").raw<TypeMap::String>();
                                 SLR::SpectrumType type;
                                 if (!strToSpectrumType(typeStr, &type)) {
                                     *err = ErrorMessage("Specified spectrum type is invalid.");
                                     return Element();
                                 }
                                 
                                 // TODO: ?? make a memory allocator selectable.
                                 SLR::DefaultAllocator &defMem = SLR::DefaultAllocator::instance();
                                 TiledImage2DRef image = Image::createTiledImage(path.c_str(), &defMem, type);
                                 return Element(TypeMap::Image2D(), image);
                             })
                    );
            stack["SpectrumTexture"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"spectrum", Type::Spectrum}},
                                 {{"image", Type::Image2D}},
                                 {{"procedure", Type::String}, {"params", Type::Tuple}}
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     InputSpectrumRef spectrum = args.at("spectrum").rawRef<TypeMap::Spectrum>();
                                     return Element(TypeMap::SpectrumTexture(), createShared<ConstantSpectrumTexture>(spectrum));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     const auto &image = args.at("image").rawRef<TypeMap::Image2D>();
                                     return Element(TypeMap::SpectrumTexture(), createShared<ImageSpectrumTexture>(image));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     std::string procedure = args.at("procedure").raw<TypeMap::String>();
                                     const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                     if (procedure == "checker board") {
                                         const static Function configFunc{
                                             0, {
                                                 {"c0", Type::Spectrum}, {"c1", Type::Spectrum},
                                             },
                                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                 const InputSpectrumRef c0 = args.at("c0").rawRef<TypeMap::Spectrum>();
                                                 const InputSpectrumRef c1 = args.at("c1").rawRef<TypeMap::Spectrum>();
                                                 return Element(TypeMap::SpectrumTexture(), createShared<CheckerBoardSpectrumTexture>(c0, c1));
                                             }
                                         };
                                         return configFunc(params, context, err);
                                     }
                                     *err = ErrorMessage("Specified procedure is invalid.");
                                     return Element();
                                 }
                             })
                    );
            stack["FloatTexture"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"value", Type::RealNumber}},
                                 {{"image", Type::Image2D}}
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     auto value = args.at("value").raw<TypeMap::RealNumber>();
                                     return Element(TypeMap::FloatTexture(), createShared<ConstantFloatTexture>(value));
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     const auto &image = args.at("image").rawRef<TypeMap::Image2D>();
                                     return Element(TypeMap::FloatTexture(), createShared<ImageFloatTexture>(image));
                                 }
                             })
                    );
            stack["createSurfaceMaterial"] =
            Element(TypeMap::Function(),
                    Function(1, {{"type", Type::String}, {"params", Type::Tuple}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::string type = args.at("type").raw<TypeMap::String>();
                                 const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                 if (type == "matte") {
                                     const static Function configFunc{
                                         0, {
                                             {"reflectance", Type::SpectrumTexture},
                                             {"sigma", Type::FloatTexture, Element(TypeMap::FloatTexture(), nullptr)}
                                         },
                                         [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                             SpectrumTextureRef reflectance = args.at("reflectance").rawRef<TypeMap::SpectrumTexture>();
                                             FloatTextureRef sigma = args.at("sigma").rawRef<TypeMap::FloatTexture>();
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createMatte(reflectance, sigma));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createMetal(coeffR, eta, k));
                                         }
                                     };
                                     return configFunc(params, context, err);
                                 }
                                 else if (type == "glass") {
                                     const static Function configFunc{
                                         0, {
                                             {"coeffR", Type::SpectrumTexture}, {"coeffT", Type::SpectrumTexture},
                                             {"etaExt", Type::SpectrumTexture}, {"etaInt", Type::SpectrumTexture}
                                         },
                                         [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                             SpectrumTextureRef coeffR = args.at("coeffR").rawRef<TypeMap::SpectrumTexture>();
                                             SpectrumTextureRef coeffT = args.at("coeffT").rawRef<TypeMap::SpectrumTexture>();
                                             SpectrumTextureRef etaExt = args.at("etaExt").rawRef<TypeMap::SpectrumTexture>();
                                             SpectrumTextureRef etaInt = args.at("etaInt").rawRef<TypeMap::SpectrumTexture>();
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createGlass(coeffR, coeffT, etaExt, etaInt));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createModifiedWardDur(R, anisoX, anisoY));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createAshikhminShirley(Rd, Rs, nx, ny));
                                         }
                                     };
                                     return configFunc(params, context, err);
                                 }
                                 else if (type == "inverse") {
                                     const static Function configFunc{
                                         0, {
                                             {"base", Type::SurfaceMaterial}
                                         },
                                         [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                             SurfaceMaterialRef base = args.at("base").rawRef<TypeMap::SurfaceMaterial>();
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createInverseMaterial(base));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createEmitterSurfaceMaterial(scatter, emitter));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createMixedMaterial(mat0, mat1, factor));
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
                                             return Element(TypeMap::SurfaceMaterial(), SurfaceMaterial::createSummedMaterial(mat0, mat1));
                                         }
                                     };
                                     return configFunc(params, context, err);
                                 }
                                 *err = ErrorMessage("Specified material type is invalid.");
                                 return Element();
                             })
                    );
            stack["createEmitterSurfaceProperty"] =
            Element(TypeMap::Function(),
                    Function(1, {{"type", Type::String}, {"params", Type::Tuple}},
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
                                             return Element(TypeMap::EmitterSurfaceProperty(), SurfaceMaterial::createDiffuseEmitter(emittance));
                                         }
                                     };
                                     return configFunc(params, context, err);
                                 }
                                 *err = ErrorMessage("Specified material type is invalid.");
                                 return Element();
                             })
                    );
            stack["createMesh"] =
            Element(TypeMap::Function(),
                    Function(1, {{"vertices", Type::Tuple}, {"faces", Type::Tuple}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 const std::vector<Element> &vertices = args.at("vertices").raw<TypeMap::Tuple>().unnamed;
                                 const std::vector<Element> &faces = args.at("faces").raw<TypeMap::Tuple>().unnamed;
                                 
                                 struct Face {
                                     uint32_t v0, v1, v2;
                                     SurfaceMaterialRef mat;
                                 };
                                 static const Function sigFace{
                                     1, {{"indices", Type::Tuple}, {"material", Type::SurfaceMaterial}}
                                 };
                                 static const auto procFace = [](const std::map<std::string, Element> &arg) {
                                     std::map<std::string, Element> indices;
                                     mapParamsToArgs(arg.at("indices").raw<TypeMap::Tuple>(),
                                                     {{"v0", Type::Integer}, {"v1", Type::Integer}, {"v2", Type::Integer}},
                                                     &indices);
                                     uint32_t v0 = indices.at("v0").raw<TypeMap::Integer>();
                                     uint32_t v1 = indices.at("v1").raw<TypeMap::Integer>();
                                     uint32_t v2 = indices.at("v2").raw<TypeMap::Integer>();
                                     return Face{v0, v1, v2, arg.at("material").rawRef<TypeMap::SurfaceMaterial>()};
                                 };
                                 
                                 TriangleMeshNodeRef lightMesh = createShared<TriangleMeshNode>();
                                 for (int i = 0; i < vertices.size(); ++i) {
                                     if (vertices[i].type == Type::Vertex) {
                                         lightMesh->addVertex(vertices[i].raw<TypeMap::Vertex>());
                                     }
                                     else {
                                         const Function &CreateVertex = context.stackVariables["createVertex"].raw<TypeMap::Function>();
                                         Element vtx = CreateVertex(vertices[i].raw<TypeMap::Tuple>(), context, err);
                                         if (err->error)
                                             return Element();
                                         lightMesh->addVertex(vtx.raw<TypeMap::Vertex>());
                                     }
                                 }
                                 for (int i = 0; i < faces.size(); ++i) {
                                     Face f = sigFace.perform<Face>(procFace, faces[i].raw<TypeMap::Tuple>(), Face(), err);
                                     if (err->error)
                                         return Element();
                                     lightMesh->addTriangle(f.v0, f.v1, f.v2, f.mat);
                                 }
                                 
                                 return Element(TypeMap::Mesh(), lightMesh);
                             })
                    );
            stack["createNode"] =
            Element(TypeMap::Function(),
                    Function(1, {},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 return Element(TypeMap::Node(), InternalNode());
                             })
                    );
            stack["createReferenceNode"] =
            Element(TypeMap::Function(),
                    Function(1, {{"node", Type::Node}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 InternalNodeRef node = args.at("node").rawRef<TypeMap::Node>();
                                 ReferenceNodeRef refNode = createShared<ReferenceNode>(node);
                                 return Element(TypeMap::ReferenceNode(), refNode);
                             })
                    );
            stack["setTransform"] =
            Element(TypeMap::Function(),
                    Function(1, {{"node", Type::Node}, {"transform", Type::Transform}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 InternalNodeRef node = args.at("node").rawRef<TypeMap::Node>();
                                 TransformRef tf = args.at("transform").rawRef<TypeMap::Transform>();
                                 node->setTransform(tf);
                                 return Element();
                             })
                    );
            stack["addChild"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {{"parent", Type::Node}, {"child", Type::Node}},
                                 {{"parent", Type::Node}, {"child", Type::ReferenceNode}},
                                 {{"parent", Type::Node}, {"child", Type::Mesh}},
                                 {{"parent", Type::Node}, {"child", Type::Camera}}
                             },
                             {
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     InternalNodeRef parent = args.at("parent").rawRef<TypeMap::Node>();
                                     InternalNodeRef child = args.at("child").rawRef<TypeMap::Node>();
                                     parent->addChildNode(child);
                                     return Element();
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     InternalNodeRef parent = args.at("parent").rawRef<TypeMap::Node>();
                                     ReferenceNodeRef child = args.at("child").rawRef<TypeMap::ReferenceNode>();
                                     parent->addChildNode(child);
                                     return Element();
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     InternalNodeRef parent = args.at("parent").rawRef<TypeMap::Node>();
                                     TriangleMeshNodeRef child = args.at("child").rawRef<TypeMap::Mesh>();
                                     parent->addChildNode(child);
                                     return Element();
                                 },
                                 [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                     InternalNodeRef parent = args.at("parent").rawRef<TypeMap::Node>();
                                     CameraNodeRef child = args.at("child").rawRef<TypeMap::Camera>();
                                     parent->addChildNode(child);
                                     return Element();
                                 }
                             })
                    );
            stack["load3DModel"] =
            Element(TypeMap::Function(),
                    Function(1, {{"path", Type::String}, {"matProc", Type::Function, Element(TypeMap::Function(), nullptr)}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::string path = context.absFileDirPath + args.at("path").raw<TypeMap::String>();
                                 std::string pathPrefix = path.substr(0, path.find_last_of("/") + 1);
                                 auto matProcRef = args.at("matProc").rawRef<TypeMap::Function>();
                                 
                                 CreateMaterialFunction nativeMatProc = createMaterialDefaultFunction;
                                 if (matProcRef) {
                                     const Function &matProc = *matProcRef.get();
                                     nativeMatProc = [&pathPrefix, &matProc, &context, &err](const aiMaterial* aiMat, const std::string &pathPrefix, SLR::Allocator* mem) {
                                         using namespace SLR;
                                         aiString aiStr;
                                         float color[3];
                                         
                                         aiMat->Get(AI_MATKEY_NAME, aiStr);
                                         Element matName = Element(std::string(aiStr.C_Str()));
                                         Element matAttrs = Element(TypeMap::Tuple(), ParameterList());
                                         
                                         auto &attrs = matAttrs.raw<TypeMap::Tuple>();
#ifdef SLR_Defs_Windows
                                         FuncGetPathElement getPathElement{ pathPrefix };
#else
                                         auto getPathElement = [&pathPrefix](const aiString &str) {
                                             return Element(pathPrefix + std::string(str.C_Str()));
                                         };
#endif
                                         auto getRGBElement = [](const float* RGB) {
                                             ParameterListRef values = createShared<ParameterList>();
                                             values->add("", Element(RGB[0]));
                                             values->add("", Element(RGB[1]));
                                             values->add("", Element(RGB[2]));
                                             return Element(TypeMap::Tuple(), values);
                                         };
                                         
                                         Element diffuseTexturesElement = Element(TypeMap::Tuple(), ParameterList());
                                         auto &diffuseTextures = diffuseTexturesElement.raw<TypeMap::Tuple>();
                                         for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_DIFFUSE); ++i) {
                                             if (aiMat->Get(AI_MATKEY_TEXTURE_DIFFUSE(i), aiStr) == aiReturn_SUCCESS)
                                                 diffuseTextures.add("", getPathElement(aiStr));
                                         }
                                         attrs.add("diffuse textures", diffuseTexturesElement);
                                         
                                         Element specularTexturesElement = Element(TypeMap::Tuple(), ParameterList());
                                         auto &specularTextures = specularTexturesElement.raw<TypeMap::Tuple>();
                                         for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_SPECULAR); ++i) {
                                             if (aiMat->Get(AI_MATKEY_TEXTURE_SPECULAR(i), aiStr) == aiReturn_SUCCESS)
                                                 specularTextures.add("", getPathElement(aiStr));
                                         }
                                         attrs.add("specular textures", specularTexturesElement);
                                         
                                         Element emissiveTexturesElement = Element(TypeMap::Tuple(), ParameterList());
                                         auto &emissiveTextures = emissiveTexturesElement.raw<TypeMap::Tuple>();
                                         for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_EMISSIVE); ++i) {
                                             if (aiMat->Get(AI_MATKEY_TEXTURE_EMISSIVE(i), aiStr) == aiReturn_SUCCESS)
                                                 emissiveTextures.add("", getPathElement(aiStr));
                                         }
                                         attrs.add("emissive textures", emissiveTexturesElement);
                                         
                                         Element heightTexturesElement = Element(TypeMap::Tuple(), ParameterList());
                                         auto &heightTextures = heightTexturesElement.raw<TypeMap::Tuple>();
                                         for (int i = 0; i < aiMat->GetTextureCount(aiTextureType_HEIGHT); ++i) {
                                             if (aiMat->Get(AI_MATKEY_TEXTURE_HEIGHT(i), aiStr) == aiReturn_SUCCESS)
                                                 heightTextures.add("", getPathElement(aiStr));
                                         }
                                         attrs.add("height textures", heightTexturesElement);
                                         
                                         Element normalTexturesElement = Element(TypeMap::Tuple(), ParameterList());
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
                                         Element result = matProc(params, context, err);
                                         return result.rawRef<TypeMap::SurfaceMaterial>();
                                     };
                                 }
                                 
                                 InternalNodeRef modelNode;
                                 construct(path, modelNode, nativeMatProc);
                                 if (!modelNode) {
                                     *err = ErrorMessage("Some errors occur during loading a 3D model.");
                                     return Element();
                                 }
                                 modelNode->setName(path);
                                 
                                 return Element(TypeMap::Node(), modelNode);
                             })
                    );
            stack["scanXZFromYPlus"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {"node", Type::Node},
                                 {"numX", Type::Integer},
                                 {"numY", Type::Integer},
                                 {"callback", Type::Function}
                             },
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 using namespace SLR;
                                 InternalNodeRef node = args.at("node").rawRef<TypeMap::Node>();
                                 uint32_t numX = args.at("numX").raw<TypeMap::Integer>();
                                 uint32_t numY = args.at("numY").raw<TypeMap::Integer>();
                                 const Function &callback = args.at("callback").raw<TypeMap::Function>();
                                 
                                 RenderingData buildData;
                                 SLR::ArenaAllocator mem;
                                 node->getRenderingData(mem, nullptr, &buildData);
                                 auto aggregate = createUnique<SurfaceObjectAggregate>(buildData.surfObjs);
                                 
                                 SLR::BoundingBox3D bounds = aggregate->bounds();
                                 for (int i = 0; i < numY; ++i) {
                                     for (int j = 0; j < numX; ++j) {
                                         Ray ray(Point3D(bounds.minP.x + (bounds.maxP.x - bounds.minP.x) * (j + 0.5f) / numX,
                                                         bounds.maxP.y * 1.5f,
                                                         bounds.minP.z + (bounds.maxP.z - bounds.minP.z) * (i + 0.5f) / numY),
                                                 Vector3D(0, -1, 0), 0.0f);
                                         Intersection isect;
                                         aggregate->intersect(ray, &isect);
                                         SurfacePoint surfPt;
                                         isect.getSurfacePoint(&surfPt);
                                         
                                         ParameterList params;
                                         
                                         Element elPosition = Element(surfPt.p);
                                         Element elNormal = Element(Normal3D(surfPt.shadingFrame.z));
                                         Element elTangent = Element(surfPt.shadingFrame.x);
                                         Element elBitangent = Element(surfPt.shadingFrame.y);
                                         params.add("", elPosition);
                                         params.add("", elTangent);
                                         params.add("", elBitangent);
                                         params.add("", elNormal);
                                         
                                         callback(params, context, err);
                                     }
                                 }
                                 
                                 return Element();
                             })
                    );
            stack["createPerspectiveCamera"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
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
                                 
                                 return Element(TypeMap::Camera(), createShared<PerspectiveCameraNode>(sensitivity, aspect, fovY, radius, imgDist, objDist));
                             })
                    );
            stack["setRenderer"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
                                 {"method", Type::String}, {"config", Type::Tuple, Element(TypeMap::Tuple(), ParameterList())}
                             },
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::string method = args.at("method").raw<TypeMap::String>();
                                 const ParameterList &config = args.at("config").raw<TypeMap::Tuple>();
                                 if (method == "PT") {
                                     const static Function configPT{
                                         0, {{"samples", Type::Integer, Element(8)}},
                                         [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                             uint32_t spp = args.at("samples").raw<TypeMap::Integer>();
                                             context.renderingContext->renderer = createUnique<SLR::PathTracingRenderer>(spp);
                                             return Element();
                                         }
                                     };
                                     return configPT(config, context, err);
                                 }
                                 else {
                                     *err = ErrorMessage("Unknown method is specified.");
                                 }
                                 return Element();
                             })
                    );
            stack["setRenderSettings"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {
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
                             })
                    );
            stack["setEnvironment"] =
            Element(TypeMap::Function(),
                    Function(1,
                             {{"path", Type::String}, {"scale", Type::RealNumber, Element(1.0)}},
                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                 std::string path = context.absFileDirPath + args.at("path").raw<TypeMap::String>();
                                 float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                 SLR::DefaultAllocator &defMem = SLR::DefaultAllocator::instance();
                                 
                                 // TODO: make memory allocator selectable.
                                 TiledImage2DRef img = Image::createTiledImage(path, &defMem, SLR::SpectrumType::Illuminant);
                                 SpectrumTextureRef IBLTex = createShared<ImageSpectrumTexture>(img);
                                 std::weak_ptr<Scene> sceneWRef = context.scene;
                                 InfiniteSphereNodeRef infSphere = createShared<InfiniteSphereNode>(sceneWRef, IBLTex, scale);
                                 
                                 context.scene->setEnvNode(infSphere);
                                 
                                 return Element();
                             })
                    );
        }
        
        SceneParsingDriver parser;
//        parser.traceParsing = true;
        StatementsRef statements = parser.parse(filePath);
        if (!statements) {
            printf("Failed to parse scene file: %s", filePath.c_str());
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
        
#ifdef Use_Spectral_Representation
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
            return createShared<UpsampledContinuousSpectrum>(spType, space, e0, e1, e2);
        }
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
            return createShared<RegularContinuousSpectrum>(minLambda, maxLambda, values, numSamples);
        }
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
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
        
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2) {
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
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples) {
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
        SLR_SCENEGRAPH_API InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples) {
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
        
        SLR_SCENEGRAPH_API std::shared_ptr<SLR::TiledImage2D> createTiledImage(const std::string &filepath, SLR::Allocator *mem, SLR::SpectrumType spType, bool gammaCorrection) {
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
