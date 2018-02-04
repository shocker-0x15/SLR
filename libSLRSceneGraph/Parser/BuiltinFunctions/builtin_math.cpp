//
//  builtin_math.cpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "builtin_math.h"

#include <libSLR/BasicTypes/Point3D.h>
#include <libSLR/BasicTypes/Vector3D.h>

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Math {
            const Element abs = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::fabs(x));
                                               });
            const Element min = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                                   float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                                   return Element(std::min(x0, x1));
                                               });
            const Element max = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                                   float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                                   return Element(std::max(x0, x1));
                                               });
            const Element clamp = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}, {"min", Type::RealNumber}, {"max", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   float min = args.at("min").raw<TypeMap::RealNumber>();
                                                   float max = args.at("max").raw<TypeMap::RealNumber>();
                                                   return Element(std::clamp(x, min, max));
                                               });
            const Element sqrt = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::sqrt(x));
                                               });
            const Element pow = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}, {"e", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   float e = args.at("e").raw<TypeMap::RealNumber>();
                                                   return Element(std::pow(x, e));
                                               });
            const Element exp = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::exp(x));
                                               });
            const Element ln = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::log(x));
                                               });
            const Element log2 = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::log2(x));
                                               });
            const Element log10 = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::log10(x));
                                               });
            const Element sin = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::sin(x));
                                               });
            const Element cos = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::cos(x));
                                               });
            const Element tan = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::tan(x));
                                               });
            const Element asin = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::asin(x));
                                               });
            const Element acos = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<ArgInfo>{{"x", Type::RealNumber}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   float x = args.at("x").raw<TypeMap::RealNumber>();
                                                   return Element(std::acos(x));
                                               });
            const Element atan = 
            Element::create<TypeMap::Function>(1,
                                               std::vector<std::vector<ArgInfo>>{
                                                   {{"x", Type::RealNumber}},
                                                   {{"y", Type::RealNumber}, {"x", Type::RealNumber}}
                                               },
                                               std::vector<Function::Procedure>{
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       float x = args.at("x").raw<TypeMap::RealNumber>();
                                                       return Element(std::atan(x));
                                                   },
                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                       float y = args.at("y").raw<TypeMap::RealNumber>();
                                                       float x = args.at("x").raw<TypeMap::RealNumber>();
                                                       return Element(std::atan2(y, x));
                                                   }
                                               });
            const Element normalize = 
            Element::create<TypeMap::Function>(1, 
                                               std::vector<ArgInfo>{{"v", Type::Vector}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const SLR::Vector3D &v = args.at("v").raw<TypeMap::Vector>();
                                                   return Element(SLR::normalize(v));
                                               });
            const Element dot = 
            Element::create<TypeMap::Function>(1, 
                                               std::vector<ArgInfo>{{"v0", Type::Vector}, {"v1", Type::Vector}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const SLR::Vector3D &v0 = args.at("v0").raw<TypeMap::Vector>();
                                                   const SLR::Vector3D &v1 = args.at("v1").raw<TypeMap::Vector>();
                                                   return Element(SLR::dot(v0, v1));
                                               });
            const Element cross = 
            Element::create<TypeMap::Function>(1, 
                                               std::vector<ArgInfo>{{"v0", Type::Vector}, {"v1", Type::Vector}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const SLR::Vector3D &v0 = args.at("v0").raw<TypeMap::Vector>();
                                                   const SLR::Vector3D &v1 = args.at("v1").raw<TypeMap::Vector>();
                                                   return Element(SLR::cross(v0, v1));
                                               });
            const Element distance = 
            Element::create<TypeMap::Function>(1, 
                                               std::vector<ArgInfo>{{"p0", Type::Point}, {"p1", Type::Point}},
                                               [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                   const SLR::Point3D &p0 = args.at("p0").raw<TypeMap::Point>();
                                                   const SLR::Point3D &p1 = args.at("p1").raw<TypeMap::Point>();
                                                   return Element(SLR::distance(p0, p1));
                                               });
        }
    }
}
