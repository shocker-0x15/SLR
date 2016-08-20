//
//  builtin_math.cpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "builtin_math.hpp"
#include <libSLR/BasicTypes/Vector3.h>

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Math {
            const Function min = Function(1, {{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                              float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                              return Element(std::min(x0, x1));
                                          });
            const Function max = Function(1, {{"x0", Type::RealNumber}, {"x1", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x0 = args.at("x0").raw<TypeMap::RealNumber>();
                                              float x1 = args.at("x1").raw<TypeMap::RealNumber>();
                                              return Element(std::max(x0, x1));
                                          });
            const Function clamp = Function(1, {{"x", Type::RealNumber}, {"min", Type::RealNumber}, {"max", Type::RealNumber}},
                                            [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                float x = args.at("x").raw<TypeMap::RealNumber>();
                                                float min = args.at("min").raw<TypeMap::RealNumber>();
                                                float max = args.at("max").raw<TypeMap::RealNumber>();
                                                return Element(std::clamp(x, min, max));
                                            });
            const Function sqrt = Function(1, {{"x", Type::RealNumber}},
                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                               float x = args.at("x").raw<TypeMap::RealNumber>();
                                               return Element(std::sqrt(x));
                                           });
            const Function pow = Function(1, {{"x", Type::RealNumber}, {"e", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x = args.at("x").raw<TypeMap::RealNumber>();
                                              float e = args.at("e").raw<TypeMap::RealNumber>();
                                              return Element(std::pow(x, e));
                                          });
            const Function sin = Function(1, {{"x", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x = args.at("x").raw<TypeMap::RealNumber>();
                                              return Element(std::sin(x));
                                          });
            const Function cos = Function(1, {{"x", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x = args.at("x").raw<TypeMap::RealNumber>();
                                              return Element(std::cos(x));
                                          });
            const Function tan = Function(1, {{"x", Type::RealNumber}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              float x = args.at("x").raw<TypeMap::RealNumber>();
                                              return Element(std::tan(x));
                                          });
            const Function asin = Function(1, {{"x", Type::RealNumber}},
                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                               float x = args.at("x").raw<TypeMap::RealNumber>();
                                               return Element(std::asin(x));
                                           });
            const Function acos = Function(1, {{"x", Type::RealNumber}},
                                           [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                               float x = args.at("x").raw<TypeMap::RealNumber>();
                                               return Element(std::acos(x));
                                           });
            const Function atan = Function(1,
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
                                           });
            const Function dot = Function(1, {{"v0", Type::Vector}, {"v1", Type::Vector}},
                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                              const SLR::Vector3D v0 = args.at("v0").raw<TypeMap::Vector>();
                                              const SLR::Vector3D v1 = args.at("v1").raw<TypeMap::Vector>();
                                              return Element(SLR::dot(v0, v1));
                                          });
            const Function cross = Function(1, {{"v0", Type::Vector}, {"v1", Type::Vector}},
                                            [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                const SLR::Vector3D v0 = args.at("v0").raw<TypeMap::Vector>();
                                                const SLR::Vector3D v1 = args.at("v1").raw<TypeMap::Vector>();
                                                return Element(SLR::cross(v0, v1));
                                            });
        }
    }
}
