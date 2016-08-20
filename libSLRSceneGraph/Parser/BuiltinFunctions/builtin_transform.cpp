//
//  builtin_transform.cpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "builtin_transform.hpp"
#include <libSLR/Core/Transform.h>

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Transform {
            const Function translate = Function(1,
                                                {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}},
                                                [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                    float tx = args.at("x").raw<TypeMap::RealNumber>();
                                                    float ty = args.at("y").raw<TypeMap::RealNumber>();
                                                    float tz = args.at("z").raw<TypeMap::RealNumber>();
                                                    
                                                    return Element(TypeMap::Matrix(), SLR::translate(tx, ty, tz));
                                                });
            const Function rotate = Function(1, {{"angle", Type::RealNumber}, {"axis", Type::Vector}},
                                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                 float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                                 SLR::Vector3D axis = args.at("axis").raw<TypeMap::Vector>();
                                                 return Element(TypeMap::Matrix(), SLR::rotate(angle, axis));
                                             });
            const Function rotateX = Function(1, {{"angle", Type::RealNumber}},
                                              [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                  float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                                  return Element(TypeMap::Matrix(), SLR::rotateX(angle));
                                              });
            const Function rotateY = Function(1, {{"angle", Type::RealNumber}},
                                              [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                  float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                                  return Element(TypeMap::Matrix(), SLR::rotateY(angle));
                                              });
            const Function rotateZ = Function(1, {{"angle", Type::RealNumber}},
                                              [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                  float angle = args.at("angle").raw<TypeMap::RealNumber>();
                                                  return Element(TypeMap::Matrix(), SLR::rotateZ(angle));
                                              });
            const Function scale = Function(1,
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
                                            });
            const Function lookAt = Function(1, {{"eye", Type::Tuple}, {"target", Type::Tuple}, {"up", Type::Tuple}},
                                             [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                 static const Function sigEye = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                 static const Function sigTarget = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                 static const Function sigUp = Function(1, {{"x", Type::RealNumber}, {"y", Type::RealNumber}, {"z", Type::RealNumber}});
                                                 static const auto procEye = [](const std::map<std::string, Element> &arg) {
                                                     return SLR::Point3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                 };
                                                 static const auto procTarget = [](const std::map<std::string, Element> &arg) {
                                                     return SLR::Point3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                 };
                                                 static const auto procUp = [](const std::map<std::string, Element> &arg) {
                                                     return SLR::Vector3D(arg.at("x").raw<TypeMap::RealNumber>(), arg.at("y").raw<TypeMap::RealNumber>(), arg.at("z").raw<TypeMap::RealNumber>());
                                                 };
                                                 SLR::Matrix4x4 matRawLookAt = SLR::lookAt(sigEye.perform<SLR::Point3D>(procEye, args.at("eye").raw<TypeMap::Tuple>()),
                                                                                           sigTarget.perform<SLR::Point3D>(procTarget, args.at("target").raw<TypeMap::Tuple>()),
                                                                                           sigUp.perform<SLR::Vector3D>(procUp, args.at("up").raw<TypeMap::Tuple>()));
                                                 SLR::Matrix4x4 mat = SLR::invert(matRawLookAt) * SLR::rotateY((float)M_PI);
                                                 return Element(TypeMap::Matrix(), mat);
                                             });
            const Function AnimatedTransform = Function(1, {{"tfStart", Type::Matrix}, {"tfEnd", Type::Matrix}, {"tBegin", Type::RealNumber}, {"tEnd", Type::RealNumber}},
                                                        [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                            const SLR::Matrix4x4 &tfStart = args.at("tfStart").raw<TypeMap::Matrix>();
                                                            const SLR::Matrix4x4 &tfEnd = args.at("tfEnd").raw<TypeMap::Matrix>();
                                                            float tBegin = args.at("tBegin").raw<TypeMap::RealNumber>();
                                                            float tEnd = args.at("tEnd").raw<TypeMap::RealNumber>();
                                                            TransformRef tf = createShared<SLR::AnimatedTransform>(tfStart, tfEnd, tBegin, tEnd);
                                                            return Element(TypeMap::Transform(), tf);
                                                        });
        }
    }
}

