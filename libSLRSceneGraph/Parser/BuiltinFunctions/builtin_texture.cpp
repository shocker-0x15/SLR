//
//  builtin_texture.cpp
//
//  Created by 渡部 心 on 2016/08/20.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "builtin_texture.hpp"
#include "textures.hpp"

namespace SLRSceneGraph {
    namespace BuiltinFunctions {
        namespace Texture {
            const Function Texture2DMapping = Function(1,
                                                       {
                                                           {"type", Type::String, Element("texcoord 2D")}, {"params", Type::Tuple, Element(TypeMap::Tuple(), ParameterList())}
                                                       },
                                                       [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                           std::string type = args.at("type").raw<TypeMap::String>();
                                                           const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                           if (type == "texcoord 2D") {
                                                               return Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef());
                                                           }
                                                           *err = ErrorMessage("Specified type is invalid.");
                                                           return Element();
                                                       });
            const Function Texture3DMapping = Function(1,
                                                       {
                                                           {"type", Type::String, Element("texcoord 2D")}, {"params", Type::Tuple, Element(TypeMap::Tuple(), ParameterList())}
                                                       },
                                                       [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                           std::string type = args.at("type").raw<TypeMap::String>();
                                                           const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                           if (type == "texcoord 2D") {
                                                               return Element(TypeMap::Texture3DMapping(), Texture3DMapping::sharedInstanceRef());
                                                           }
                                                           else if (type == "world pos") {
                                                               return Element(TypeMap::Texture3DMapping(), WorldPosition3DMapping::sharedInstanceRef());
                                                           }
                                                           *err = ErrorMessage("Specified type is invalid.");
                                                           return Element();
                                                       });
            const Function SpectrumTexture = Function(1,
                                                      {
                                                          {{"spectrum", Type::Spectrum}},
                                                          {{"image", Type::Image2D}, {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}},
                                                          {{"procedure", Type::String}, {"params", Type::Tuple}}
                                                      },
                                                      {
                                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                              InputSpectrumRef spectrum = args.at("spectrum").rawRef<TypeMap::Spectrum>();
                                                              return Element(TypeMap::SpectrumTexture(), createShared<ConstantSpectrumTexture>(spectrum));
                                                          },
                                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                              const auto &image = args.at("image").rawRef<TypeMap::Image2D>();
                                                              const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                              return Element(TypeMap::SpectrumTexture(), createShared<ImageSpectrumTexture>(mapping, image));
                                                          },
                                                          [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                              std::string procedure = args.at("procedure").raw<TypeMap::String>();
                                                              const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                              if (procedure == "checker board") {
                                                                  const static Function configFunc{
                                                                      0, {
                                                                          {"c0", Type::Spectrum}, {"c1", Type::Spectrum}, {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}
                                                                      },
                                                                      [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                          const InputSpectrumRef c0 = args.at("c0").rawRef<TypeMap::Spectrum>();
                                                                          const InputSpectrumRef c1 = args.at("c1").rawRef<TypeMap::Spectrum>();
                                                                          const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                                          return Element(TypeMap::SpectrumTexture(), createShared<CheckerBoardSpectrumTexture>(mapping, c0, c1));
                                                                      }
                                                                  };
                                                                  return configFunc(params, context, err);
                                                              }
                                                              else if (procedure == "voronoi") {
                                                                  const static Function configFunc{
                                                                      0, {
                                                                          {"scale", Type::RealNumber}, {"brightness", Type::RealNumber, Element(TypeMap::RealNumber(), 0.8f)}, {"mapping", Type::Texture3DMapping, Element(TypeMap::Texture3DMapping(), WorldPosition3DMapping::sharedInstanceRef())}
                                                                      },
                                                                      [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                          float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                                                          float brightness = args.at("brightness").raw<TypeMap::RealNumber>();
                                                                          const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture3DMapping>();
                                                                          return Element(TypeMap::SpectrumTexture(), createShared<VoronoiSpectrumTexture>(mapping, scale, brightness));
                                                                      }
                                                                  };
                                                                  return configFunc(params, context, err);
                                                              }
                                                              *err = ErrorMessage("Specified procedure is invalid.");
                                                              return Element();
                                                          }
                                                      });
            const Function NormalTexture = Function(1,
                                                    {
                                                        {{"image", Type::Image2D}, {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}},
                                                        {{"procedure", Type::String}, {"params", Type::Tuple}}
                                                    },
                                                    {
                                                        [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                            const auto &image = args.at("image").rawRef<TypeMap::Image2D>();
                                                            const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                            return Element(TypeMap::NormalTexture(), createShared<ImageNormal3DTexture>(mapping, image));
                                                        },
                                                        [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                            std::string procedure = args.at("procedure").raw<TypeMap::String>();
                                                            const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                            if (procedure == "checker board") {
                                                                const static Function configFunc{
                                                                    0, {
                                                                        {"stepWidth", Type::RealNumber, Element(0.05)},
                                                                        {"reverse", Type::Bool, Element(false)},
                                                                        {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}
                                                                    },
                                                                    [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                        float stepWidth = args.at("stepWidth").raw<TypeMap::RealNumber>();
                                                                        bool reverse = args.at("reverse").raw<TypeMap::Bool>();
                                                                        const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                                        return Element(TypeMap::NormalTexture(), createShared<CheckerBoardNormal3DTexture>(mapping, stepWidth, reverse));
                                                                    }
                                                                };
                                                                return configFunc(params, context, err);
                                                            }
                                                            else if (procedure == "voronoi") {
                                                                const static Function configFunc{
                                                                    0, {
                                                                        {"scale", Type::RealNumber},
                                                                        {"thetaMax", Type::RealNumber, Element(M_PI / 6)},
                                                                        {"mapping", Type::Texture3DMapping, Element(TypeMap::Texture3DMapping(), WorldPosition3DMapping::sharedInstanceRef())}
                                                                    },
                                                                    [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                        float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                                                        float thetaMax = args.at("thetaMax").raw<TypeMap::RealNumber>();
                                                                        const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture3DMapping>();
                                                                        return Element(TypeMap::NormalTexture(), createShared<VoronoiNormal3DTexture>(mapping, scale, thetaMax));
                                                                    }
                                                                };
                                                                return configFunc(params, context, err);
                                                            }
                                                            *err = ErrorMessage("Specified procedure is invalid.");
                                                            return Element();
                                                        }
                                                    });
            const Function FloatTexture = Function(1,
                                                   {
                                                       {{"value", Type::RealNumber}},
                                                       {{"image", Type::Image2D}, {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}},
                                                       {{"procedure", Type::String}, {"params", Type::Tuple}}
                                                   },
                                                   {
                                                       [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                           auto value = args.at("value").raw<TypeMap::RealNumber>();
                                                           return Element(TypeMap::FloatTexture(), createShared<ConstantFloatTexture>(value));
                                                       },
                                                       [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                           const auto &image = args.at("image").rawRef<TypeMap::Image2D>();
                                                           const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                           return Element(TypeMap::FloatTexture(), createShared<ImageFloatTexture>(mapping, image));
                                                       },
                                                       [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                           std::string procedure = args.at("procedure").raw<TypeMap::String>();
                                                           const ParameterList &params = args.at("params").raw<TypeMap::Tuple>();
                                                           if (procedure == "checker board") {
                                                               const static Function configFunc{
                                                                   0, {
                                                                       {"c0", Type::RealNumber}, {"c1", Type::RealNumber}, {"mapping", Type::Texture2DMapping, Element(TypeMap::Texture2DMapping(), Texture2DMapping::sharedInstanceRef())}
                                                                   },
                                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                       float c0 = args.at("c0").raw<TypeMap::RealNumber>();
                                                                       float c1 = args.at("c1").raw<TypeMap::RealNumber>();
                                                                       const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture2DMapping>();
                                                                       return Element(TypeMap::FloatTexture(), createShared<CheckerBoardFloatTexture>(mapping, c0, c1));
                                                                   }
                                                               };
                                                               return configFunc(params, context, err);
                                                           }
                                                           else if (procedure == "voronoi") {
                                                               const static Function configFunc{
                                                                   0, {
                                                                       {"scale", Type::RealNumber},
                                                                       {"valueScale", Type::RealNumber, Element(1.0)},
                                                                       {"flat", Type::Bool, Element(true)},
                                                                       {"mapping", Type::Texture3DMapping, Element(TypeMap::Texture3DMapping(), WorldPosition3DMapping::sharedInstanceRef())}
                                                                   },
                                                                   [](const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) {
                                                                       float scale = args.at("scale").raw<TypeMap::RealNumber>();
                                                                       float valueScale = args.at("valueScale").raw<TypeMap::RealNumber>();
                                                                       bool flat = args.at("flat").raw<TypeMap::Bool>();
                                                                       const auto &mapping = args.at("mapping").rawRef<TypeMap::Texture3DMapping>();
                                                                       return Element(TypeMap::FloatTexture(), createShared<VoronoiFloatTexture>(mapping, scale, valueScale, flat));
                                                                   }
                                                               };
                                                               return configFunc(params, context, err);
                                                           }
                                                           *err = ErrorMessage("Specified procedure is invalid.");
                                                           return Element();
                                                       }
                                                   });
        }
    }
}