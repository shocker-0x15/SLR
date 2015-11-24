//
//  API.hpp
//
//  Created by 渡部 心 on 2015/10/06.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#ifndef SLRSceneGraph_API_hpp
#define SLRSceneGraph_API_hpp

#include <libSLR/defines.h>
#include "references.h"

#include <libSLR/BasicTypes/Spectrum.h>

enum class Type : uint32_t {
    Integer = 0,
    RealNumber,
    String,
    Matrix,
    Transform, 
    Spectrum,
    SpectrumTexture,
    NormalTexture,
    FloatTexture,
    SurfaceMaterial,
    EmitterSurfaceProperty,
    Mesh, 
    Sensor,
    Camera,
    Node,
    Tuple, 
    Void,
    NumTypes
};

struct Element {
    Type type;
    std::shared_ptr<void> valueRef;
    Element(Type t = Type::Void, const std::shared_ptr<void> &vr = nullptr) : type(t), valueRef(vr) { };
    
    operator bool() const { return type != Type::Void; };
    template <typename T>
    const T& as() const { return *(T*)valueRef.get(); };
    template <typename T>
    std::shared_ptr<T> asRef() const { return std::static_pointer_cast<T>(valueRef); };
    bool isConvertibleTo(Type toType) const;
    Element convertTo(Type toType) const;
};

std::ostream &operator<<(std::ostream &out, const Element &elem);

struct TypeInfo {
    typedef Element (*convertFunction)(const Element &);
    std::map<Type, convertFunction> convertFunctions;
    
    bool isConvertibleTo(Type toType) const;

    static bool initialized;
    static TypeInfo infos[(uint32_t)Type::NumTypes];
    static void init();
};

struct Parameter {
    std::string name;
    Element elem;
    Parameter() : name("") { };
    Parameter(const std::string &n, const Element &e) : name(n), elem(e) { };
};

struct ParameterList {
    std::map<std::string, Element> named;
    std::vector<Element> unnamed;
    
    bool add(const Parameter &param) {
        if (param.name.empty())
            unnamed.push_back(param.elem);
        else {
            if (named.count(param.name))
                return false;
                
            named[param.name] = param.elem;
        }
        return true;
    };
    
    uint32_t numParams() const {
        return uint32_t(named.size() + unnamed.size());
    };
    
    Element operator()(const std::string &key, Type expectedType) const {
        if (named.count(key) && named.at(key).isConvertibleTo(expectedType))
            return named.at(key).convertTo(expectedType);
        else
            return Element();
    };
    
    Element operator()(uint32_t idx, Type expectedType) const {
        if (idx < unnamed.size() && unnamed[idx].isConvertibleTo(expectedType))
            return unnamed[idx].convertTo(expectedType);
        else
            return Element();
    };
};

typedef std::shared_ptr<ParameterList> ParameterListRef;

struct ErrorMessage {
    bool error;
    std::string message;
    ErrorMessage() : error(false) { };
    ErrorMessage(const std::string &msg) : error(true), message(msg) { };
};

namespace SLRSceneGraph {
    bool readScene(const std::string &filePath, Scene* scene);
    
    Element mulMatrix4x4(const Element &lm, const Element &rm);
    
    Element Translate(const ParameterList &params, ErrorMessage* err);
    Element RotateX(const ParameterList &params, ErrorMessage* err);
    Element RotateY(const ParameterList &params, ErrorMessage* err);
    Element RotateZ(const ParameterList &params, ErrorMessage* err);
    Element Scale(const ParameterList &params, ErrorMessage* err);
    
    Element CreateSpectrum(const ParameterList &params, ErrorMessage* err);
    Element CreateSpectrumTexture(const ParameterList &params, ErrorMessage* err);
    Element CreateMatte(const ParameterList &params, ErrorMessage* err);
    Element CreateDiffuseEmitter(const ParameterList &params, ErrorMessage* err);
    Element CreateEmitterSurfaceMaterial(const ParameterList &params, ErrorMessage* err);
    
    Element CreateMesh(const ParameterList &params, ErrorMessage* err);
    Element CreateNode(const ParameterList &params, ErrorMessage* err);
    Element SetTransform(const ParameterList &params, ErrorMessage* err);
    Element AddChild(const ParameterList &params, ErrorMessage* err);
    Element Load3DModel(const ParameterList &params, ErrorMessage* err);
    
    Element SetRenderer(const ParameterList &params, RenderingContext* context, ErrorMessage* err);
    
    namespace Spectrum {
        using namespace SLR;
        
        InputSpectrumRef create(SpectrumType spType, ColorSpace space, SpectrumFloat e0, SpectrumFloat e1, SpectrumFloat e2);
        InputSpectrumRef create(SpectrumType spType, SpectrumFloat minLambda, SpectrumFloat maxLambda, const SpectrumFloat* values, uint32_t numSamples);
        InputSpectrumRef create(SpectrumType spType, const SpectrumFloat* lambdas, const SpectrumFloat* values, uint32_t numSamples);
    }
    
    namespace Image {
        extern std::map<std::string, Image2DRef> s_imageDB;
        
        std::shared_ptr<SLR::TiledImage2D> createTiledImage(const std::string &filepath, SLR::Allocator *mem, SLR::SpectrumType spType, bool gammaCorrection = false);
    }
}

#endif /* API_hpp */
