//
//  SceneParser.cpp
//
//  Created by 渡部 心 on 2015/12/16.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "SceneParser.hpp"
#include "../API.hpp"

#include <libSLR/Core/Transform.h>

namespace SLRSceneGraph {
    std::ostream &operator<<(std::ostream &out, Type t) {
        switch (t) {
            case Type::Bool:
                out << "Bool";
                break;
            case Type::Integer:
                out << "Integer";
                break;
            case Type::RealNumber:
                out << "RealNumber";
                break;
            case Type::String:
                out << "String";
                break;
            case Type::Matrix:
                out << "Matrix";
                break;
            case Type::Vertex:
                out << "Vertex";
                break;
            case Type::Transform:
                out << "Transform";
                break;
            case Type::Spectrum:
                out << "Spectrum";
                break;
            case Type::SpectrumTexture:
                out << "SpectrumTexture";
                break;
            case Type::NormalTexture:
                out << "NormalTexture";
                break;
            case Type::FloatTexture:
                out << "FloatTexture";
                break;
            case Type::SurfaceMaterial:
                out << "SurfaceMaterial";
                break;
            case Type::EmitterSurfaceProperty:
                out << "EmitterSurfaceProperty";
                break;
            case Type::Mesh:
                out << "Mesh";
                break;
            case Type::Camera:
                out << "Camera";
                break;
            case Type::Node:
                out << "Node";
                break;
            case Type::ReferenceNode:
                out << "ReferenceNode";
                break;
            case Type::Tuple:
                out << "Tuple";
                break;
            case Type::Function:
                out << "Function";
                break;
            case Type::Any:
                out << "Any";
                break;
            case Type::Void:
                out << "Void";
                break;
            case Type::Error:
                out << "Error";
                break;
            case Type::NumTypes:
                out << "NumTypes";
                break;
            default:
                break;
        }
        return out;
    }
    
    std::ostream &operator<<(std::ostream &out, const Element &elem) {
        switch (elem.type) {
            case Type::Bool:
                out << elem.raw<TypeMap::Bool>();
                break;
            case Type::Integer:
                out << elem.raw<TypeMap::Integer>();
                break;
            case Type::RealNumber:
                out << elem.raw<TypeMap::RealNumber>();
                break;
            case Type::String:
                out << "\"" << elem.raw<TypeMap::String>() << "\"";
                break;
            case Type::Matrix:
                out << "Matrix";
                break;
            case Type::Vertex:
                out << "Vertex";
                break;
            case Type::Transform:
                out << "Transform";
                break;
            case Type::Spectrum:
                out << "Spectrum";
                break;
            case Type::SpectrumTexture:
                out << "SpectrumTexture";
                break;
            case Type::NormalTexture:
                out << "NormalTexture";
                break;
            case Type::FloatTexture:
                out << "FloatTexture";
                break;
            case Type::SurfaceMaterial:
                out << "SurfaceMaterial";
                break;
            case Type::EmitterSurfaceProperty:
                out << "EmitterSurfaceProperty";
                break;
            case Type::Mesh:
                out << "Mesh";
                break;
            case Type::Camera:
                out << "Camera";
                break;
            case Type::Node:
                out << "Node";
                break;
            case Type::ReferenceNode:
                out << "ReferenceNode";
                break;
            case Type::Tuple:
                out << elem.raw<TypeMap::Tuple>();
                break;
            case Type::Function:
                out << "Function";
                break;
            case Type::Any:
                out << "Any";
                break;
            case Type::Void:
                out << "Void";
                break;
            case Type::Error:
                out << "Error";
                break;
            case Type::NumTypes:
                out << "NumTypes";
                break;
            default:
                out << "Unknown";
                break;
        }
        return out;
    }
    
    std::ostream &operator<<(std::ostream &out, const ParameterList &params) {
        out << "{";
        for (auto it = params.named.begin(); it != params.named.end(); ++it) {
            out << "\"" << it->first << "\"" << ": " << it->second;
            if (std::distance(params.named.begin(), it) + 1 < params.named.size())
                out << ", ";
        }
        out << "}, ";
        out << "(";
        for (auto it = params.unnamed.begin(); it != params.unnamed.end(); ++it) {
            out << *it;
            if (std::distance(params.unnamed.begin(), it) + 1 < params.unnamed.size())
                out << ", ";
        }
        out << ")";
        return out;
    }
    
    Element::Element(const SLR::Point3D &val) : type(Type::Point), valueRef(createShared<SLR::Point3D>(val)) { }
    Element::Element(const SLR::Vector3D &val) : type(Type::Vector), valueRef(createShared<SLR::Vector3D>(val)) { }
    Element::Element(const SLR::Normal3D &val) : type(Type::Normal), valueRef(createShared<SLR::Normal3D>(val)) { }
    
    bool Element::isConvertibleTo(Type toType) const {
        return TypeInfo::infos[(uint32_t)type].isConvertibleTo(toType);
    };
    
    Element Element::as(Type toType) const {
        SLRAssert(isConvertibleTo(toType), "Specified type is invalid.");
        TypeInfo::convertFunction func = TypeInfo::infos[(uint32_t)type].convertFunctions.at(toType);
        return func(*this);
    }
    
    Element Element::operator++(int) {
        const TypeInfo &thisType = TypeInfo::infos[(uint32_t)type];
        if (!thisType.postIncOperator)
            return Element(ErrorMessage("This type does not have the postfix ++ operator definition."));
        return thisType.postIncOperator(*this);
    }
    
    Element Element::operator--(int) {
        const TypeInfo &thisType = TypeInfo::infos[(uint32_t)type];
        if (!thisType.postDecOperator)
            return Element(ErrorMessage("This type does not have the postfix -- operator definition."));
        return thisType.postDecOperator(*this);
    }
    
    Element &Element::operator++() {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.preIncOperator)
            leftType.preIncOperator(*this);
        else
            *this = Element(ErrorMessage("Left type does not have the prefix ++ operator definition."));
        return *this;
    }
    
    Element &Element::operator--() {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.preDecOperator)
            leftType.preDecOperator(*this);
        else
            *this = Element(ErrorMessage("Left type does not have the prefix -- operator definition."));
        return *this;
    }
    
    Element Element::operator+() const {
        const TypeInfo &thisType = TypeInfo::infos[(uint32_t)type];
        if (!thisType.affOperator)
            return Element(ErrorMessage("This type does not have the unary + operator definition."));
        return thisType.affOperator(*this);
    }
    
    Element Element::operator-() const {
        const TypeInfo &thisType = TypeInfo::infos[(uint32_t)type];
        if (!thisType.negOperator)
            return Element(ErrorMessage("This type does not have the unary - operator definition."));
        return thisType.negOperator(*this);
    }
    
    Element Element::operator!() const {
        const TypeInfo &thisType = TypeInfo::infos[(uint32_t)type];
        if (!thisType.logicNotOperator)
            return Element(ErrorMessage("This type does not have the unary ! operator definition."));
        return thisType.logicNotOperator(*this);
    }
    
    Element Element::operator*(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.mulOperator)
            return Element(ErrorMessage("Left type does not have the * operator definition."));
        return leftType.mulOperator(*this, rElem);
    }
    
    Element Element::operator/(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.divOperator)
            return Element(ErrorMessage("Left type does not have the / operator definition."));
        return leftType.divOperator(*this, rElem);
    }
    
    Element Element::operator%(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.remOperator)
            return Element(ErrorMessage("Left type does not have the \% operator definition."));
        return leftType.remOperator(*this, rElem);
    }
    
    Element Element::operator+(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.addOperator)
            return Element(ErrorMessage("Left type does not have the + operator definition."));
        return leftType.addOperator(*this, rElem);
    }
    
    Element Element::operator-(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.subOperator)
            return Element(ErrorMessage("Left type does not have the - operator definition."));
        return leftType.subOperator(*this, rElem);
    }
    
    Element Element::operator<(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.lessOperator)
            return Element(ErrorMessage("Left type does not have the < operator definition."));
        return leftType.lessOperator(*this, rElem);
    }
    
    Element Element::operator>(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.greaterOperator)
            return Element(ErrorMessage("Left type does not have the > operator definition."));
        return leftType.greaterOperator(*this, rElem);
    }
    
    Element Element::operator<=(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.lessEqOperator)
            return Element(ErrorMessage("Left type does not have the <= operator definition."));
        return leftType.lessEqOperator(*this, rElem);
    }
    
    Element Element::operator>=(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.greaterEqOperator)
            return Element(ErrorMessage("Left type does not have the >= operator definition."));
        return leftType.greaterEqOperator(*this, rElem);
    }
    
    Element Element::operator==(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.eqOperator)
            return Element(ErrorMessage("Left type does not have the == operator definition."));
        return leftType.eqOperator(*this, rElem);
    }
    
    Element Element::operator!=(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.neqOperator)
            return Element(ErrorMessage("Left type does not have the != operator definition."));
        return leftType.neqOperator(*this, rElem);
    }
    
    Element Element::operator&&(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.logicAndOperator)
            return Element(ErrorMessage("Left type does not have the && operator definition."));
        return leftType.logicAndOperator(*this, rElem);
    }
    
    Element Element::operator||(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.logicOrOperator)
            return Element(ErrorMessage("Left type does not have the || operator definition."));
        return leftType.logicOrOperator(*this, rElem);
    }
    
    Element &Element::substitute(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.substOperator)
            leftType.substOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the = operator definition."));
        return *this;
    }
    
    Element &Element::operator+=(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.addSubstOperator)
            leftType.addSubstOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the += operator definition."));
        return *this;
    }
    
    Element &Element::operator-=(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.subSubstOperator)
            leftType.subSubstOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the -= operator definition."));
        return *this;
    }
    
    Element &Element::operator*=(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.mulSubstOperator)
            leftType.mulSubstOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the *= operator definition."));
        return *this;
    }
    
    Element &Element::operator/=(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.divSubstOperator)
            leftType.divSubstOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the /= operator definition."));
        return *this;
    }
    
    Element &Element::operator%=(const Element &rElem) {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (leftType.remSubstOperator)
            leftType.remSubstOperator(*this, rElem);
        else
            *this = Element(ErrorMessage("Left type does not have the %= operator definition."));
        return *this;
    }
    
    bool mapParamsToArgs(const ParameterList &params, const std::vector<ArgInfo> signature, std::map<std::string, Element>* args) {
        args->clear();
        size_t numArgs = signature.size();
        std::vector<bool> assigned(numArgs, false);
        
        // find key-matched arguments defined in the function signature.
        for (auto namedParam : params.named) {
            const std::string key = namedParam.first;
            const Element value = namedParam.second;
            size_t idx = std::distance(std::begin(signature),
                                       std::find_if(std::begin(signature), std::end(signature),
                                                    [&key, &value](const ArgInfo &arg) {
                                                        return arg.name == key && (value.isConvertibleTo(arg.expectedType) || arg.expectedType == Type::Any);
                                                    }));
            // If params has a key which does not defined in the signatures, the mapping fails.
            if (idx == numArgs) {
                args->clear();
                return false;
            }
            const ArgInfo &argInfo = signature[idx];
            (*args)[argInfo.name] = argInfo.expectedType == Type::Any ? value : value.as(argInfo.expectedType);
            assigned[idx] = true;
        }
        // find arguments which have not assigned a value.
        for (auto it = std::begin(params.unnamed); it != std::end(params.unnamed); ++it) {
            const Element &value = *it;
            size_t idx = std::distance(std::begin(signature),
                                       std::find_if(std::begin(signature), std::end(signature),
                                                    [&signature, &value, &assigned](const ArgInfo &arg) {
                                                        size_t curIdx = std::distance(&signature[0], &arg);
                                                        return assigned[curIdx] == false && (value.isConvertibleTo(arg.expectedType) || arg.expectedType == Type::Any);
                                                    }));
            // If params has an extra parameters, the mapping fails.
            if (idx == numArgs) {
                args->clear();
                return false;
            }
            const ArgInfo &argInfo = signature[idx];
            (*args)[argInfo.name] = argInfo.expectedType == Type::Any ? value : value.as(argInfo.expectedType);
            assigned[idx] = true;
        }
        // If there are arguments they have not assigned yet and does not have default values, the mapping fails.
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
    }
    
    Element Function::operator()(const ParameterList &params, ExecuteContext &context, ErrorMessage* errMsg) const {
        std::map<std::string, Element> args;
        for (int i = 0; i < m_signatures.size(); ++i) {
            std::vector<ArgInfo> signature = m_signatures[i];
            const StatementRef &stmt = m_stmts[i];
            const Procedure &nativeProc = m_nativeProcs[i];
            
            if (mapParamsToArgs(params, signature, &args)) {
                context.stackVariables.current().saveFrom(m_depth);
                context.stackVariables.push();
                LocalVariables &current = context.stackVariables.current();
                for (auto &keyValue : args)
                    current[keyValue.first] = keyValue.second;
                
                Element retValue;
                context.returnFlag = false;
                if (stmt) {
                    stmt->perform(context, errMsg);
                    retValue = context.returnValue;
                }
                else {
                    retValue = nativeProc(args, context, errMsg);
                }
                context.returnValue = Element();
                context.returnFlag = false;
                
                context.stackVariables.pop();
                context.stackVariables.current().restore();
                return retValue;
            }
        }
        *errMsg = ErrorMessage("Parameters are invalid.");
        return Element();
    }
    
    bool TypeInfo::isConvertibleTo(Type toType) const {
        return convertFunctions.count(toType) > 0;
    }
    
    bool TypeInfo::initialized = false;
    TypeInfo TypeInfo::infos[(uint32_t)Type::NumTypes];
    
    void TypeInfo::init() {
        if (!initialized) {
            for (int i = 0; i < (int)Type::NumTypes; ++i) {
                TypeInfo &info = infos[i];
                info.postIncOperator = nullptr;
                info.postDecOperator = nullptr;
                info.preIncOperator = nullptr;
                info.preDecOperator = nullptr;
                info.affOperator = nullptr;
                info.negOperator = nullptr;
                info.logicNotOperator = nullptr;
                info.mulOperator = nullptr;
                info.divOperator = nullptr;
                info.remOperator = nullptr;
                info.addOperator = nullptr;
                info.subOperator = nullptr;
                info.lessOperator = nullptr;
                info.greaterOperator = nullptr;
                info.lessEqOperator = [](const Element &v0, const Element &v1) {
                    const TypeInfo &info = TypeInfo::infos[(uint32_t)v0.type];
                    return Element(info.eqOperator(v0, v1).raw<TypeMap::Bool>() || info.lessOperator(v0, v1).raw<TypeMap::Bool>());
                };
                info.greaterEqOperator = [](const Element &v0, const Element &v1) {
                    const TypeInfo &info = TypeInfo::infos[(uint32_t)v0.type];
                    return Element(info.eqOperator(v0, v1).raw<TypeMap::Bool>() || info.greaterOperator(v0, v1).raw<TypeMap::Bool>());
                };
                info.eqOperator = nullptr;
                info.neqOperator = [](const Element &v0, const Element &v1) {
                    const TypeInfo &info = TypeInfo::infos[(uint32_t)v0.type];
                    return Element(!info.eqOperator(v0, v1).raw<TypeMap::Bool>());
                };
                info.logicAndOperator = [](const Element &v0, const Element &v1) {
                    if (!v0.isConvertibleTo<TypeMap::Bool>())
                        return Element(ErrorMessage("Left operand cannot be converted to a bool value."));
                    if (!v1.isConvertibleTo<TypeMap::Bool>())
                        return Element(ErrorMessage("Right operand cannot be converted to a bool value."));
                    return Element(v0.asRaw<TypeMap::Bool>() && v1.asRaw<TypeMap::Bool>());
                };
                info.logicOrOperator = [](const Element &v0, const Element &v1) {
                    if (!v0.isConvertibleTo<TypeMap::Bool>())
                        return Element(ErrorMessage("Left operand cannot be converted to a bool value."));
                    if (!v1.isConvertibleTo<TypeMap::Bool>())
                        return Element(ErrorMessage("Right operand cannot be converted to a bool value."));
                    return Element(v0.asRaw<TypeMap::Bool>() || v1.asRaw<TypeMap::Bool>());
                };
                info.substOperator = [](Element &l, const Element &r) { l = r; };
                info.mulSubstOperator = [](Element &l, const Element &r) { l = l * r; };
                info.divSubstOperator = [](Element &l, const Element &r) { l = l / r; };
                info.remSubstOperator = [](Element &l, const Element &r) { l = l % r; };
                info.addSubstOperator = [](Element &l, const Element &r) { l = l + r; };
                info.subSubstOperator = [](Element &l, const Element &r) { l = l - r; };
                info.convertFunctions[(Type)i] = [](const Element &v) { return v; };
            }
            
            // Bool
            {
                TypeInfo &info = infos[(uint32_t)Type::Bool];
                
                info.affOperator = [](const Element &v) { return Element(v.raw<TypeMap::Bool>()); };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::Bool>()); };
                info.logicNotOperator = [](const Element &v) { return Element(!v.raw<TypeMap::Bool>()); };
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal * v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal * v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                info.divOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal / v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal / v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("/ operator does not support the right operand type."));
                };
                info.remOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal % v1.asRaw<TypeMap::Integer>());
                    return Element(ErrorMessage("% operator does not support the right operand type."));
                };
                info.addOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal + v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal + v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal - v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal - v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.lessOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal < v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal < v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("< operator does not support the right operand type."));
                };
                info.greaterOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal > v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal > v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("> operator does not support the right operand type."));
                };
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Bool>();
                    if (v1.isConvertibleTo<TypeMap::Bool>())
                        return Element(lVal == v1.asRaw<TypeMap::Bool>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Integer] = [](const Element &elemFrom) { return Element(TypeMap::Integer(), elemFrom.raw<TypeMap::Bool>()); };
                info.convertFunctions[Type::RealNumber] = [](const Element &elemFrom) { return Element(TypeMap::RealNumber(), elemFrom.raw<TypeMap::Bool>()); };
            }
            
            // Integer
            {
                TypeInfo &info = infos[(uint32_t)Type::Integer];
                
                info.postIncOperator = [](Element &v) {
                    Element old = v;
                    v = Element(v.raw<TypeMap::Integer>() + 1);
                    return old;
                };
                info.postDecOperator = [](Element &v) {
                    Element old = v;
                    v = Element(v.raw<TypeMap::Integer>() - 1);
                    return old;
                };
                info.preIncOperator = [](Element &v) { v = Element(v.raw<TypeMap::Integer>() + 1); };
                info.preDecOperator = [](Element &v) { v = Element(v.raw<TypeMap::Integer>() - 1); };
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::Integer>()); };
                info.logicNotOperator = [](const Element &v) { return Element(!v.raw<TypeMap::Integer>()); };
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal * v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal * v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                info.divOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal / v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal / v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("/ operator does not support the right operand type."));
                };
                info.remOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal % v1.asRaw<TypeMap::Integer>());
                    return Element(ErrorMessage("% operator does not support the right operand type."));
                };
                info.addOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal + v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal + v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal - v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal - v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.lessOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal < v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal < v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("< operator does not support the right operand type."));
                };
                info.greaterOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal > v1.asRaw<TypeMap::Integer>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal > v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("> operator does not support the right operand type."));
                };
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Integer>();
                    if (v1.isConvertibleTo<TypeMap::Integer>())
                        return Element(lVal == v1.asRaw<TypeMap::Bool>());
                    else if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal == v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Bool] = [](const Element &elemFrom) { return Element(TypeMap::Bool(), elemFrom.raw<TypeMap::Integer>()); };
                info.convertFunctions[Type::RealNumber] = [](const Element &elemFrom) { return Element(TypeMap::RealNumber(), elemFrom.raw<TypeMap::Integer>()); };
            }
            
            // Real Number
            {
                TypeInfo &info = infos[(uint32_t)Type::RealNumber];
                
                info.postIncOperator = [](Element &v) {
                    Element old = v;
                    v = Element(v.raw<TypeMap::RealNumber>() + 1);
                    return old;
                };
                info.postDecOperator = [](Element &v) {
                    Element old = v;
                    v = Element(v.raw<TypeMap::RealNumber>() - 1);
                    return old;
                };
                info.preIncOperator = [](Element &v) { v = Element(v.raw<TypeMap::RealNumber>() + 1); };
                info.preDecOperator = [](Element &v) { v = Element(v.raw<TypeMap::RealNumber>() - 1); };
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::RealNumber>()); };
                info.logicNotOperator = [](const Element &v) { return Element(!v.raw<TypeMap::RealNumber>()); };
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal * v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                info.divOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal / v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("/ operator does not support the right operand type."));
                };
                info.addOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal + v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal - v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.lessOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal < v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("< operator does not support the right operand type."));
                };
                info.greaterOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal > v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("> operator does not support the right operand type."));
                };
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::RealNumber>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(lVal == v1.asRaw<TypeMap::RealNumber>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
            }
            
            // String
            {
                TypeInfo &info = infos[(uint32_t)Type::String];
                
                info.addOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::String>();
                    if (v1.isConvertibleTo<TypeMap::String>())
                        return Element(lVal + v1.asRaw<TypeMap::String>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::String>();
                    if (v1.isConvertibleTo<TypeMap::String>())
                        return Element(lVal == v1.asRaw<TypeMap::String>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
            }
            
            // Point
            {
                TypeInfo &info = infos[(uint32_t)Type::Point];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::Point>()); };
                
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Point>();
                    if (v1.isConvertibleTo<TypeMap::Point>())
                        return Element(lVal == v1.asRaw<TypeMap::Point>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
            }
            
            // Vector
            {
                TypeInfo &info = infos[(uint32_t)Type::Vector];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::Vector>()); };
                
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Vector>();
                    if (v1.isConvertibleTo<TypeMap::Vector>())
                        return Element(lVal == v1.asRaw<TypeMap::Vector>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
            }
            
            // Normal
            {
                TypeInfo &info = infos[(uint32_t)Type::Normal];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.raw<TypeMap::Normal>()); };
                
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Normal>();
                    if (v1.isConvertibleTo<TypeMap::Normal>())
                        return Element(lVal == v1.asRaw<TypeMap::Normal>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Vector] = [](const Element &elemFrom) { return Element(SLR::Vector3D(elemFrom.raw<TypeMap::Normal>())); };
            }
            
            // Matrix
            {
                TypeInfo &info = infos[(uint32_t)Type::Matrix];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(TypeMap::Matrix(), -v.raw<TypeMap::Matrix>()); };
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Matrix>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(TypeMap::Matrix(), lVal * v1.asRaw<TypeMap::RealNumber>());
                    else if (v1.isConvertibleTo<TypeMap::Vertex>()) {
                        auto rVal = v1.asRaw<TypeMap::Vertex>();
                        SLR::Vertex vtx;
                        vtx.position = lVal * rVal.position;
                        vtx.normal = lVal * rVal.normal;
                        vtx.tangent = lVal * rVal.tangent;
                        vtx.texCoord = rVal.texCoord;
                        return Element(TypeMap::Vertex(), vtx);
                    }
                    else if (v1.isConvertibleTo<TypeMap::Matrix>())
                        return Element(TypeMap::Matrix(), lVal * v1.asRaw<TypeMap::Matrix>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                info.eqOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.raw<TypeMap::Matrix>();
                    if (v1.isConvertibleTo<TypeMap::Matrix>())
                        return Element(lVal == v1.asRaw<TypeMap::Matrix>());
                    return Element(ErrorMessage("== operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Transform] = [](const Element &elemFrom) {
                    return Element(TypeMap::Transform(), createShared<SLR::StaticTransform>(elemFrom.raw<TypeMap::Matrix>()));
                };
            }
            
            // Spectrum
            {
                TypeInfo &info = infos[(uint32_t)Type::Spectrum];
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    auto lVal = v0.rawRef<TypeMap::Spectrum>();
                    if (v1.isConvertibleTo<TypeMap::RealNumber>())
                        return Element(TypeMap::Spectrum(), InputSpectrumRef(lVal->createScaled(v1.asRaw<TypeMap::RealNumber>())));
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
            }
            initialized = true;
        }
    }
    
    BlockStatement::BlockStatement(const StatementsRef &statements) {
        for (int i = 0; i < statements->size(); ++i) {
            StatementRef &statement = statements->at(i);
            m_statements.push_back(statement);
        }
    }
    
    bool BlockStatement::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        context.stackVariables.current().pushDepth();
        for (int i = 0; i < m_statements.size(); ++i) {
            if (!m_statements[i]->perform(context, errMsg))
                return false;
            if (context.returnFlag)
                break;
        }
        context.stackVariables.current().popDepth();
        
        return true;
    }
    
    bool IfElseStatement::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        if (!m_condExpr->perform(context, errMsg))
            return false;
        if (!m_condExpr->result().isConvertibleTo<TypeMap::Bool>()) {
            *errMsg = ErrorMessage("Must provide a boolean value.");
            return false;
        }
        bool condition = m_condExpr->result().asRaw<TypeMap::Bool>();
        if (condition) {
            if (!m_trueBlockStmt->perform(context, errMsg))
                return false;
        }
        else if (m_falseBlockStmt) {
            if (!m_falseBlockStmt->perform(context, errMsg))
                return false;
        }
        
        return true;
    }
    
    bool ForStatement::perform(ExecuteContext &context, ErrorMessage* errMsg) const {
        if (!m_preExpr->perform(context, errMsg))
            return false;
        
        if (!m_condExpr->perform(context, errMsg))
            return false;
        if (!m_condExpr->result().isConvertibleTo<TypeMap::Bool>()) {
            *errMsg = ErrorMessage("Must provide a boolean value.");
            return false;
        }
        bool condition = m_condExpr->result().asRaw<TypeMap::Bool>();
        while (condition) {
            if (!m_blockStmt->perform(context, errMsg))
                return false;
            
            if (!m_postExpr->perform(context, errMsg))
                return false;
            
            if (!m_condExpr->perform(context, errMsg))
                return false;
            condition = m_condExpr->result().asRaw<TypeMap::Bool>();
        }
        
        return true;
    }
    
    bool FunctionDefinitionStatement::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        std::vector<ArgInfo> args;
        for (int i = 0; i < m_argDefs->size(); ++i) {
            args.emplace_back();
            
            ArgumentDefinitionRef argDef = m_argDefs->at(i);
            argDef->perform(context, errMsg);
            argDef->getArgInfo(&args.back());
        }
        
        LocalVariables &current = context.stackVariables.current();
        current[m_funcName] = Element(TypeMap::Function(), createShared<Function>(context.stackVariables.current().depth(), args, m_blockStmt));
        return true;
    }
    
    bool ReturnStatement::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        if (!m_expr) {
            context.returnValue = Element();
            context.returnFlag = true;
            return true;
        }
        if (!m_expr->perform(context, errMsg))
            return false;
        context.returnValue = m_expr->result();
        context.returnFlag = true;
        return true;
    }
    
    bool BinaryExpression::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_left->perform(context, errMsg))
            return false;
        if (!m_right->perform(context, errMsg))
            return false;
        
        if (m_op == "+")
            m_result = m_left->result() + m_right->result();
        else if (m_op == "-")
            m_result = m_left->result() - m_right->result();
        else if (m_op == "<")
            m_result = m_left->result() < m_right->result();
        else if (m_op == ">")
            m_result = m_left->result() > m_right->result();
        else if (m_op == "<=")
            m_result = m_left->result() <= m_right->result();
        else if (m_op == ">=")
            m_result = m_left->result() >= m_right->result();
        else if (m_op == "==")
            m_result = m_left->result() == m_right->result();
        else if (m_op == "!=")
            m_result = m_left->result() != m_right->result();
        else if (m_op == "&&")
            m_result = m_left->result() && m_right->result();
        else if (m_op == "||")
            m_result = m_left->result() || m_right->result();
        if (m_result.type == Type::Error) {
            *errMsg = m_result.raw<TypeMap::Error>();
            return false;
        }
        
        return true;
    }

    bool SubstitutionExpression::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_right->perform(context, errMsg))
            return false;
        if (m_op != "=" && context.stackVariables.exists(m_varName) == false) {
            *errMsg = ErrorMessage("Undefined variable: %s", m_varName.c_str());
            return false;
        }
        LocalVariables &current = context.stackVariables.current();
        Element &var = current[m_varName];
        
        if (m_op == "=")
            var.substitute(m_right->result());
        else if (m_op == "+=")
            var += m_right->result();
        else if (m_op == "-=")
            var -= m_right->result();
        else if (m_op == "*=")
            var *= m_right->result();
        else if (m_op == "/=")
            var /= m_right->result();
        else if (m_op == "%=")
            var %= m_right->result();
        if (var.type == Type::Error) {
            *errMsg = var.raw<TypeMap::Error>();
            return false;
        }
        
        m_result = var;
        return true;
    }
    
    bool UnaryTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_term->perform(context, errMsg))
            return false;
        if (m_op == "+")
            m_result = +m_term->result();
        else if (m_op == "-")
            m_result = -m_term->result();
        else if (m_op == "!")
            m_result = !m_term->result();
        return true;
    }
    
    bool UnarySubstitutionTerm::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        if (context.stackVariables.exists(m_varName) == false) {
            *errMsg = ErrorMessage("Undefined variable: %s", m_varName.c_str());
            return false;
        }
        LocalVariables &current = context.stackVariables.current();
        Element &var = current.at(m_varName);
        if (m_op == "++*")
            m_result = ++var;
        else if (m_op == "--*")
            m_result = --var;
        else if (m_op == "*++")
            m_result = var++;
        else if (m_op == "*--")
            m_result = var--;
        return true;
    }
    
    bool BinaryTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_left->perform(context, errMsg))
            return false;
        if (!m_right->perform(context, errMsg))
            return false;

        if (m_op == "*")
            m_result = m_left->result() * m_right->result();
        else if (m_op == "/")
            m_result = m_left->result() / m_right->result();
        else if (m_op == "%")
            m_result = m_left->result() % m_right->result();
        if (m_result.type == Type::Error) {
            *errMsg = m_result.raw<TypeMap::Error>();
            return false;
        }
        return true;
    }
    
    bool FunctionCallSingleTerm::perform(ExecuteContext &context, ErrorMessage* errMsg) const {
        if (context.stackVariables.exists(m_funcID)) {
            const Element &funcElem = context.stackVariables.at(m_funcID);
            if (funcElem.type == Type::Function) {
                ParameterList params;
                for (int i = 0; i < m_args->size(); ++i) {
                    ParameterRef arg = m_args->at(i);
                    if (!arg->perform(context, errMsg))
                        return false;
                    Element key, value;
                    arg->getKeyAndValue(&key, &value);
                    if (key.type != Type::String) {
                        *errMsg = ErrorMessage("Key expression must results in string type.");
                        return false;
                    }
                    params.add(key.raw<TypeMap::String>(), value);
                }
                
                const Function &func = funcElem.raw<TypeMap::Function>();
                m_result = func(params, context, errMsg);
                if (errMsg->error)
                    return false;
                
                return true;
            }
        }
        *errMsg = ErrorMessage("Function %s is not defined.", m_funcID.c_str());
        return false;
    }
    
    bool EnclosedSingleTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_expr->perform(context, errMsg))
            return false;
        m_result = m_expr->result();
        return true;
    }
    
    bool TupleElementSingleTerm::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        if (!m_tuple->perform(context, errMsg))
            return false;
        if (!m_idxExpr->perform(context, errMsg))
            return false;
        const Element &tuple = m_tuple->result();
        const Element &idx = m_idxExpr->result();
        if (!tuple.isConvertibleTo<TypeMap::Tuple>()) {
            *errMsg = ErrorMessage("Element access operator [] cannot be used to non tuple value.");
            return false;
        }
        
        const ParameterList &paramList = tuple.raw<TypeMap::Tuple>();
        if (idx.isConvertibleTo<TypeMap::Integer>()) {
            m_result = paramList(idx.raw<TypeMap::Integer>());
            if (m_result.type == Type::Void) {
                *errMsg = ErrorMessage("Index value is out or range.");
                return false;
            }
            return true;
        }
        else if (idx.isConvertibleTo<TypeMap::String>()) {
            m_result = paramList(idx.raw<TypeMap::String>());
            if (m_result.type == Type::Void) {
                *errMsg = ErrorMessage("Index value is invalid.");
                return false;
            }
            return true;
        }
        *errMsg = ErrorMessage("Index value must be integer or string compatible type.");
        return false;
    }
    
    bool TupleValue::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        ParameterListRef params = createShared<ParameterList>();
        for (int i = 0; i < m_elements->size(); ++i) {
            ParameterRef arg = m_elements->at(i);
            if (!arg->perform(context, errMsg))
                return false;
            Element key, value;
            arg->getKeyAndValue(&key, &value);
            if (key.type != Type::String) {
                *errMsg = ErrorMessage("Key expression must results in string type.");
                return false;
            }
            params->add(key.raw<TypeMap::String>(), value);
        }
        m_result = Element(TypeMap::Tuple(), params);
        return true;
    }
    
    bool VariableValue::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!context.stackVariables.exists(m_varName)) {
            *errMsg = ErrorMessage("Undefined variable is used.");
            return false;
        }
        m_result = context.stackVariables.at(m_varName);
        return true;
    }
    
    bool ArgumentDefinition::perform(SLRSceneGraph::ExecuteContext &context, SLRSceneGraph::ErrorMessage *errMsg) const {
        if (m_defaultValueExpr && !m_defaultValueExpr->perform(context, errMsg))
            return false;
        return true;
    }
    
    void ArgumentDefinition::getArgInfo(ArgInfo* info) const {
        info->name = m_name;
        info->defaultValue = m_defaultValueExpr ? m_defaultValueExpr->result() : Element();
        info->expectedType = Type::Any;
    }
    
    bool Parameter::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (m_keyExpr && !m_keyExpr->perform(context, errMsg))
            return false;
        if (!m_valueExpr->perform(context, errMsg))
            return false;
        return true;
    }
    
    void Parameter::getKeyAndValue(Element* key, Element* value) const {
        *key = m_keyExpr ? m_keyExpr->result() : Element(TypeMap::String(), "");
        *value = m_valueExpr->result();
    }
}
