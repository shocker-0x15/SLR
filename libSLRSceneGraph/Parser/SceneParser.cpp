//
//  SceneParser.cpp
//  SLR
//
//  Created by 渡部 心 on 2015/12/16.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#include "SceneParser.hpp"
#include "API.hpp"

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
            case Type::Tuple:
                out << "Tuple";
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
            case Type::Tuple:
                out << elem.raw<TypeMap::Tuple>();
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
                    return Element(info.eqOperator(v0, v1).raw<TypeMap::Bool>() && info.lessOperator(v0, v1).raw<TypeMap::Bool>());
                };
                info.greaterEqOperator = [](const Element &v0, const Element &v1) {
                    const TypeInfo &info = TypeInfo::infos[(uint32_t)v0.type];
                    return Element(info.eqOperator(v0, v1).raw<TypeMap::Bool>() && info.greaterOperator(v0, v1).raw<TypeMap::Bool>());
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
            initialized = true;
        }
    }
    
    ForStatement::ForStatement(const ExpressionRef &preExpr, const ExpressionRef &condExpr, const ExpressionRef &postExpr, const StatementsRef &statementList) :
    m_preExpr(preExpr), m_condExpr(condExpr), m_postExpr(postExpr) {
        for (int i = 0; i < statementList->size(); ++i)
            m_block.push_back(statementList->at(i));
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
            context.stackVariables.pushDepth();
            for (int i = 0; i < m_block.size(); ++i) {
                if (!m_block[i]->perform(context, errMsg))
                    return false;
            }
            context.stackVariables.popDepth();
            if (!m_postExpr->perform(context, errMsg))
                return false;
            
            if (!m_condExpr->perform(context, errMsg))
                return false;
            condition = m_condExpr->result().asRaw<TypeMap::Bool>();
        }
        
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
        Element &var = context.stackVariables[m_varName];
        
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
        else if (m_op == "\%=")
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
        Element &var = context.stackVariables[m_varName];
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
        else if (m_op == "\%")
            m_result = m_left->result() % m_right->result();
        if (m_result.type == Type::Error) {
            *errMsg = m_result.raw<TypeMap::Error>();
            return false;
        }
        return true;
    }
    
    bool FunctionSingleTerm::perform(ExecuteContext &context, ErrorMessage* errMsg) const {
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
        
        switch (m_funcID) {
            case API::AddItem:
                m_result = AddItem(params, errMsg);
                break;
            case API::Translate:
                m_result = Translate(params, errMsg);
                break;
            case API::RotateX:
                m_result = RotateX(params, errMsg);
                break;
            case API::RotateY:
                m_result = RotateY(params, errMsg);
                break;
            case API::RotateZ:
                m_result = RotateZ(params, errMsg);
                break;
            case API::Scale:
                m_result = Scale(params, errMsg);
                break;
            case API::CreateVertex:
                m_result = CreateVertex(params, errMsg);
                break;
            case API::Spectrum:
                m_result = CreateSpectrum(params, errMsg);
                break;
            case API::SpectrumTexture:
                m_result = CreateSpectrumTexture(params, errMsg);
                break;
            case API::CreateMatte:
                m_result = CreateMatte(params, errMsg);
                break;
            case API::CreateDiffuseEmitter:
                m_result = CreateDiffuseEmitter(params, errMsg);
                break;
            case API::CreateEmitterSurfaceMaterial:
                m_result = CreateEmitterSurfaceMaterial(params, errMsg);
                break;
            case API::CreateMesh:
                m_result = CreateMesh(params, errMsg);
                break;
            case API::CreateNode:
                m_result = CreateNode(params, errMsg);
                break;
            case API::SetTransform:
                m_result = SetTransform(params, errMsg);
                break;
            case API::AddChild:
                m_result = AddChild(params, errMsg);
                break;
            case API::SetRenderer:
                m_result = SetRenderer(params, context.renderingContext, errMsg);
                break;
            case API::SetRenderSettings:
                m_result = SetRenderSettings(params, context.renderingContext, errMsg);
                break;
            case API::CreatePerspectiveCamera:
                m_result = CreatePerspectiveCamera(params, errMsg);
                break;
            case API::Load3DModel:
                m_result = Load3DModel(params, errMsg);
                break;
            case API::LoadImage:
                break;
            case API::SetEnvironment:
                break;
            default:
                break;
        }
        return !errMsg->error;
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
                *errMsg = ErrorMessage("index value is out or range.");
                return false;
            }
            return true;
        }
        else if (idx.isConvertibleTo<TypeMap::String>()) {
            m_result = paramList(idx.raw<TypeMap::String>());
            if (m_result.type == Type::Void) {
                *errMsg = ErrorMessage("index value is invalid.");
                return false;
            }
            return true;
        }
        *errMsg = ErrorMessage("index value must be integer or string compatible type.");
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
        m_result = context.stackVariables[m_varName];
        return true;
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
