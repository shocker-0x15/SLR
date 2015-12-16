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
            case Type::Void:
                out << "Void";
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
            case Type::Integer:
                out << elem.as<int32_t>();
                break;
            case Type::RealNumber:
                out << elem.as<double>();
                break;
            case Type::String:
                out << "\"" << elem.as<std::string>() << "\"";
                break;
            case Type::Matrix:
                out << "Matrix";
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
                out << elem.as<ParameterList>();
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
    }
    
    Element Element::convertTo(Type toType) const {
        SLRAssert(isConvertibleTo(toType), "Specified type is invalid.");
        TypeInfo::convertFunction func = TypeInfo::infos[(uint32_t)type].convertFunctions.at(toType);
        return func(*this);
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
    
    Element Element::operator*(const Element &rElem) const {
        const TypeInfo &leftType = TypeInfo::infos[(uint32_t)type];
        if (!leftType.mulOperator)
            return Element(ErrorMessage("Left type does not have the * operator definition."));
        return leftType.mulOperator(*this, rElem);
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
    
    bool TypeInfo::isConvertibleTo(Type toType) const {
        return convertFunctions.count(toType) > 0;
    }
    
    bool TypeInfo::initialized = false;
    TypeInfo TypeInfo::infos[(uint32_t)Type::NumTypes];
    
    void TypeInfo::init() {
        if (!initialized) {
            for (int i = 0; i < (int)Type::NumTypes; ++i) {
                TypeInfo &info = infos[i];
                info.affOperator = nullptr;
                info.negOperator = nullptr;
                info.addOperator = nullptr;
                info.subOperator = nullptr;
                info.mulOperator = nullptr;
                info.convertFunctions[(Type)i] = [](const Element &v) { return v; };
            }
            
            {
                TypeInfo &info = infos[(uint32_t)Type::Bool];
                
                info.affOperator = [](const Element &v) { return Element(v.as<bool>()); };
                info.negOperator = [](const Element &v) { return Element(-(int32_t)v.as<bool>()); };
                
                info.addOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<bool>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal + v1.convertTo<TypeMap::Integer>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal + v1.convertTo<TypeMap::RealNumber>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<bool>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal - v1.convertTo(Type::Integer).as<int32_t>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal - v1.convertTo(Type::RealNumber).as<int32_t>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<bool>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal * v1.convertTo(Type::Integer).as<int32_t>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal * v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Integer] = [](const Element &elemFrom) { return Element((int32_t)elemFrom.as<bool>()); };
                info.convertFunctions[Type::RealNumber] = [](const Element &elemFrom) { return Element((double)elemFrom.as<bool>()); };
            }
            
            {
                TypeInfo &info = infos[(uint32_t)Type::Integer];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.as<int32_t>()); };
                
                info.addOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<int32_t>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal + v1.convertTo(Type::Integer).as<int32_t>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal + v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<int32_t>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal - v1.convertTo(Type::Integer).as<int32_t>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal - v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    int32_t lVal = v0.as<int32_t>();
                    if (v1.isConvertibleTo(Type::Integer))
                        return Element(lVal * v1.convertTo(Type::Integer).as<int32_t>());
                    else if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal * v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Bool] = [](const Element &elemFrom) { return Element((bool)elemFrom.as<int32_t>()); };
                info.convertFunctions[Type::RealNumber] = [](const Element &elemFrom) { return Element((double)elemFrom.as<int32_t>()); };
            }
            {
                TypeInfo &info = infos[(uint32_t)Type::RealNumber];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(-v.as<double>()); };
                
                info.addOperator = [](const Element &v0, const Element &v1) {
                    double lVal = v0.as<double>();
                    if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal + v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("+ operator does not support the right operand type."));
                };
                info.subOperator = [](const Element &v0, const Element &v1) {
                    double lVal = v0.as<double>();
                    if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal - v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("- operator does not support the right operand type."));
                };
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    double lVal = v0.as<double>();
                    if (v1.isConvertibleTo(Type::RealNumber))
                        return Element(lVal * v1.convertTo(Type::RealNumber).as<double>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
            }
            {
                TypeInfo &info = infos[(uint32_t)Type::Matrix];
                
                info.affOperator = [](const Element &v) { return v; };
                info.negOperator = [](const Element &v) { return Element(TypeMap::Matrix(), -v.as<SLR::Matrix4x4>()); };
                
                info.mulOperator = [](const Element &v0, const Element &v1) {
                    SLR::Matrix4x4 lVal = v0.as<SLR::Matrix4x4>();
                    if (v1.isConvertibleTo(Type::Matrix))
                        return Element(TypeMap::Matrix(), lVal * v1.convertTo<TypeMap::Matrix>());
                    return Element(ErrorMessage("* operator does not support the right operand type."));
                };
                
                info.convertFunctions[Type::Transform] = [](const Element &elemFrom) {
                    return Element(TypeMap::Transform(), createShared<SLR::StaticTransform>(elemFrom.as<SLR::Matrix4x4>()));
                };
            }
            initialized = true;
        }
    }
    
    bool Argument::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (m_keyExpr && !m_keyExpr->perform(context, errMsg))
            return false;
        if (!m_valueExpr->perform(context, errMsg))
            return false;
        m_result = Parameter(m_keyExpr ? m_keyExpr->result().as<std::string>() : "", m_valueExpr->result());
        return true;
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
        if (!m_condExpr->result().isConvertibleTo(Type::Bool)) {
            *errMsg = ErrorMessage("Must provide a boolean value.");
            return false;
        }
        bool condition = m_condExpr->result().as<bool>();
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
            condition = m_condExpr->result().as<bool>();
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
        if (m_result.type == Type::Error) {
            *errMsg = m_result.as<ErrorMessage>();
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
        
        if (m_op == "=")
            m_result = m_right->result();
        else if (m_op == "+=")
            m_result = context.stackVariables.at(m_varName) + m_right->result();
        else if (m_op == "-=")
            m_result = context.stackVariables.at(m_varName) - m_right->result();
        else if (m_op == "*=")
            m_result = context.stackVariables.at(m_varName) * m_right->result();
        if (m_result.type == Type::Error) {
            *errMsg = m_result.as<ErrorMessage>();
            return false;
        }
        
        context.stackVariables[m_varName] = m_result;
        return true;
    }
    
    bool FunctionTerm::perform(ExecuteContext &context, ErrorMessage* errMsg) const {
        ParameterList params;
        for (int i = 0; i < m_args->size(); ++i) {
            ArgumentRef arg = m_args->at(i);
            if (!arg->perform(context, errMsg))
                return false;
            params.add(arg->result());
        }
        
        switch (m_funcID) {
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
    
    bool UnaryTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_term->perform(context, errMsg))
            return false;
        if (m_op == "+")
            m_result = +m_term->result();
        else if (m_op == "-")
            m_result = -m_term->result();
        return true;
    }
    
    bool BinaryTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_left->perform(context, errMsg))
            return false;
        if (!m_right->perform(context, errMsg))
            return false;
        
        if (m_op == "*")
            m_result = m_left->result() * m_right->result();
        if (m_result.type == Type::Error) {
            *errMsg = m_result.as<ErrorMessage>();
            return false;
        }
        return true;
    }
    
    bool EnclosedTerm::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!m_expr->perform(context, errMsg))
            return false;
        m_result = m_expr->result();
        return true;
    }
    
    bool TupleValue::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        ParameterListRef params = createShared<ParameterList>();
        for (int i = 0; i < m_elements->size(); ++i) {
            ArgumentRef arg = m_elements->at(i);
            if (!arg->perform(context, errMsg))
                return false;
            params->add(arg->result());
        }
        m_result = Element(TypeMap::Tuple(), params);
        return true;
    }
    
    bool VariableValue::perform(ExecuteContext &context, ErrorMessage *errMsg) const {
        if (!context.stackVariables.exists(m_varName))
            return false;
        m_result = context.stackVariables[m_varName];
        return true;
    }    
}
