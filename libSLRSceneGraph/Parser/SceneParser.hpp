//
//  SceneParser.hpp
//
//  Created by 渡部 心 on 2015/12/16.
//  Copyright © 2015年 渡部 心. All rights reserved.
//

#ifndef SceneParser_hpp
#define SceneParser_hpp

#include <libSLR/defines.h>
#include "references.h"

#include "Scene.h"

namespace SLRSceneGraph {
    enum class Type : uint32_t {
        Bool = 0,
        Integer,
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
        Camera,
        Node,
        Tuple,
        Void,
        Error,
        NumTypes
    };
    std::ostream &operator<<(std::ostream &out, Type t);
    
    namespace TypeMap {
#define TypeMapDef(TypeName, InternalTypeName) struct TypeName { enum : uint32_t { Value = (uint32_t)Type::TypeName }; typedef InternalTypeName InternalType; };
        TypeMapDef(Bool, bool);
        TypeMapDef(Integer, int32_t);
        TypeMapDef(RealNumber, double);
        TypeMapDef(String, std::string);
        TypeMapDef(Matrix, SLR::Matrix4x4);
        TypeMapDef(Transform, SLR::Transform);
        TypeMapDef(Spectrum, SLR::InputSpectrum);
        TypeMapDef(SpectrumTexture, SLRSceneGraph::SpectrumTexture);
        TypeMapDef(NormalTexture, Normal3DTexture);
        TypeMapDef(FloatTexture, FloatTexture);
        TypeMapDef(SurfaceMaterial, SLRSceneGraph::SurfaceMaterial);
        TypeMapDef(EmitterSurfaceProperty, SLRSceneGraph::EmitterSurfaceProperty);
        TypeMapDef(Mesh, SLRSceneGraph::TriangleMeshNode);
        TypeMapDef(Camera, SLRSceneGraph::CameraNode);
        TypeMapDef(Node, SLRSceneGraph::InternalNode);
        TypeMapDef(Tuple, SLRSceneGraph::ParameterList);
        TypeMapDef(Void, void);
        TypeMapDef(Error, SLRSceneGraph::ErrorMessage);
    }
    
    enum class API : uint32_t {
        Translate,
        RotateX,
        RotateY,
        RotateZ,
        Scale,
        Spectrum,
        SpectrumTexture,
        CreateMatte,
        CreateDiffuseEmitter,
        CreateEmitterSurfaceMaterial,
        CreateMesh,
        CreateNode,
        SetTransform,
        AddChild,
        CreatePerspectiveCamera,
        SetRenderer,
        SetRenderSettings,
        SetEnvironment,
        LoadImage,
        Load3DModel,
    };
    
    struct ErrorMessage {
        bool error;
        std::string message;
        ErrorMessage() : error(false) { }
        ErrorMessage(const std::string &msg) : error(true), message(msg) { }
        ErrorMessage(const char* fmt, ...) : error(true) {
            va_list args;
            va_start(args, fmt);
            char str[256];
            vsprintf(str, fmt, args);
            message = str;
            va_end(args);
        }
    };
    
    struct Element {
        Type type;
        std::shared_ptr<void> valueRef;
        Element(bool val) : type(Type::Bool), valueRef(createShared<bool>(val)) { }
        Element(int32_t val) : type(Type::Integer), valueRef(createShared<int32_t>(val)) { }
        Element(double val) : type(Type::RealNumber), valueRef(createShared<double>(val)) { }
        Element(const std::string &val) : type(Type::String), valueRef(createShared<std::string>(val)) { }
        Element(const ErrorMessage &val) : type(Type::Error), valueRef(createShared<ErrorMessage>(val)) { }
        
        Element() : type(Type::Void), valueRef(nullptr) { }
        template <typename Map>
        Element(Map dummy, const typename Map::InternalType &val) : type((Type)Map::Value), valueRef(createShared<typename Map::InternalType>(val)) { }
        template <typename Map>
        Element(Map dummy, const std::shared_ptr<typename Map::InternalType> &valRef = nullptr) : type((Type)Map::Value), valueRef(valRef) { }
//        Element(Type t = Type::Void, const std::shared_ptr<void> &vr = nullptr) : type(t), valueRef(vr) { };
        
        template <typename T>
        const T& as() const { return *(T*)valueRef.get(); }
        template <typename T>
        std::shared_ptr<T> asRef() const { return std::static_pointer_cast<T>(valueRef); }
        bool isConvertibleTo(Type toType) const;
        Element convertTo(Type toType) const;
        template <typename Map>
        typename Map::InternalType convertTo() const;
        
        Element operator++(int);
        Element operator--(int);
        Element &operator++();
        Element &operator--();
        Element operator+() const;
        Element operator-() const;
        Element operator!() const;
        
        Element operator*(const Element &rElem) const;
        Element operator/(const Element &rElem) const;
        Element operator%(const Element &rElem) const;
        Element operator+(const Element &rElem) const;
        Element operator-(const Element &rElem) const;
        
        Element operator<(const Element &rElem) const;
        Element operator>(const Element &rElem) const;
        Element operator<=(const Element &rElem) const;
        Element operator>=(const Element &rElem) const;
        Element operator==(const Element &rElem) const;
        Element operator!=(const Element &rElem) const;
        
        Element operator&&(const Element &rElem) const;
        Element operator||(const Element &rElem) const;
        
        Element &substitute(const Element &rElem);
        Element &operator+=(const Element &rElem);
        Element &operator-=(const Element &rElem);
        Element &operator*=(const Element &rElem);
        Element &operator/=(const Element &rElem);
        Element &operator%=(const Element &rElem);
    };
    std::ostream &operator<<(std::ostream &out, const Element &elem);
    
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
    std::ostream &operator<<(std::ostream &out, const ParameterList &params);
    
    struct TypeInfo {
        Element &operator%=(const Element &rElem) const;
        
        typedef Element (*convertFunction)(const Element &);
        typedef Element (*unOpFunction)(const Element &);
        typedef Element (*postSubstUnOpFunction)(Element &);
        typedef void (*preSubstUnOpFunction)(Element &);
        typedef Element (*binOpFunction)(const Element &, const Element &);
        typedef void (*substBinOpFunction)(Element &, const Element &);
        postSubstUnOpFunction postIncOperator;
        postSubstUnOpFunction postDecOperator;
        preSubstUnOpFunction preIncOperator;
        preSubstUnOpFunction preDecOperator;
        unOpFunction affOperator;
        unOpFunction negOperator;
        unOpFunction logicNotOperator;
        binOpFunction mulOperator;
        binOpFunction divOperator;
        binOpFunction remOperator;
        binOpFunction addOperator;
        binOpFunction subOperator;
        binOpFunction lessOperator;
        binOpFunction greaterOperator;
        binOpFunction lessEqOperator;
        binOpFunction greaterEqOperator;
        binOpFunction eqOperator;
        binOpFunction neqOperator;
        binOpFunction logicAndOperator;
        binOpFunction logicOrOperator;
        substBinOpFunction substOperator;
        substBinOpFunction mulSubstOperator;
        substBinOpFunction divSubstOperator;
        substBinOpFunction remSubstOperator;
        substBinOpFunction addSubstOperator;
        substBinOpFunction subSubstOperator;
        std::map<Type, convertFunction> convertFunctions;
        
        bool isConvertibleTo(Type toType) const;
        
        static bool initialized;
        static TypeInfo infos[(uint32_t)Type::NumTypes];
        static void init();
    };
    
    template <typename Map>
    typename Map::InternalType Element::convertTo() const {
        SLRAssert(isConvertibleTo((Type)Map::Value), "Specified type is invalid.");
        TypeInfo::convertFunction func = TypeInfo::infos[Map::Value].convertFunctions.at((Type)Map::Value);
        return func(*this).as<typename Map::InternalType>();
    };
    
    class LocalVariables {
        std::map<std::string, Element> m_primaryMap;
        std::stack<std::set<std::string>> m_depthwiseVars;
    public:
        LocalVariables() { m_depthwiseVars.emplace(); }
        void pushDepth() { m_depthwiseVars.emplace(); }
        void popDepth() {
            for (const std::string &name : m_depthwiseVars.top())
                m_primaryMap.erase(name);
            m_depthwiseVars.pop();
        }
        bool exists(const std::string &key) {
            return m_primaryMap.count(key) > 0;
        }
        Element &operator[](const std::string &key) {
            if (m_primaryMap.count(key) == 0)
                m_depthwiseVars.top().insert(key);
            return m_primaryMap[key];
        }
        const Element &at(const std::string &key) {
            return m_primaryMap.at(key);
        }
        const Element &operator[](const std::string &key) const { return m_primaryMap.at(key); }
    };
    
    struct ExecuteContext {
        LocalVariables stackVariables;
        Scene* scene;
        RenderingContext* renderingContext;
    };
    
    
    class Statement {
    public:
        virtual bool perform(ExecuteContext &context, ErrorMessage* errMsg) const = 0;
    };
    
    class Expression : public Statement {
    protected:
        mutable Element m_result;
    public:
        const Element &result() const { return m_result; };
    };
    
    class Term : public Expression {
    };
    
    class Value : public Term {
    };
    
    class Argument {
        ExpressionRef m_keyExpr;
        ExpressionRef m_valueExpr;
        mutable Parameter m_result;
    public:
        Argument(const ExpressionRef &keyExpr, const ExpressionRef &valueExpr) : m_keyExpr(keyExpr), m_valueExpr(valueExpr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const;
        const Parameter &result() const { return m_result; };
    };
    
    class ForStatement : public Statement {
        ExpressionRef m_preExpr;
        ExpressionRef m_condExpr;
        ExpressionRef m_postExpr;
        std::vector<StatementRef> m_block;
    public:
        ForStatement(const ExpressionRef &preExpr, const ExpressionRef &condExpr, const ExpressionRef &postExpr, const StatementsRef &statementList);
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class BinaryExpression : public Expression {
        TermRef m_left;
        std::string m_op;
        TermRef m_right;
    public:
        BinaryExpression(const TermRef &left, const std::string &op, const TermRef &right) : m_left(left), m_op(op), m_right(right) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class SubstitutionExpression : public Expression {
        std::string m_varName;
        std::string m_op;
        ExpressionRef m_right;
    public:
        SubstitutionExpression(const std::string &varName, const std::string &op, const ExpressionRef &right) : m_varName(varName), m_op(op), m_right(right) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class FunctionTerm : public Term {
        API m_funcID;
        ArgumentsRef m_args;
    public:
        FunctionTerm(API funcID, const ArgumentsRef &args) : m_funcID(funcID), m_args(args) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class UnaryTerm : public Term {
        std::string m_op;
        TermRef m_term;
    public:
        UnaryTerm(const std::string &op, const TermRef &term) : m_op(op), m_term(term) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class BinaryTerm : public Term {
        TermRef m_left;
        std::string m_op;
        TermRef m_right;
    public:
        BinaryTerm(const TermRef &left, const std::string &op, const TermRef &right) : m_left(left), m_op(op), m_right(right) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class EnclosedTerm : public Term {
        ExpressionRef m_expr;
    public:
        EnclosedTerm(const ExpressionRef &expr) : m_expr(expr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class ImmediateValue : public Value {
    public:
        ImmediateValue(const Element &value) { m_result = value; }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override { return true; }
    };
    
    class TupleValue : public Value {
        ArgumentsRef m_elements;
    public:
        TupleValue(const ArgumentsRef &elements) : m_elements(elements) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class VariableValue : public Value {
        std::string m_varName;
    public:
        VariableValue(const std::string &varName) : m_varName(varName) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
}

#endif /* SceneParser_hpp */
