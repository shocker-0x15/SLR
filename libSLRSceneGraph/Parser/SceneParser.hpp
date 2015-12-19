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
        Vertex,
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
        Any,
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
        TypeMapDef(Vertex, SLR::Vertex);
        TypeMapDef(Transform, SLR::Transform);
        TypeMapDef(Spectrum, SLR::InputSpectrum);
        TypeMapDef(SpectrumTexture, SLRSceneGraph::SpectrumTexture);
        TypeMapDef(NormalTexture, SLRSceneGraph::Normal3DTexture);
        TypeMapDef(FloatTexture, SLRSceneGraph::FloatTexture);
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
        AddItem,
        Translate,
        RotateX,
        RotateY,
        RotateZ,
        Scale,
        CreateVertex,
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
        
        template <typename Map>
        const typename Map::InternalType& raw() const { return *(typename Map::InternalType*)valueRef.get(); }
        template <typename Map>
        std::shared_ptr<typename Map::InternalType> rawRef() const { return std::static_pointer_cast<typename Map::InternalType>(valueRef); }
        bool isConvertibleTo(Type toType) const;
        template <typename Map>
        bool isConvertibleTo() const;
        Element as(Type toType) const;
        template <typename Map>
        Element as() const;
        template <typename Map>
        typename Map::InternalType asRaw() const;
        
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
    
    struct ParameterList {
        std::map<std::string, Element> named;
        std::vector<Element> unnamed;
        
        bool add(const std::string &key, const Element &value) {
            if (key.empty())
                unnamed.push_back(value);
            else {
                if (named.count(key))
                    return false;
                
                named[key] = value;
            }
            return true;
        };
        
        uint32_t numParams() const {
            return uint32_t(named.size() + unnamed.size());
        };
        
        Element operator()(const std::string &key, Type expectedType = Type::Any) const {
            if (named.count(key) == 0)
                return Element();
            
            const Element &value = named.at(key);
            if (expectedType == Type::Any)
                return value;
            else if (value.isConvertibleTo(expectedType))
                return value.as(expectedType);
            else
                return Element();
        };
        
        Element operator()(uint32_t idx, Type expectedType = Type::Any) const {
            if (idx >= unnamed.size())
                return Element();
            
            const Element &value = unnamed.at(idx);
            if (expectedType == Type::Any)
                return value;
            else if (value.isConvertibleTo(expectedType))
                return value.as(expectedType);
            else
                return Element();
        };
    };
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
    bool Element::isConvertibleTo() const {
        return TypeInfo::infos[(uint32_t)type].isConvertibleTo((Type)Map::Value);
    };
    
    template <typename Map>
    Element Element::as() const {
        SLRAssert(isConvertibleTo<Map>(), "Specified type is invalid.");
        TypeInfo::convertFunction func = TypeInfo::infos[(uint32_t)type].convertFunctions.at((Type)Map::Value);
        return func(*this);
    }
    
    template <typename Map>
    typename Map::InternalType Element::asRaw() const {
        SLRAssert(isConvertibleTo<Map>(), "Specified type is invalid.");
        TypeInfo::convertFunction func = TypeInfo::infos[(uint32_t)type].convertFunctions.at((Type)Map::Value);
        return func(*this).raw<Map>();
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
        
        void printVarList() const {
            auto depthwiseVars = m_depthwiseVars;
            while (!depthwiseVars.empty()) {
                std::set<std::string> vars = depthwiseVars.top();
                size_t depth = depthwiseVars.size();
                printf("%2lu:", depth);
                for (auto it = std::begin(vars); it != std::end(vars); ++it)
                    printf(" %s", it->c_str());
                printf("\n");
                depthwiseVars.pop();
            }
            printf("\n");
            depthwiseVars = m_depthwiseVars;
            while (!depthwiseVars.empty()) {
                std::set<std::string> vars = depthwiseVars.top();
                size_t depth = depthwiseVars.size();
                printf("%2lu:\n", depth);
                for (auto it = std::begin(vars); it != std::end(vars); ++it) {
                    const std::string &varName = *it;
                    printf(" %s: \n", varName.c_str());
                    std::cout << "  " << m_primaryMap.at(varName) << std::endl;
                }
                printf("\n");
                depthwiseVars.pop();
            }
        }
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
    
    class SingleTerm : public Term {
    };
    
    class Value : public SingleTerm {
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
    
    class UnaryTerm : public Term {
        std::string m_op;
        TermRef m_term;
    public:
        UnaryTerm(const std::string &op, const TermRef &term) : m_op(op), m_term(term) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class UnarySubstitutionTerm : public Term {
        std::string m_op;
        std::string m_varName;
    public:
        UnarySubstitutionTerm(const std::string &op, const std::string &varName) : m_op(op), m_varName(varName) { }
        
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
    
    class FunctionSingleTerm : public SingleTerm {
        API m_funcID;
        ParameterVecRef m_args;
    public:
        FunctionSingleTerm(API funcID, const ParameterVecRef &args) : m_funcID(funcID), m_args(args) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class EnclosedSingleTerm : public SingleTerm {
        ExpressionRef m_expr;
    public:
        EnclosedSingleTerm(const ExpressionRef &expr) : m_expr(expr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class TupleElementSingleTerm : public SingleTerm {
        SingleTermRef m_tuple;
        ExpressionRef m_idxExpr;
    public:
        TupleElementSingleTerm(const SingleTermRef &tuple, const ExpressionRef &idxExpr) : m_tuple(tuple), m_idxExpr(idxExpr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class ImmediateValue : public Value {
    public:
        ImmediateValue(const Element &value) { m_result = value; }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override { return true; }
    };
    
    class TupleValue : public Value {
        ParameterVecRef m_elements;
    public:
        TupleValue(const ParameterVecRef &elements) : m_elements(elements) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class VariableValue : public Value {
        std::string m_varName;
    public:
        VariableValue(const std::string &varName) : m_varName(varName) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class Parameter {
        ExpressionRef m_keyExpr;
        ExpressionRef m_valueExpr;
    public:
        Parameter(const ExpressionRef &keyExpr, const ExpressionRef &valueExpr) : m_keyExpr(keyExpr), m_valueExpr(valueExpr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const;
        void getKeyAndValue(Element* key, Element* value) const;
    };
}

#endif /* SceneParser_hpp */
