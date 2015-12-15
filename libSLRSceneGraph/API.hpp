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

#include "Scene.h"
#include <libSLR/BasicTypes/Spectrum.h>

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
        NumTypes
    };
    std::ostream &operator<<(std::ostream &out, Type t);
    
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
    
    struct Element {
        Type type;
        std::shared_ptr<void> valueRef;
        Element(bool val) : type(Type::Bool), valueRef(createShared<bool>(val)) { };
        Element(int32_t val) : type(Type::Integer), valueRef(createShared<int32_t>(val)) { };
        Element(double val) : type(Type::RealNumber), valueRef(createShared<double>(val)) { };
        Element(const std::string &val) : type(Type::String), valueRef(createShared<std::string>(val)) { };
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
        typedef Element (*convertFunction)(const Element &);
        typedef Element (*unOpFunction)(const Element &);
        typedef Element (*binOpFunction)(const Element &, const Element &);
        unOpFunction negOperator;
        binOpFunction addOperator;
        binOpFunction subOperator;
        binOpFunction mulOperator;
        std::map<Type, convertFunction> convertFunctions;
        
        bool isConvertibleTo(Type toType) const;
        
        static bool initialized;
        static TypeInfo infos[(uint32_t)Type::NumTypes];
        static void init();
    };
    
    struct ErrorMessage {
        bool error;
        std::string message;
        ErrorMessage() : error(false) { };
        ErrorMessage(const std::string &msg) : error(true), message(msg) { };
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
        static const std::map<std::string, std::function<Element(const Element &, const Element &)>> s_functions;
        std::string m_op;
        TermRef m_left;
        TermRef m_right;
    public:
        BinaryExpression(const std::string &op, const TermRef &left, const TermRef &right) : m_op(op), m_left(left), m_right(right) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class SubstitutionExpression : public Expression {
        std::string m_varName;
        ExpressionRef m_right;
    public:
        SubstitutionExpression(const std::string &varName, const ExpressionRef &right) : m_varName(varName), m_right(right) { }
        
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
        static const std::map<std::string, std::function<Element(const Element &)>> s_functions;
        std::string m_op;
        TermRef m_term;
    public:
        UnaryTerm(const std::string &op, const TermRef &term) : m_op(op), m_term(term) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class BinaryTerm : public Term {
        static const std::map<std::string, std::function<Element(const Element &, const Element &)>> s_functions;
        std::string m_op;
        TermRef m_left;
        TermRef m_right;
    public:
        BinaryTerm(const std::string &op, const TermRef &left, const TermRef &right) : m_op(op), m_left(left), m_right(right) { }
        
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
    
    
    
    bool readScene(const std::string &filePath, Scene* scene, RenderingContext* context);
    
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
    
    Element CreatePerspectiveCamera(const ParameterList &params, ErrorMessage* err);
    
    Element SetRenderer(const ParameterList &params, RenderingContext* context, ErrorMessage* err);
    Element SetRenderSettings(const ParameterList &params, RenderingContext* context, ErrorMessage* err);
    
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
