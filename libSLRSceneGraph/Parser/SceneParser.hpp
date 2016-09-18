//
//  SceneParser.hpp
//
//  Created by 渡部 心 on 2015/12/16.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SceneParser_hpp
#define SceneParser_hpp

#include <libSLR/defines.h>
#include "../references.h"

#include "../Scene.h"

namespace SLRSceneGraph {
    enum class Type : uint32_t {
        Bool = 0,
        Integer,
        RealNumber,
        String,
        Point,
        Vector,
        Normal,
        Matrix,
        Vertex,
        Transform,
        Spectrum,
        Image2D,
        Texture2DMapping,
        Texture3DMapping,
        SpectrumTexture,
        NormalTexture,
        FloatTexture,
        SurfaceMaterial,
        EmitterSurfaceProperty,
        Mesh,
        Camera,
        Node,
        ReferenceNode,
        Tuple,
        Function, 
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
        TypeMapDef(Point, SLR::Point3D);
        TypeMapDef(Vector, SLR::Vector3D);
        TypeMapDef(Normal, SLR::Normal3D);
        TypeMapDef(Matrix, SLR::Matrix4x4);
        TypeMapDef(Vertex, SLR::Vertex);
        TypeMapDef(Transform, SLR::Transform);
        TypeMapDef(Spectrum, SLR::InputSpectrum);
        TypeMapDef(Image2D, SLR::TiledImage2D);
        TypeMapDef(Texture2DMapping, SLRSceneGraph::Texture2DMapping);
        TypeMapDef(Texture3DMapping, SLRSceneGraph::Texture3DMapping);
        TypeMapDef(SpectrumTexture, SLRSceneGraph::SpectrumTexture);
        TypeMapDef(NormalTexture, SLRSceneGraph::Normal3DTexture);
        TypeMapDef(FloatTexture, SLRSceneGraph::FloatTexture);
        TypeMapDef(SurfaceMaterial, SLRSceneGraph::SurfaceMaterial);
        TypeMapDef(EmitterSurfaceProperty, SLRSceneGraph::EmitterSurfaceProperty);
        TypeMapDef(Mesh, SLRSceneGraph::TriangleMeshNode);
        TypeMapDef(Camera, SLRSceneGraph::CameraNode);
        TypeMapDef(Node, SLRSceneGraph::InternalNode);
        TypeMapDef(ReferenceNode, SLRSceneGraph::ReferenceNode);
        TypeMapDef(Tuple, SLRSceneGraph::ParameterList);
        TypeMapDef(Function, SLRSceneGraph::Function);
        TypeMapDef(Void, void);
        TypeMapDef(Error, SLRSceneGraph::ErrorMessage);
#undef TypeMapDef
    }
    
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
        Element(const SLR::Point3D &val);
        Element(const SLR::Vector3D &val);
        Element(const SLR::Normal3D &val);
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
        typename Map::InternalType& raw() { return *(typename Map::InternalType*)valueRef.get(); }
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
        
        size_t numNamed() const { return named.size(); }
        size_t numUnnamed() const { return unnamed.size(); }
        size_t numParams() const { return named.size() + unnamed.size(); };
        
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
    
    struct ArgInfo {
        std::string name;
        Type expectedType;
        Element defaultValue;
    };
    
    bool mapParamsToArgs(const ParameterList &params, const std::vector<ArgInfo> signature, std::map<std::string, Element>* args);
    
    class Function {
        typedef std::function<Element(const std::map<std::string, Element> &, ExecuteContext &, ErrorMessage*)> Procedure;
        
		// workaround: MSVC(VS2015 Update1) fails to compile the 3rd constructor.
#ifdef SLR_Platform_Windows_MSVC
		const uint32_t m_depth;
		const std::vector<std::vector<ArgInfo>> m_signatures;
		std::vector<Procedure> m_nativeProcs;
		const std::vector<StatementRef> m_stmts;

		static Element NoOpProcedure(const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) { SLRAssert_NotImplemented(); return Element(); };
	public:
		Function(uint32_t depth, const std::vector<ArgInfo> &sig, const StatementRef &stmt = nullptr) :
			m_depth(depth), m_signatures{ sig }, m_nativeProcs{ NoOpProcedure }, m_stmts{ stmt } { }
		Function(uint32_t depth, const std::vector<ArgInfo> &sig, const Procedure &proc) :
			m_depth(depth), m_signatures{ sig }, m_nativeProcs{ proc }, m_stmts{ nullptr } { }
		Function(uint32_t depth, const std::vector<std::vector<ArgInfo>> &sigs, const std::vector<Procedure> &procs) :
			m_depth(depth), m_signatures{ sigs }, m_stmts{ m_signatures.size(), nullptr } {
			for (int i = 0; i < procs.size(); ++i)
				m_nativeProcs.push_back(procs[i]);
		}
#else
        const uint32_t m_depth;
        const std::vector<std::vector<ArgInfo>> m_signatures;
        const std::vector<Procedure> m_nativeProcs;
        const std::vector<StatementRef> m_stmts;
        
        static Element NoOpProcedure(const std::map<std::string, Element> &args, ExecuteContext &context, ErrorMessage* err) { SLRAssert_NotImplemented(); return Element(); };
    public:
        Function(uint32_t depth, const std::vector<ArgInfo> &sig, const StatementRef &stmt = nullptr) :
        m_depth(depth), m_signatures{sig}, m_nativeProcs{NoOpProcedure}, m_stmts{stmt} { }
        Function(uint32_t depth, const std::vector<ArgInfo> &sig, const Procedure &proc) :
        m_depth(depth), m_signatures{sig}, m_nativeProcs{proc}, m_stmts{nullptr} { }
        Function(uint32_t depth, const std::vector<std::vector<ArgInfo>> &sigs, const std::vector<Procedure> &procs) :
        m_depth(depth), m_signatures{sigs}, m_nativeProcs{procs}, m_stmts{m_signatures.size(), nullptr} { }
#endif
        
        Element operator()(const ParameterList &params, ExecuteContext &context, ErrorMessage* errMsg) const;
        
        template <typename T>
        T perform(const std::function<T(const std::map<std::string, Element> &)> &proc, const ParameterList &params, const T &failValue = T(), ErrorMessage* errMsg = nullptr) const {
            std::map<std::string, Element> args;

            std::vector<ArgInfo> signature = m_signatures[0];
            if (mapParamsToArgs(params, signature, &args)) {
                if (errMsg)
                    *errMsg = ErrorMessage();
                return proc(args);
            }
            if (errMsg)
                *errMsg = ErrorMessage("Parameters are invalid.");
            return failValue;
        };
    };
    
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
        std::vector<std::map<std::string, Element>> m_depthwiseVariables;
        uint32_t m_maxDepth;
    public:
        LocalVariables() : m_maxDepth(UINT32_MAX) { m_depthwiseVariables.emplace_back(); }
        void pushDepth() { m_depthwiseVariables.emplace_back(); }
        void popDepth() { m_depthwiseVariables.pop_back(); }
        bool exists(const std::string &key) const {
            uint32_t maxDepth = std::min(uint32_t(m_depthwiseVariables.size()), m_maxDepth) - 1;
            for (int i = maxDepth; i >=0 ; --i)
                if (m_depthwiseVariables[i].count(key) > 0)
                    return true;
            return false;
        }
        uint32_t depth() const { return (uint32_t)m_depthwiseVariables.size(); }
        
        void saveFrom(uint32_t depth) { m_maxDepth = depth; }
        void restore() { m_maxDepth = UINT32_MAX; }
        
        const Element &operator[](const std::string &key) const {
            uint32_t maxDepth = std::min(uint32_t(m_depthwiseVariables.size()), m_maxDepth) - 1;
            for (int i = maxDepth; i >=0 ; --i)
                if (m_depthwiseVariables[i].count(key) > 0)
                    return m_depthwiseVariables[i].at(key);
            return m_depthwiseVariables.back().at(key);
        }
        Element &operator[](const std::string &key) {
            uint32_t maxDepth = std::min(uint32_t(m_depthwiseVariables.size()), m_maxDepth) - 1;
            for (int i = maxDepth; i >=0 ; --i)
                if (m_depthwiseVariables[i].count(key) > 0)
                    return m_depthwiseVariables[i].at(key);
            return m_depthwiseVariables.back()[key];
        }
        const Element &at(const std::string &key) const {
            return (*this)[key];
        }
        Element &at(const std::string &key) {
            uint32_t maxDepth = std::min(uint32_t(m_depthwiseVariables.size()), m_maxDepth) - 1;
            for (int i = maxDepth; i >=0 ; --i)
                if (m_depthwiseVariables[i].count(key) > 0)
                    return m_depthwiseVariables[i].at(key);
            return m_depthwiseVariables.back().at(key);
        }
    };
    
    class StackVariables {
        std::vector<LocalVariables> m_stack;
    public:
        StackVariables() { m_stack.emplace_back(); }
        LocalVariables &current() { return m_stack.back(); }
        
        void push() { m_stack.emplace_back(); }
        void pop() { m_stack.pop_back(); }
        
        bool exists(const std::string &key) const {
            for (int i = int(m_stack.size() - 1); i >= 0; --i) {
                const LocalVariables &layer = m_stack[i];
                if (layer.exists(key))
                    return true;
            }
            return false;
        }
        
        const Element &operator[](const std::string &key) const {
            for (int i = int(m_stack.size() - 1); i >= 0; --i) {
                const LocalVariables &layer = m_stack[i];
                if (layer.exists(key))
                    return layer.at(key);
            }
            return m_stack.back().at(key);
        }
        const Element &at(const std::string &key) const {
            for (int i = int(m_stack.size() - 1); i >= 0; --i) {
                const LocalVariables &layer = m_stack[i];
                if (layer.exists(key))
                    return layer.at(key);
            }
            return m_stack.back().at(key);
        }
    };
    
    struct ExecuteContext {
        StackVariables stackVariables;
        Element returnValue;
        bool returnFlag;
        
        std::string absFileDirPath;
        SceneRef scene;
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
    
    class BlockStatement : public Statement {
        std::vector<StatementRef> m_statements;
    public:
        BlockStatement(const StatementsRef &statements);
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class IfElseStatement : public Statement {
        ExpressionRef m_condExpr;
        StatementRef m_trueBlockStmt;
        StatementRef m_falseBlockStmt;
    public:
        IfElseStatement(const ExpressionRef &condExpr, const StatementRef &trueBlockStmt, const StatementRef &falseBlockStmt = nullptr) :
        m_condExpr(condExpr), m_trueBlockStmt(trueBlockStmt), m_falseBlockStmt(falseBlockStmt) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class ForStatement : public Statement {
        ExpressionRef m_preExpr;
        ExpressionRef m_condExpr;
        ExpressionRef m_postExpr;
        StatementRef m_blockStmt;
    public:
        ForStatement(const ExpressionRef &preExpr, const ExpressionRef &condExpr, const ExpressionRef &postExpr, const StatementRef &blockStmt) :
        m_preExpr(preExpr), m_condExpr(condExpr), m_postExpr(postExpr), m_blockStmt(blockStmt) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class FunctionDefinitionStatement : public Statement {
        std::string m_funcName;
        ArgumentDefinitionVecRef m_argDefs;
        StatementRef m_blockStmt;
    public:
        FunctionDefinitionStatement(const std::string &funcName, const ArgumentDefinitionVecRef &argDefs, const StatementRef &blockStmt) :
        m_funcName(funcName), m_argDefs(argDefs), m_blockStmt(blockStmt) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class ReturnStatement : public Statement {
        ExpressionRef m_expr;
    public:
        ReturnStatement(const ExpressionRef &expr = nullptr) : m_expr(expr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class BinaryExpression : public Expression {
        ExpressionRef m_left;
        std::string m_op;
        ExpressionRef m_right;
    public:
        BinaryExpression(const ExpressionRef &left, const std::string &op, const ExpressionRef &right) : m_left(left), m_op(op), m_right(right) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const override;
    };
    
    class SubstitutionExpression : public Expression {
        std::string m_varName;
        std::string m_op;
        ExpressionRef m_right;
    public:
        SubstitutionExpression(const std::string &varName, const std::string &op, const ExpressionRef &right) :
        m_varName(varName), m_op(op), m_right(right) { }
        
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
    
    class FunctionCallSingleTerm : public SingleTerm {
        std::string m_funcID;
        ParameterVecRef m_args;
    public:
        FunctionCallSingleTerm(const std::string &funcID, const ParameterVecRef &args) : m_funcID(funcID), m_args(args) { }
        
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
    
    class ArgumentDefinition {
        std::string m_name;
        ExpressionRef m_defaultValueExpr;
    public:
        ArgumentDefinition(const std::string &name, const ExpressionRef &defaultValueExpr = nullptr) : m_name(name), m_defaultValueExpr(defaultValueExpr) { }
        
        bool perform(ExecuteContext &context, ErrorMessage* errMsg) const;
        void getArgInfo(ArgInfo* info) const;
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
