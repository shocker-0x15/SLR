%skeleton "lalr1.cc"
%require "3.0.4"
%defines
%define api.namespace {SLRSceneGraph}
%define parser_class_name {SceneParser}
%define api.token.constructor
%define api.value.type variant
%define parse.assert

%define parse.trace
%define parse.error verbose

%param { SceneParsingDriver &driver }

%code requires {
    #include "../API.hpp"

    #if DEBUG
    #define DSTMT(stmt) stmt;
    #define DPRINTF(fmt, ...) printf(fmt, ##__VA_ARGS__)
    #else
    #define DSTMT(stmt)
    #define DPRINTF(fmt, ...)
    #endif
}

%locations
%initial-action {
    // Initialize the initial location.
    @$.begin.filename = @$.end.filename = &driver.file;
};

%code {
    #include "SceneParsingDriver.h"
}

%define api.token.prefix {TOKEN_}
// %leftなどでもトークンの宣言にはなるが、リテラル文字列との関連付けは%tokenで行う必要がある。
%token
    EOF 0 "end of file"
    L_PAREN "("
    R_PAREN ")"
    L_BRACE "{"
    R_BRACE "}"
    L_ANGLE "<"
    R_ANGLE ">"
    L_BRACK "["
    R_BRACK "]"

    PLUS "+"
    MINUS "-"
    AST "*"
    SLASH "/"
    PERC "%"
    EXC "!"
    L_ANGLE_EQ "<="
    R_ANGLE_EQ ">="
    EQ_EQ "=="
    EXC_EQ "!="
    AND_AND "&&"
    VBAR_VBAR "||"
    EQ "="
    PLUS_EQ "+="
    MINUS_EQ "-="
    AST_EQ "*="
    SLASH_EQ "/="
    PERC_EQ "%="
    PLUS_PLUS "++"
    MINUS_MINUS "--"

    COLON ":"
    COMMA ","
    SEMICOLON ";"
;
%token<char> CHAR
%token<bool> BOOL
%token<int32_t> INTEGER
%token<double> REALNUMBER
%token<std::string> STRING
%token<std::string> ID
%token IF ELSE FOR FUNCTION RETURN

%type<StatementsRef> Statements
%type<StatementRef> Statement
%type<ExpressionRef> Expression
%type<TermRef> Term
%type<SingleTermRef> SingleTerm
%type<ValueRef> Value
%type<ValueRef> ImmValue
%type<ValueRef> TupleValue
%type<ArgumentDefinitionRef> ArgumentDefinition
%type<ArgumentDefinitionVecRef> ArgumentDefinitions
%type<ParameterRef> Parameter
%type<ParameterVecRef> Elements
%type<ParameterVecRef> Arguments

%printer { /*yyoutput << $$;*/ } <*>;

%nonassoc SEMICOLON
%left COMMA
%right ":"
%right PREC_SUBST "=" "+=" "-=" "*=" "/=" "%="
%left PREC_LOGIC_OR "||"
%left PREC_LOGIC_AND "&&"
%left PREC_EQ_REL "==" "!="
%left PREC_INEQ_REL "<" ">" "<=" ">="
%left PREC_ADD "+" "-"
%left PREC_MUL "*" "/" "%"
%right PREC_PRE_INC "++" "--" "!"
%left PREC_POST_INC

%%

Statements:
Statement {
    $$ = createShared<std::vector<StatementRef>>();
    if ($1)
        $$->push_back($1);
    driver.statements = $$;
} | 
Statements Statement {
    $$ = $1;
    if ($2)
        $$->push_back($2); 
    driver.statements = $$;
}
;

Statement:
Expression ";" { $$ = $1; } |
"{" Statements "}" { $$ = createShared<BlockStatement>($2); } |
IF "(" Expression ")" Statement { $$ = createShared<IfElseStatement>($3, $5); } |
IF "(" Expression ")" Statement ELSE Statement { $$ = createShared<IfElseStatement>($3, $5, $7); } |
FOR "(" Expression ";" Expression ";" Expression ")" Statement { $$ = createShared<ForStatement>($3, $5, $7, $9); } |
FUNCTION ID "(" ArgumentDefinitions ")" Statement { $$ = createShared<FunctionDefinitionStatement>($2, $4, $6); } |
RETURN ";" { $$ = createShared<ReturnStatement>(); } |
RETURN Expression ";" { $$ = createShared<ReturnStatement>($2); } |
error {
    printf("Parsing aborted.\n");
    YYABORT;
}
;

Expression:
Term { $$ = $1; } |
Expression "+" Term { $$ = createShared<BinaryExpression>($1, "+", $3); } |
Expression "-" Term { $$ = createShared<BinaryExpression>($1, "-", $3); } |
Expression "<" Expression { $$ = createShared<BinaryExpression>($1, "<", $3); } |
Expression ">" Expression { $$ = createShared<BinaryExpression>($1, ">", $3); } |
Expression "<=" Expression { $$ = createShared<BinaryExpression>($1, "<=", $3); } |
Expression ">=" Expression { $$ = createShared<BinaryExpression>($1, ">=", $3); } |
Expression "==" Expression { $$ = createShared<BinaryExpression>($1, "==", $3); } |
Expression "!=" Expression { $$ = createShared<BinaryExpression>($1, "!=", $3); } |
Expression "&&" Expression { $$ = createShared<BinaryExpression>($1, "&&", $3); } |
Expression "||" Expression { $$ = createShared<BinaryExpression>($1, "||", $3); } |
ID "=" Expression { $$ = createShared<SubstitutionExpression>($1, "=", $3); } |
ID "+=" Expression { $$ = createShared<SubstitutionExpression>($1, "+=", $3); } |
ID "-=" Expression { $$ = createShared<SubstitutionExpression>($1, "-=", $3); } |
ID "*=" Expression { $$ = createShared<SubstitutionExpression>($1, "*=", $3); } |
ID "/=" Expression { $$ = createShared<SubstitutionExpression>($1, "/=", $3); } |
ID "%=" Expression { $$ = createShared<SubstitutionExpression>($1, "%=", $3); }
;

Term:
SingleTerm { $$ = $1; } | 
"+" SingleTerm %prec PREC_PRE_INC { $$ = createShared<UnaryTerm>("+", $2); } |
"-" SingleTerm %prec PREC_PRE_INC { $$ = createShared<UnaryTerm>("-", $2); } |
"!" SingleTerm { $$ = createShared<UnaryTerm>("!", $2); } |
"++" ID { $$ = createShared<UnarySubstitutionTerm>("++*", $2); } |
"--" ID { $$ = createShared<UnarySubstitutionTerm>("--*", $2); } |
ID "++" %prec PREC_POST_INC { $$ = createShared<UnarySubstitutionTerm>("*++", $1); } |
ID "--" %prec PREC_POST_INC { $$ = createShared<UnarySubstitutionTerm>("*--", $1); } |
Term "*" Term { $$ = createShared<BinaryTerm>($1, "*", $3); } |
Term "/" Term { $$ = createShared<BinaryTerm>($1, "/", $3); } |
Term "%" Term { $$ = createShared<BinaryTerm>($1, "%", $3); }
;

SingleTerm:
Value { $$ = $1; } |
ID "(" Arguments ")" { $$ = createShared<FunctionCallSingleTerm>($1, $3); } |
"(" Expression ")" { $$ = createShared<EnclosedSingleTerm>($2); } |
SingleTerm "[" Expression "]" { $$ = createShared<TupleElementSingleTerm>($1, $3); }
;

Value:
ImmValue { $$ = $1; } |
TupleValue { $$ = $1; } |
ID { $$ = createShared<VariableValue>($1); }
;

ImmValue:
BOOL { $$ = createShared<ImmediateValue>(Element::create<TypeMap::Bool>($1)); } |
INTEGER { $$ = createShared<ImmediateValue>(Element::create<TypeMap::Integer>($1)); } |
REALNUMBER { $$ = createShared<ImmediateValue>(Element::create<TypeMap::RealNumber>($1)); } |
STRING { $$ = createShared<ImmediateValue>(Element::create<TypeMap::String>($1)); }
;

TupleValue:
"(" "," ")" {
    $$ = createShared<TupleValue>(createShared<std::vector<ParameterRef>>());
} | 
"(" Parameter "," ")" {
    ParameterVecRef elem = createShared<std::vector<ParameterRef>>();
    elem->push_back($2);
    $$ = createShared<TupleValue>(elem);
} |
"(" Elements ")" {
    $$ = createShared<TupleValue>($2);
}
;

ArgumentDefinition:
ID { $$ = createShared<ArgumentDefinition>($1); } |
ID "=" Expression { $$ = createShared<ArgumentDefinition>($1, $3); }
;

ArgumentDefinitions:
/* empty */ {
    $$ = createShared<std::vector<ArgumentDefinitionRef>>();
} |
ArgumentDefinition {
    $$ = createShared<std::vector<ArgumentDefinitionRef>>();
    $$->push_back($1);
} |
ArgumentDefinitions "," ArgumentDefinition {
    $$ = $1;
    $$->push_back($3);
}
;

Parameter:
Expression { $$ = createShared<Parameter>(nullptr, $1); } |
Expression ":" Expression { $$ = createShared<Parameter>($1, $3); }
;

Elements:
Parameter "," Parameter {
    $$ = createShared<std::vector<ParameterRef>>();
    $$->push_back($1);
    $$->push_back($3);
} |
Elements "," Parameter {
    $$ = $1;
    $$->push_back($3);
}
;

Arguments:
/* empty */ {
    $$ = createShared<std::vector<ParameterRef>>();
} |
Parameter {
    $$ = createShared<std::vector<ParameterRef>>();
    $$->push_back($1);
} |
Arguments "," Parameter {
    $$ = $1;
    $$->push_back($3);
}
;

%%

namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
