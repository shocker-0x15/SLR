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

    #if 1
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
    LPAR "("
    RPAR ")"
    LBRC "{"
    RBRC "}"
    COMMA ","
    COLON ":"
    SEMICOLON ";"
;
%token<char> CHAR
%token<API> API
%token<bool> BOOL
%token<int32_t> INTEGER
%token<double> REALNUMBER
%token<std::string> STRING
%token<std::string> ID
%token FOR

%type<StatementsRef> Statements
%type<StatementRef> Statement
%type<ExpressionRef> Expression
%type<TermRef> Term
%type<ValueRef> Value
%type<ValueRef> ImmValue
%type<ValueRef> TupleValue
%type<ArgumentsRef> Elements
%type<TermRef> Function
%type<ArgumentsRef> Arguments
%type<ArgumentRef> Argument

%printer { /*yyoutput << $$;*/ } <*>;

%left<char> COMMA
%nonassoc<char> COLON
%nonassoc<char> SEMICOLON
%right<std::string> SUBSTITUTION
%left<std::string> OPERATOR0
%left<std::string> OPERATOR1
%nonassoc NEG
%nonassoc<std::string> INC_DEC

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
Expression ";" {
    $$ = $1;
} |
FOR "(" Expression ";" Expression ";" Expression ")" "{" Statements "}" {
    $$ = createShared<ForStatement>($3, $5, $7, $10);
} |
error {
    printf("Parsing aborted.\n");
    YYABORT;
}
;

Expression:
Term {
    $$ = $1;
} |
Term OPERATOR0 Term {
    $$ = createShared<BinaryExpression>($1, $2, $3);
} |
ID SUBSTITUTION Expression {
    $$ = createShared<SubstitutionExpression>($1, $2, $3);
}
;

Term:
Value {
    $$ = $1;
} |
Function {
    $$ = $1;
} |
OPERATOR0 Term %prec NEG {
    $$ = createShared<UnaryTerm>($1, $2);
} |
INC_DEC ID {
    
} |
ID INC_DEC {

} |
Term OPERATOR1 Term {
    $$ = createShared<BinaryTerm>($1, $2, $3);
} |
"(" Expression ")" {
    $$ = createShared<EnclosedTerm>($2);
}
;

Value:
ImmValue {
    $$ = $1;
} |
TupleValue {
    $$ = $1;
} |
ID {
    $$ = createShared<VariableValue>($1);
}
;

ImmValue:
BOOL {
    $$ = createShared<ImmediateValue>(Element($1));
} |
INTEGER {
    $$ = createShared<ImmediateValue>(Element($1));
} |
REALNUMBER {
    $$ = createShared<ImmediateValue>(Element($1));
} |
STRING {
    $$ = createShared<ImmediateValue>(Element($1));
}
;

TupleValue:
"(" "," ")" {
    $$ = createShared<TupleValue>(createShared<std::vector<ArgumentRef>>());
} | 
"(" Argument "," ")" {
    ArgumentsRef elem = createShared<std::vector<ArgumentRef>>();
    elem->push_back($2);
    $$ = createShared<TupleValue>(elem);
} |
"(" Elements ")" {
    $$ = createShared<TupleValue>($2);
}
;

Elements:
Argument "," Argument {
    $$ = createShared<std::vector<ArgumentRef>>();
    $$->push_back($1);
    $$->push_back($3);
} |
Elements "," Argument {
    $$ = $1;
    $$->push_back($3);
}
;

Function:
API "(" Arguments ")" {
    $$ = createShared<FunctionTerm>($1, $3);
}
;

Arguments:
/* empty */ {
    $$ = createShared<std::vector<ArgumentRef>>();
} |
Argument {
    $$ = createShared<std::vector<ArgumentRef>>();
    $$->push_back($1);
} |
Arguments "," Argument {
    $$ = $1;
    $$->push_back($3);
}
;

Argument: 
Expression {
    $$ = createShared<Argument>(nullptr, $1);
} |
Expression ":" Expression {
    $$ = createShared<Argument>($1, $3);
}
;

%%

namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
