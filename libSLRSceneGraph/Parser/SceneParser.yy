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

    namespace SLRSceneGraph {
        class SceneParsingDriver;

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
    }
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
%token
    END 0 "end of file"
    LPAR "("
    RPAR ")"
    COMMA ","
    SUBSTITUTION "="
    MUL "*"
    COLON ":"
;
%token<char> CHAR
%token<API> API
%token<int32_t> INTEGER
%token<double> REALNUMBER
%token<std::string> STRING
%token<std::string> ID

%type<Element> expression
%type<Element> value
%type<Element> function_call
%type<ParameterListRef> arguments
%type<Parameter> argument

%printer { yyoutput << $$; } <*>;
%printer { printf(""); } <API>;
%printer { printf(""); } <Parameter>;
%printer { printf(""); } <ParameterListRef>;

%left ","
%left ":"
%left "="
%left "*"

%%

input:
/* empty */ |
input statement
;

statement:
expression {
    printf("statement done.\n");
} |
ID "=" expression {
    if (!$3) {
        error(@3, "The expression does not return a value.");
        YYERROR;
    }
    driver.variables[$1] = $3;
    DPRINTF("substitution statement done.\n");
} |
error {
    printf("Parsing aborted.\n");
    YYABORT;
}
;

expression:
value {
    $$ = $1;
} |
function_call {
    $$ = $1;
} |
expression "*" expression {
    if ($1.type != Type::Matrix || $3.type != Type::Matrix) {
        error(@3, "'*' operator is only valid for matrix multiplication.");
        YYERROR;
    }
    $$ = SLRSceneGraph::mulMatrix4x4($1, $3);
}
;

value:
INTEGER {
    // DSTMT(std::cout << "value Integer: " << @1 << std::endl;)
    $$ = Element(Type::Integer, createShared<int32_t>($1));
} |
REALNUMBER {
    // DSTMT(std::cout << "value RealNumber: " << @1 << std::endl;)
    $$ = Element(Type::RealNumber, createShared<double>($1));
} |
STRING {
    // DSTMT(std::cout << "value String: " << @1 << std::endl;)
    $$ = Element(Type::String, createShared<std::string>($1));
} |
"(" arguments ")" {
    // DSTMT(std::cout << "value Array: " << @1 << std::endl;)
    $$ = Element(Type::Tuple, $2);
} |
ID {
    // DSTMT(std::cout << "value ID: " << @1 << std::endl;)
    if (driver.variables.count($1) == 0) {
        error(@$, "undefined variable: " + $1);
        YYERROR;
    }
    else {
        $$ = driver.variables.at($1);
    }
}
;

function_call:
API "(" arguments ")" {
    // DSTMT(std::cout << "function_call ID: " << @1 << " Args: " << @3 << std::endl;)
    // DPRINTF("%u params (named: %u, unnamed: %u)\n", 
    //         $3->numParams(), $3->named.size(), $3->unnamed.size());
    DSTMT(
        for (auto it = $3->named.begin(); it != $3->named.end(); ++it)
            std::cout << it->first << ": " << it->second << std::endl;
        for (auto it = $3->unnamed.begin(); it != $3->unnamed.end(); ++it)
            std::cout << *it << std::endl;
        )

    ParameterList &params = *$3.get();
    ErrorMessage errMsg;
    $$ = Element();
    switch($1) {
    case API::Translate:
        printf("Call Translate\n");
        $$ = SLRSceneGraph::Translate(params, &errMsg);
        break;
    case API::RotateX:
        printf("Call RotateX\n");
        $$ = SLRSceneGraph::RotateX(params, &errMsg);
        break;
    case API::RotateY:
        printf("Call RotateY\n");
        $$ = SLRSceneGraph::RotateY(params, &errMsg);
        break;
    case API::RotateZ:
        printf("Call RotateZ\n");
        $$ = SLRSceneGraph::RotateZ(params, &errMsg);
        break;
    case API::Scale:
        printf("Call Scale\n");
        $$ = SLRSceneGraph::Scale(params, &errMsg);
        break;
    case API::Spectrum:
        printf("Call Spectrum\n");
        $$ = SLRSceneGraph::CreateSpectrum(params, &errMsg);
        break;
    case API::SpectrumTexture:
        printf("Call SpectrumTexture\n");
        $$ = SLRSceneGraph::CreateSpectrumTexture(params, &errMsg);
        break;
    case API::CreateMatte:
        printf("Call CreateMatte\n");
        $$ = SLRSceneGraph::CreateMatte(params, &errMsg);
        break;
    case API::CreateDiffuseEmitter:
        printf("Call CreateDiffuseEmitter\n");
        $$ = SLRSceneGraph::CreateDiffuseEmitter(params, &errMsg);
        break;
    case API::CreateEmitterSurfaceMaterial:
        printf("Call CreateEmitterSurfaceMaterial\n");
        $$ = SLRSceneGraph::CreateEmitterSurfaceMaterial(params, &errMsg);
        break;
    case API::CreateMesh:
        printf("Call CreateMesh\n");
        $$ = SLRSceneGraph::CreateMesh(params, &errMsg);
        break;
    case API::CreateNode:
        printf("Call CreateNode\n");
        $$ = SLRSceneGraph::CreateNode(params, &errMsg);
        break;
    case API::SetTransform:
        printf("Call SetTransform\n");
        $$ = SLRSceneGraph::SetTransform(params, &errMsg);
        break;
    case API::AddChild:
        printf("Call AddChild\n");
        $$ = SLRSceneGraph::AddChild(params, &errMsg);
        break;
    case API::SetRenderer:
        printf("Call SetRenderer\n");
        $$ = SLRSceneGraph::SetRenderer(*$3.get(), &driver.context, &errMsg);
        break;
    case API::SetRenderSettings:
        printf("Call SetRenderSettings\n");
        $$ = SLRSceneGraph::SetRenderSettings(*$3.get(), &driver.context, &errMsg);
        break;
    case API::CreatePerspectiveCamera:
        printf("Call CreatePerspectiveCamera\n");
        $$ = SLRSceneGraph::CreatePerspectiveCamera(params, &errMsg);
        break;
    case API::Load3DModel:
        printf("Call Load3DModel\n");
        $$ = SLRSceneGraph::Load3DModel(params, &errMsg);
        break;
    case API::LoadImage:
        printf("Call LoadImage\n");
        break;
    case API::SetEnvironment:
        printf("Call SetEnvironment\n");
        break;
    default:
        break;
    }

    if (errMsg.error) {
        error(@$, errMsg.message);
        YYERROR;
    }
}
;

arguments:
/* empty */ {
    DSTMT(std::cout << "empty: " << @$ << std::endl;)
    $$ = createShared<ParameterList>();
} |
argument {
    $$ = createShared<ParameterList>();
    $$->add($1);
} |
arguments "," argument {
    // DSTMT(std::cout << @$ << " (" << @1 << ", " << @3 << ")" << std::endl;)
    $$ = $1;
    $$->add($3);
}
;

argument: 
expression {
    $$ = Parameter("", $1);
} |
STRING ":" expression {
    $$ = Parameter($1, $3);
}
;

%%

namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
