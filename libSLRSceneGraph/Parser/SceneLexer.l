%{
#include <cstdio>
#include <cstdlib>
#include <string>
#include "SceneParsingDriver.h"
#include "SceneParser.tab.hh"

// Work around an incompatibility in flex (at least versions
// 2.5.31 through 2.5.33): it generates code that does
// not conform to C89.  See Debian bug 333231
// <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=333231>.
#undef yywrap
#define yywrap() 1
%}

%option noyywrap nounput batch debug noinput

%x STRING
%x COMMENT

NONZERO         [1-9]
DIGIT           [[:digit:]]
FLOAT           ({NONZERO}{DIGIT}*\.?|0\.|\.{DIGIT}){DIGIT}*|0
EXPONENT        [eE](-|\+?){DIGIT}+

INTEGER         [\+-]?{NONZERO}{DIGIT}*|0
REALNUMBER      [\+-]?{FLOAT}{EXPONENT}?
ID              [_[:alpha:]][_[:alnum:]]*
OPERATORS		=|\*|,|:
PARENTHESES		\(|\)
WSPACE          [[:blank:]]+

%{
    // Code run each time a pattern is matched.
    #define YY_USER_ACTION \
    /* printf("b: %u.%u, e: %u.%u, l: %u\n", loc.begin.line, loc.begin.column, loc.end.line, loc.end.column, yyleng); */ \
    loc.columns(yyleng);
%}

%%

%{
    // Code run each time yylex is called.
    //printf("yylex()\n");
    using namespace SLRSceneGraph;
    location &loc = driver.currentLocation;
    loc.step();
%}

<<EOF>> {
    DPRINTF("End of File\n");
    return SceneParser::make_END(loc);
}

{INTEGER} {
    int32_t value = std::atoi(yytext);
    DPRINTF("Integer: %d\n", value);
    return SceneParser::make_INTEGER(value, loc);
}

{REALNUMBER} {
    double value = std::atof(yytext);
    DPRINTF("Real Number: %g\n", value);
    return SceneParser::make_REALNUMBER(value, loc);
}

Translate {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::Translate, loc);
}
RotateX {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::RotateX, loc);
}
RotateY {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::RotateY, loc);
}
RotateZ {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::RotateZ, loc);
}
Scale {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::Scale, loc);
}
Spectrum {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::Spectrum, loc);
}
CreateNode {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateNode, loc);
}
SetRenderer {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::SetRenderer, loc);
}
CreateSensor {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateSensor, loc);
}
CreatePerspectiveCamera {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreatePerspectiveCamera, loc);
}
LoadImage {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::LoadImage, loc);
}
SpectrumTexture {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::SpectrumTexture, loc);
}
SetEnvironment {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::SetEnvironment, loc);
}
Load3DModel {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::Load3DModel, loc);
}
SetTransform {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::SetTransform, loc);
}
AddChild {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::AddChild, loc);   
}
CreateMatte {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateMatte, loc);
}
CreateDiffuseEmitter {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateDiffuseEmitter, loc);
}
CreateEmitterSurfaceMaterial {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateEmitterSurfaceMaterial, loc);
}
CreateMesh {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_API(API::CreateMesh, loc);   
}

{ID} {
    std::string value = yytext;
    DPRINTF("ID: %s\n", value.c_str());
    return SceneParser::make_ID(value, loc);
}

{PARENTHESES} {
    char c = yytext[0];
    DPRINTF("Parenthesis: %c\n", c);
    if (c == '(')
        return SceneParser::make_LPAR(loc);
    else if (c == ')')
        return SceneParser::make_RPAR(loc);
}

{OPERATORS}	{
    char c = yytext[0];
    DPRINTF("Operator: %c\n", c);
    if (c == ',')
        return SceneParser::make_COMMA(loc);
    else if (c == '=')
        return SceneParser::make_SUBSTITUTION(loc);
    else if (c == '*')
        return SceneParser::make_MUL(loc);
    else if (c == ':')
        return SceneParser::make_COLON(loc);
}

{WSPACE} {
    loc.step();
}

\n {
    loc.lines(yyleng);
    loc.step();
}

"//".*$ {
    printf("single line comment(s): (%s)\n", yytext);
}

"/*" {
    BEGIN(COMMENT);
    yymore();
}
<COMMENT>"*/" {
    printf("Comment(s): (%s)\n", yytext);
    BEGIN(INITIAL);
}
<COMMENT>"*" {
    yymore();
}
<COMMENT>[^*]+ {
    yymore();
}

\" {
    BEGIN(STRING);
    yymore();
}
<STRING>\n {
    driver.error(loc, std::string("Irregal literal: ") + yytext);
    BEGIN(INITIAL);
}
<STRING>([^\"\n]|\\\"|\\\n)+ {
    loc.end = loc.begin + yyleng;
    yymore();
}
<STRING>\" {
    BEGIN(INITIAL);
    loc.end = loc.begin + yyleng;
    std::string value = yytext;
    value = value.substr(1, value.length() - 2);
    DPRINTF("String: %s\n", value.c_str());
    return SceneParser::make_STRING(value, loc);
}

. {
    driver.error(loc, std::string("Irregal character: (") + yytext + ")\n");
    return SceneParser::make_CHAR(yytext[0], loc);
}

%%

namespace SLRSceneGraph {
    void SceneParsingDriver::beginScan() {
        yy_flex_debug = traceScanning;
        if (file.empty() || file == "-")
            yyin = stdin;
        else if (!(yyin = fopen(file.c_str(), "r"))) {
            error("cannot open " + file + ": " + strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    void SceneParsingDriver::endScan() {
        fclose(yyin);
    }
}
