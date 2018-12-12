%{
#include "SceneParsingDriver.h"
#include "SceneParser.tab.hh"
#ifdef SLR_Platform_Windows
#include <io.h>
#else
#include <unistd.h>
#endif

// Work around an incompatibility in flex (at least versions
// 2.5.31 through 2.5.33): it generates code that does
// not conform to C89.  See Debian bug 333231
// <http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=333231>.
#undef yywrap
#define yywrap() 1
%}

%option noyywrap nounput batch debug noinput nounistd

%x STRING
%x COMMENT

NONZERO         [1-9]
DIGIT           [[:digit:]]
FLOAT           ({NONZERO}{DIGIT}*\.?|0\.|\.{DIGIT}){DIGIT}*|0
EXPONENT        [eE](-|\+?){DIGIT}+

BOOL            True|False|true|false
INTEGER         {NONZERO}{DIGIT}*|0
REALNUMBER      {FLOAT}{EXPONENT}?
ID              [_[:alpha:]][_[:alnum:]]*
OPERATORS       <=|>=|<|>|==|!=|=|&&|\|\||\+=|-=|\*=|\/=|%=|\+\+|--|\+|-|\*|\/|%|!
DELIMITERS      ,|:|;
PARENTHESES     \(|\)|\{|\}|\[|\]
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
    return SceneParser::make_EOF(loc);
}

{BOOL} {
    std::string str = yytext;
    bool value = str == "True" || str == "true";
    DPRINTF("Bool: %s\n", value ? "true" : "false");
    return SceneParser::make_BOOL(value, loc);
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

if {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_IF(loc);
}

else {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_ELSE(loc);
}

for {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_FOR(loc);
}

function {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_FUNCTION(loc);
}

return {
    DPRINTF("%s\n", yytext);
    return SceneParser::make_RETURN(loc);
}

{ID} {
    std::string value = yytext;
    DPRINTF("ID: %s\n", value.c_str());
    return SceneParser::make_ID(value, loc);
}

{PARENTHESES} {
    char c = yytext[0];
    DPRINTF("Parenthesis: %c\n", c);
    if (c == '(') return SceneParser::make_L_PAREN(loc);
    if (c == ')') return SceneParser::make_R_PAREN(loc);
    if (c == '{') return SceneParser::make_L_BRACE(loc);
    if (c == '}') return SceneParser::make_R_BRACE(loc);
    if (c == '[') return SceneParser::make_L_BRACK(loc);
    if (c == ']') return SceneParser::make_R_BRACK(loc);
    SLRAssert(false, "Undefined parenthesis.");
}

{OPERATORS} {
    std::string str = yytext;
    DPRINTF("Operator: %s\n", str.c_str());
    if (str == "<=") return SceneParser::make_L_ANGLE_EQ(loc);
    if (str == ">=") return SceneParser::make_R_ANGLE_EQ(loc);
    if (str == "<") return SceneParser::make_L_ANGLE(loc);
    if (str == ">") return SceneParser::make_R_ANGLE(loc);
    if (str == "==") return SceneParser::make_EQ_EQ(loc);
    if (str == "!=") return SceneParser::make_EXC_EQ(loc);
    if (str == "&&") return SceneParser::make_AND_AND(loc);
    if (str == "||") return SceneParser::make_VBAR_VBAR(loc);
    if (str == "=") return SceneParser::make_EQ(loc);
    if (str == "+=") return SceneParser::make_PLUS_EQ(loc);
    if (str == "-=") return SceneParser::make_MINUS_EQ(loc);
    if (str == "*=") return SceneParser::make_AST_EQ(loc);
    if (str == "/=") return SceneParser::make_SLASH_EQ(loc);
    if (str == "\%=") return SceneParser::make_PERC_EQ(loc);
    if (str == "++") return SceneParser::make_PLUS_PLUS(loc);
    if (str == "--") return SceneParser::make_MINUS_MINUS(loc);
    if (str == "+") return SceneParser::make_PLUS(loc);
    if (str == "-") return SceneParser::make_MINUS(loc);
    if (str == "*") return SceneParser::make_AST(loc);
    if (str == "/") return SceneParser::make_SLASH(loc);
    if (str == "\%") return SceneParser::make_PERC(loc);
    if (str == "!") return SceneParser::make_EXC(loc);
    SLRAssert(false, "Undefined operator.");
}

{DELIMITERS} {
    char c = yytext[0];
    DPRINTF("Delimiter: %c\n", c);
    if (c == ',')
        return SceneParser::make_COMMA(loc);
    else if (c == ':')
        return SceneParser::make_COLON(loc);
    else if (c == ';')
        return SceneParser::make_SEMICOLON(loc);
    SLRAssert(false, "Undefined delimiter.");
}

{WSPACE} {
    loc.step();
}

\n {
    loc.lines(yyleng);
    loc.step();
}

"//".*$ {
    DPRINTF("single line comment(s): (%s)\n", yytext);
}

"/*" {
    BEGIN(COMMENT);
    yymore();
}
<COMMENT>"*/" {
    DPRINTF("Comment(s): (%s)\n", yytext);
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
