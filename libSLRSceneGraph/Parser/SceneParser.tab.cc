// A Bison parser, made by GNU Bison 3.0.4.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.


// First part of user declarations.

#line 37 "SceneParser.tab.cc" // lalr1.cc:404

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

#include "SceneParser.tab.hh"

// User implementation prologue.

#line 51 "SceneParser.tab.cc" // lalr1.cc:412
// Unqualified %code blocks.
#line 60 "SceneParser.yy" // lalr1.cc:413

    #include "SceneParsingDriver.h"

#line 57 "SceneParser.tab.cc" // lalr1.cc:413


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (/*CONSTCOND*/ false)
# endif


// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void) (E))

// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << std::endl;                  \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE(Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void>(0)
# define YY_STACK_PRINT()                static_cast<void>(0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 4 "SceneParser.yy" // lalr1.cc:479
namespace SLRSceneGraph {
#line 143 "SceneParser.tab.cc" // lalr1.cc:479

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  SceneParser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              // Fall through.
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }


  /// Build a parser object.
  SceneParser::SceneParser (SceneParsingDriver &driver_yyarg)
    :
#if YYDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      driver (driver_yyarg)
  {}

  SceneParser::~SceneParser ()
  {}


  /*---------------.
  | Symbol types.  |
  `---------------*/



  // by_state.
  inline
  SceneParser::by_state::by_state ()
    : state (empty_state)
  {}

  inline
  SceneParser::by_state::by_state (const by_state& other)
    : state (other.state)
  {}

  inline
  void
  SceneParser::by_state::clear ()
  {
    state = empty_state;
  }

  inline
  void
  SceneParser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  inline
  SceneParser::by_state::by_state (state_type s)
    : state (s)
  {}

  inline
  SceneParser::symbol_number_type
  SceneParser::by_state::type_get () const
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  inline
  SceneParser::stack_symbol_type::stack_symbol_type ()
  {}


  inline
  SceneParser::stack_symbol_type::stack_symbol_type (state_type s, symbol_type& that)
    : super_type (s, that.location)
  {
      switch (that.type_get ())
    {
      case 10: // API
        value.move< API > (that.value);
        break;

      case 18: // expression
      case 19: // value
      case 20: // function_call
        value.move< Element > (that.value);
        break;

      case 22: // argument
        value.move< Parameter > (that.value);
        break;

      case 21: // arguments
        value.move< ParameterListRef > (that.value);
        break;

      case 9: // CHAR
        value.move< char > (that.value);
        break;

      case 12: // REALNUMBER
        value.move< double > (that.value);
        break;

      case 11: // INTEGER
        value.move< int32_t > (that.value);
        break;

      case 13: // STRING
      case 14: // ID
        value.move< std::string > (that.value);
        break;

      default:
        break;
    }

    // that is emptied.
    that.type = empty_symbol;
  }

  inline
  SceneParser::stack_symbol_type&
  SceneParser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
      switch (that.type_get ())
    {
      case 10: // API
        value.copy< API > (that.value);
        break;

      case 18: // expression
      case 19: // value
      case 20: // function_call
        value.copy< Element > (that.value);
        break;

      case 22: // argument
        value.copy< Parameter > (that.value);
        break;

      case 21: // arguments
        value.copy< ParameterListRef > (that.value);
        break;

      case 9: // CHAR
        value.copy< char > (that.value);
        break;

      case 12: // REALNUMBER
        value.copy< double > (that.value);
        break;

      case 11: // INTEGER
        value.copy< int32_t > (that.value);
        break;

      case 13: // STRING
      case 14: // ID
        value.copy< std::string > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }


  template <typename Base>
  inline
  void
  SceneParser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  SceneParser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    switch (yytype)
    {
            case 9: // CHAR

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< char > (); }
#line 380 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 10: // API

#line 88 "SceneParser.yy" // lalr1.cc:636
        { printf(""); }
#line 387 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 11: // INTEGER

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< int32_t > (); }
#line 394 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 12: // REALNUMBER

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< double > (); }
#line 401 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 13: // STRING

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< std::string > (); }
#line 408 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 14: // ID

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< std::string > (); }
#line 415 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 18: // expression

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< Element > (); }
#line 422 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 19: // value

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< Element > (); }
#line 429 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 20: // function_call

#line 87 "SceneParser.yy" // lalr1.cc:636
        { yyoutput << yysym.value.template as< Element > (); }
#line 436 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 21: // arguments

#line 90 "SceneParser.yy" // lalr1.cc:636
        { printf(""); }
#line 443 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 22: // argument

#line 89 "SceneParser.yy" // lalr1.cc:636
        { printf(""); }
#line 450 "SceneParser.tab.cc" // lalr1.cc:636
        break;


      default:
        break;
    }
    yyo << ')';
  }
#endif

  inline
  void
  SceneParser::yypush_ (const char* m, state_type s, symbol_type& sym)
  {
    stack_symbol_type t (s, sym);
    yypush_ (m, t);
  }

  inline
  void
  SceneParser::yypush_ (const char* m, stack_symbol_type& s)
  {
    if (m)
      YY_SYMBOL_PRINT (m, s);
    yystack_.push (s);
  }

  inline
  void
  SceneParser::yypop_ (unsigned int n)
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  SceneParser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  SceneParser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  SceneParser::debug_level_type
  SceneParser::debug_level () const
  {
    return yydebug_;
  }

  void
  SceneParser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  inline SceneParser::state_type
  SceneParser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  inline bool
  SceneParser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  inline bool
  SceneParser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  SceneParser::parse ()
  {
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << std::endl;


    // User initialization code.
    #line 55 "SceneParser.yy" // lalr1.cc:745
{
    // Initialize the initial location.
    yyla.location.begin.filename = yyla.location.end.filename = &driver.file;
}

#line 569 "SceneParser.tab.cc" // lalr1.cc:745

    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, yyla);

    // A new symbol was pushed on the stack.
  yynewstate:
    YYCDEBUG << "Entering state " << yystack_[0].state << std::endl;

    // Accept?
    if (yystack_[0].state == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    // Backup.
  yybackup:

    // Try to take a decision without lookahead.
    yyn = yypact_[yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
        try
          {
            symbol_type yylookahead (yylex (driver));
            yyla.move (yylookahead);
          }
        catch (const syntax_error& yyexc)
          {
            error (yyexc);
            goto yyerrlab1;
          }
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      goto yydefault;

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", yyn, yyla);
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
        switch (yyr1_[yyn])
    {
      case 10: // API
        yylhs.value.build< API > ();
        break;

      case 18: // expression
      case 19: // value
      case 20: // function_call
        yylhs.value.build< Element > ();
        break;

      case 22: // argument
        yylhs.value.build< Parameter > ();
        break;

      case 21: // arguments
        yylhs.value.build< ParameterListRef > ();
        break;

      case 9: // CHAR
        yylhs.value.build< char > ();
        break;

      case 12: // REALNUMBER
        yylhs.value.build< double > ();
        break;

      case 11: // INTEGER
        yylhs.value.build< int32_t > ();
        break;

      case 13: // STRING
      case 14: // ID
        yylhs.value.build< std::string > ();
        break;

      default:
        break;
    }


      // Compute the default @$.
      {
        slice<stack_symbol_type, stack_type> slice (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, slice, yylen);
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
      try
        {
          switch (yyn)
            {
  case 4:
#line 105 "SceneParser.yy" // lalr1.cc:859
    {
    printf("statement done.\n");
}
#line 716 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 5:
#line 108 "SceneParser.yy" // lalr1.cc:859
    {
    if (!yystack_[0].value.as< Element > ()) {
        error(yystack_[0].location, "The expression does not return a value.");
        YYERROR;
    }
    driver.variables[yystack_[2].value.as< std::string > ()] = yystack_[0].value.as< Element > ();
    DPRINTF("substitution statement done.\n");
}
#line 729 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 6:
#line 116 "SceneParser.yy" // lalr1.cc:859
    {
    printf("Parsing aborted.\n");
    YYABORT;
}
#line 738 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 7:
#line 123 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< Element > () = yystack_[0].value.as< Element > ();
}
#line 746 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 8:
#line 126 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< Element > () = yystack_[0].value.as< Element > ();
}
#line 754 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 9:
#line 129 "SceneParser.yy" // lalr1.cc:859
    {
    if (yystack_[2].value.as< Element > ().type != Type::Matrix || yystack_[0].value.as< Element > ().type != Type::Matrix) {
        error(yystack_[0].location, "'*' operator is only valid for matrix multiplication.");
        YYERROR;
    }
    yylhs.value.as< Element > () = SLRSceneGraph::mulMatrix4x4(yystack_[2].value.as< Element > (), yystack_[0].value.as< Element > ());
}
#line 766 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 10:
#line 139 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "value Integer: " << @1 << std::endl;)
    yylhs.value.as< Element > () = Element(Type::Integer, createShared<int32_t>(yystack_[0].value.as< int32_t > ()));
}
#line 775 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 11:
#line 143 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "value RealNumber: " << @1 << std::endl;)
    yylhs.value.as< Element > () = Element(Type::RealNumber, createShared<double>(yystack_[0].value.as< double > ()));
}
#line 784 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 12:
#line 147 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "value String: " << @1 << std::endl;)
    yylhs.value.as< Element > () = Element(Type::String, createShared<std::string>(yystack_[0].value.as< std::string > ()));
}
#line 793 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 13:
#line 151 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "value Array: " << @1 << std::endl;)
    yylhs.value.as< Element > () = Element(Type::Tuple, yystack_[1].value.as< ParameterListRef > ());
}
#line 802 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 14:
#line 155 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "value ID: " << @1 << std::endl;)
    if (driver.variables.count(yystack_[0].value.as< std::string > ()) == 0) {
        error(yylhs.location, "undefined variable: " + yystack_[0].value.as< std::string > ());
        YYERROR;
    }
    else {
        yylhs.value.as< Element > () = driver.variables.at(yystack_[0].value.as< std::string > ());
    }
}
#line 817 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 15:
#line 168 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << "function_call ID: " << @1 << " Args: " << @3 << std::endl;)
    // DPRINTF("%u params (named: %u, unnamed: %u)\n", 
    //         $3->numParams(), $3->named.size(), $3->unnamed.size());
    DSTMT(
        for (auto it = yystack_[1].value.as< ParameterListRef > ()->named.begin(); it != yystack_[1].value.as< ParameterListRef > ()->named.end(); ++it)
            std::cout << it->first << ": " << it->second << std::endl;
        for (auto it = yystack_[1].value.as< ParameterListRef > ()->unnamed.begin(); it != yystack_[1].value.as< ParameterListRef > ()->unnamed.end(); ++it)
            std::cout << *it << std::endl;
        )

    ParameterList &params = *yystack_[1].value.as< ParameterListRef > ().get();
    ErrorMessage errMsg;
    yylhs.value.as< Element > () = Element();
    switch(yystack_[3].value.as< API > ()) {
    case API::Translate:
        printf("Call Translate\n");
        yylhs.value.as< Element > () = SLRSceneGraph::Translate(params, &errMsg);
        break;
    case API::RotateX:
        printf("Call RotateX\n");
        yylhs.value.as< Element > () = SLRSceneGraph::RotateX(params, &errMsg);
        break;
    case API::RotateY:
        printf("Call RotateY\n");
        yylhs.value.as< Element > () = SLRSceneGraph::RotateY(params, &errMsg);
        break;
    case API::RotateZ:
        printf("Call RotateZ\n");
        yylhs.value.as< Element > () = SLRSceneGraph::RotateZ(params, &errMsg);
        break;
    case API::Scale:
        printf("Call Scale\n");
        yylhs.value.as< Element > () = SLRSceneGraph::Scale(params, &errMsg);
        break;
    case API::Spectrum:
        printf("Call Spectrum\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateSpectrum(params, &errMsg);
        break;
    case API::SpectrumTexture:
        printf("Call SpectrumTexture\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateSpectrumTexture(params, &errMsg);
        break;
    case API::CreateMatte:
        printf("Call CreateMatte\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateMatte(params, &errMsg);
        break;
    case API::CreateDiffuseEmitter:
        printf("Call CreateDiffuseEmitter\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateDiffuseEmitter(params, &errMsg);
        break;
    case API::CreateEmitterSurfaceMaterial:
        printf("Call CreateEmitterSurfaceMaterial\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateEmitterSurfaceMaterial(params, &errMsg);
        break;
    case API::CreateMesh:
        printf("Call CreateMesh\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateMesh(params, &errMsg);
        break;
    case API::CreateNode:
        printf("Call CreateNode\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreateNode(params, &errMsg);
        break;
    case API::SetTransform:
        printf("Call SetTransform\n");
        yylhs.value.as< Element > () = SLRSceneGraph::SetTransform(params, &errMsg);
        break;
    case API::AddChild:
        printf("Call AddChild\n");
        yylhs.value.as< Element > () = SLRSceneGraph::AddChild(params, &errMsg);
        break;
    case API::SetRenderer:
        printf("Call SetRenderer\n");
        yylhs.value.as< Element > () = SLRSceneGraph::SetRenderer(*yystack_[1].value.as< ParameterListRef > ().get(), &driver.context, &errMsg);
        break;
    case API::SetRenderSettings:
        printf("Call SetRenderSettings\n");
        yylhs.value.as< Element > () = SLRSceneGraph::SetRenderSettings(*yystack_[1].value.as< ParameterListRef > ().get(), &driver.context, &errMsg);
        break;
    case API::CreatePerspectiveCamera:
        printf("Call CreatePerspectiveCamera\n");
        yylhs.value.as< Element > () = SLRSceneGraph::CreatePerspectiveCamera(params, &errMsg);
        break;
    case API::Load3DModel:
        printf("Call Load3DModel\n");
        yylhs.value.as< Element > () = SLRSceneGraph::Load3DModel(params, &errMsg);
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
        error(yylhs.location, errMsg.message);
        YYERROR;
    }
}
#line 924 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 16:
#line 273 "SceneParser.yy" // lalr1.cc:859
    {
    DSTMT(std::cout << "empty: " << yylhs.location << std::endl;)
    yylhs.value.as< ParameterListRef > () = createShared<ParameterList>();
}
#line 933 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 17:
#line 277 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterListRef > () = createShared<ParameterList>();
    yylhs.value.as< ParameterListRef > ()->add(yystack_[0].value.as< Parameter > ());
}
#line 942 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 18:
#line 281 "SceneParser.yy" // lalr1.cc:859
    {
    // DSTMT(std::cout << @$ << " (" << @1 << ", " << @3 << ")" << std::endl;)
    yylhs.value.as< ParameterListRef > () = yystack_[2].value.as< ParameterListRef > ();
    yylhs.value.as< ParameterListRef > ()->add(yystack_[0].value.as< Parameter > ());
}
#line 952 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 19:
#line 289 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< Parameter > () = Parameter("", yystack_[0].value.as< Element > ());
}
#line 960 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 20:
#line 292 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< Parameter > () = Parameter(yystack_[2].value.as< std::string > (), yystack_[0].value.as< Element > ());
}
#line 968 "SceneParser.tab.cc" // lalr1.cc:859
    break;


#line 972 "SceneParser.tab.cc" // lalr1.cc:859
            default:
              break;
            }
        }
      catch (const syntax_error& yyexc)
        {
          error (yyexc);
          YYERROR;
        }
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, yylhs);
    }
    goto yynewstate;

  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;
    yyerror_range[1].location = yystack_[yylen - 1].location;
    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yyterror_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = yyn;
      yypush_ ("Shifting", error_token);
    }
    goto yynewstate;

    // Accept.
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    // Abort.
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << std::endl;
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  void
  SceneParser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what());
  }

  // Generate an error message.
  std::string
  SceneParser::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    size_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state
         merging (from LALR or IELR) and default reductions corrupt the
         expected token list.  However, the list is correct for
         canonical LR with one exception: it will still contain any
         token that will not be accepted due to an error action in a
         later state.
    */
    if (!yyla.empty ())
      {
        int yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];
        int yyn = yypact_[yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yyterror_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
        YYCASE_(0, YY_("syntax error"));
        YYCASE_(1, YY_("syntax error, unexpected %s"));
        YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    size_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char SceneParser::yypact_ninf_ = -10;

  const signed char SceneParser::yytable_ninf_ = -1;

  const signed char
  SceneParser::yypact_[] =
  {
     -10,     1,   -10,   -10,    13,     0,   -10,   -10,   -10,    11,
     -10,    -2,   -10,   -10,    10,   -10,    -2,     2,   -10,    13,
      19,    19,    19,   -10,    13,     4,    -2,   -10,    -2,   -10,
     -10
  };

  const unsigned char
  SceneParser::yydefact_[] =
  {
       2,     0,     1,     6,    16,     0,    10,    11,    12,    14,
       3,     4,     7,     8,    12,    14,    19,     0,    17,    16,
       0,     0,     0,    13,     0,     0,     5,     9,    20,    18,
      15
  };

  const signed char
  SceneParser::yypgoto_[] =
  {
     -10,   -10,   -10,    -1,   -10,   -10,    -9,    12
  };

  const signed char
  SceneParser::yydefgoto_[] =
  {
      -1,     1,    10,    16,    12,    13,    17,    18
  };

  const unsigned char
  SceneParser::yytable_[] =
  {
      11,     2,     3,    19,     4,    21,    23,    24,    30,    24,
      25,     5,     6,     7,     8,     9,     4,    20,    22,    26,
      27,    28,     4,     5,     6,     7,    14,    15,     0,     5,
       6,     7,     8,    15,     0,     0,    29
  };

  const signed char
  SceneParser::yycheck_[] =
  {
       1,     0,     1,     3,     3,     7,     4,     5,     4,     5,
      19,    10,    11,    12,    13,    14,     3,     6,     8,    20,
      21,    22,     3,    10,    11,    12,    13,    14,    -1,    10,
      11,    12,    13,    14,    -1,    -1,    24
  };

  const unsigned char
  SceneParser::yystos_[] =
  {
       0,    16,     0,     1,     3,    10,    11,    12,    13,    14,
      17,    18,    19,    20,    13,    14,    18,    21,    22,     3,
       6,     7,     8,     4,     5,    21,    18,    18,    18,    22,
       4
  };

  const unsigned char
  SceneParser::yyr1_[] =
  {
       0,    15,    16,    16,    17,    17,    17,    18,    18,    18,
      19,    19,    19,    19,    19,    20,    21,    21,    21,    22,
      22
  };

  const unsigned char
  SceneParser::yyr2_[] =
  {
       0,     2,     0,     2,     1,     3,     1,     1,     1,     3,
       1,     1,     1,     3,     1,     4,     0,     1,     3,     1,
       3
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const SceneParser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "\"(\"", "\")\"", "\",\"",
  "\"=\"", "\"*\"", "\":\"", "CHAR", "API", "INTEGER", "REALNUMBER",
  "STRING", "ID", "$accept", "input", "statement", "expression", "value",
  "function_call", "arguments", "argument", YY_NULLPTR
  };

#if YYDEBUG
  const unsigned short int
  SceneParser::yyrline_[] =
  {
       0,    99,    99,   101,   105,   108,   116,   123,   126,   129,
     139,   143,   147,   151,   155,   168,   273,   277,   281,   289,
     292
  };

  // Print the state stack on the debug stream.
  void
  SceneParser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << i->state;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  SceneParser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):" << std::endl;
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG


#line 4 "SceneParser.yy" // lalr1.cc:1167
} // SLRSceneGraph
#line 1358 "SceneParser.tab.cc" // lalr1.cc:1167
#line 297 "SceneParser.yy" // lalr1.cc:1168


namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
