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
#line 33 "SceneParser.yy" // lalr1.cc:413

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
      case 63: // ArgumentDefinition
        value.move< ArgumentDefinitionRef > (that.value);
        break;

      case 64: // ArgumentDefinitions
        value.move< ArgumentDefinitionVecRef > (that.value);
        break;

      case 57: // Expression
        value.move< ExpressionRef > (that.value);
        break;

      case 65: // Parameter
        value.move< ParameterRef > (that.value);
        break;

      case 66: // Elements
      case 67: // Arguments
        value.move< ParameterVecRef > (that.value);
        break;

      case 59: // SingleTerm
        value.move< SingleTermRef > (that.value);
        break;

      case 56: // Statement
        value.move< StatementRef > (that.value);
        break;

      case 55: // Statements
        value.move< StatementsRef > (that.value);
        break;

      case 58: // Term
        value.move< TermRef > (that.value);
        break;

      case 60: // Value
      case 61: // ImmValue
      case 62: // TupleValue
        value.move< ValueRef > (that.value);
        break;

      case 35: // BOOL
        value.move< bool > (that.value);
        break;

      case 34: // CHAR
        value.move< char > (that.value);
        break;

      case 37: // REALNUMBER
        value.move< double > (that.value);
        break;

      case 36: // INTEGER
        value.move< int32_t > (that.value);
        break;

      case 38: // STRING
      case 39: // ID
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
      case 63: // ArgumentDefinition
        value.copy< ArgumentDefinitionRef > (that.value);
        break;

      case 64: // ArgumentDefinitions
        value.copy< ArgumentDefinitionVecRef > (that.value);
        break;

      case 57: // Expression
        value.copy< ExpressionRef > (that.value);
        break;

      case 65: // Parameter
        value.copy< ParameterRef > (that.value);
        break;

      case 66: // Elements
      case 67: // Arguments
        value.copy< ParameterVecRef > (that.value);
        break;

      case 59: // SingleTerm
        value.copy< SingleTermRef > (that.value);
        break;

      case 56: // Statement
        value.copy< StatementRef > (that.value);
        break;

      case 55: // Statements
        value.copy< StatementsRef > (that.value);
        break;

      case 58: // Term
        value.copy< TermRef > (that.value);
        break;

      case 60: // Value
      case 61: // ImmValue
      case 62: // TupleValue
        value.copy< ValueRef > (that.value);
        break;

      case 35: // BOOL
        value.copy< bool > (that.value);
        break;

      case 34: // CHAR
        value.copy< char > (that.value);
        break;

      case 37: // REALNUMBER
        value.copy< double > (that.value);
        break;

      case 36: // INTEGER
        value.copy< int32_t > (that.value);
        break;

      case 38: // STRING
      case 39: // ID
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
            case 34: // CHAR

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 438 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 35: // BOOL

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 445 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 36: // INTEGER

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 452 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 37: // REALNUMBER

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 459 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 38: // STRING

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 466 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 39: // ID

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 473 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 55: // Statements

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 480 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 56: // Statement

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 487 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 57: // Expression

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 494 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 58: // Term

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 501 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 59: // SingleTerm

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 508 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 60: // Value

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 515 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 61: // ImmValue

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 522 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 62: // TupleValue

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 529 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 63: // ArgumentDefinition

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 536 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 64: // ArgumentDefinitions

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 543 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 65: // Parameter

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 550 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 66: // Elements

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 557 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 67: // Arguments

#line 97 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 564 "SceneParser.tab.cc" // lalr1.cc:636
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
    #line 28 "SceneParser.yy" // lalr1.cc:745
{
    // Initialize the initial location.
    yyla.location.begin.filename = yyla.location.end.filename = &driver.file;
}

#line 683 "SceneParser.tab.cc" // lalr1.cc:745

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
      case 63: // ArgumentDefinition
        yylhs.value.build< ArgumentDefinitionRef > ();
        break;

      case 64: // ArgumentDefinitions
        yylhs.value.build< ArgumentDefinitionVecRef > ();
        break;

      case 57: // Expression
        yylhs.value.build< ExpressionRef > ();
        break;

      case 65: // Parameter
        yylhs.value.build< ParameterRef > ();
        break;

      case 66: // Elements
      case 67: // Arguments
        yylhs.value.build< ParameterVecRef > ();
        break;

      case 59: // SingleTerm
        yylhs.value.build< SingleTermRef > ();
        break;

      case 56: // Statement
        yylhs.value.build< StatementRef > ();
        break;

      case 55: // Statements
        yylhs.value.build< StatementsRef > ();
        break;

      case 58: // Term
        yylhs.value.build< TermRef > ();
        break;

      case 60: // Value
      case 61: // ImmValue
      case 62: // TupleValue
        yylhs.value.build< ValueRef > ();
        break;

      case 35: // BOOL
        yylhs.value.build< bool > ();
        break;

      case 34: // CHAR
        yylhs.value.build< char > ();
        break;

      case 37: // REALNUMBER
        yylhs.value.build< double > ();
        break;

      case 36: // INTEGER
        yylhs.value.build< int32_t > ();
        break;

      case 38: // STRING
      case 39: // ID
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
  case 2:
#line 115 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementsRef > () = createShared<std::vector<StatementRef>>();
    if (yystack_[0].value.as< StatementRef > ())
        yylhs.value.as< StatementsRef > ()->push_back(yystack_[0].value.as< StatementRef > ());
    driver.statements = yylhs.value.as< StatementsRef > ();
}
#line 862 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 3:
#line 121 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementsRef > () = yystack_[1].value.as< StatementsRef > ();
    if (yystack_[0].value.as< StatementRef > ())
        yylhs.value.as< StatementsRef > ()->push_back(yystack_[0].value.as< StatementRef > ()); 
    driver.statements = yylhs.value.as< StatementsRef > ();
}
#line 873 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 4:
#line 130 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = yystack_[1].value.as< ExpressionRef > (); }
#line 879 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 5:
#line 131 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<BlockStatement>(yystack_[1].value.as< StatementsRef > ()); }
#line 885 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 6:
#line 132 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<IfElseStatement>(yystack_[2].value.as< ExpressionRef > (), yystack_[0].value.as< StatementRef > ()); }
#line 891 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 7:
#line 133 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<IfElseStatement>(yystack_[4].value.as< ExpressionRef > (), yystack_[2].value.as< StatementRef > (), yystack_[0].value.as< StatementRef > ()); }
#line 897 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 8:
#line 134 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<ForStatement>(yystack_[6].value.as< ExpressionRef > (), yystack_[4].value.as< ExpressionRef > (), yystack_[2].value.as< ExpressionRef > (), yystack_[0].value.as< StatementRef > ()); }
#line 903 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 9:
#line 135 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<FunctionDefinitionStatement>(yystack_[4].value.as< std::string > (), yystack_[2].value.as< ArgumentDefinitionVecRef > (), yystack_[0].value.as< StatementRef > ()); }
#line 909 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 10:
#line 136 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<ReturnStatement>(); }
#line 915 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 11:
#line 137 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< StatementRef > () = createShared<ReturnStatement>(yystack_[1].value.as< ExpressionRef > ()); }
#line 921 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 12:
#line 138 "SceneParser.yy" // lalr1.cc:859
    {
    printf("Parsing aborted.\n");
    YYABORT;
}
#line 930 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 13:
#line 145 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = yystack_[0].value.as< TermRef > (); }
#line 936 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 14:
#line 146 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "+", yystack_[0].value.as< TermRef > ()); }
#line 942 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 15:
#line 147 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "-", yystack_[0].value.as< TermRef > ()); }
#line 948 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 16:
#line 148 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "<", yystack_[0].value.as< ExpressionRef > ()); }
#line 954 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 17:
#line 149 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), ">", yystack_[0].value.as< ExpressionRef > ()); }
#line 960 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 18:
#line 150 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "<=", yystack_[0].value.as< ExpressionRef > ()); }
#line 966 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 19:
#line 151 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), ">=", yystack_[0].value.as< ExpressionRef > ()); }
#line 972 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 20:
#line 152 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "==", yystack_[0].value.as< ExpressionRef > ()); }
#line 978 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 21:
#line 153 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "!=", yystack_[0].value.as< ExpressionRef > ()); }
#line 984 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 22:
#line 154 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "&&", yystack_[0].value.as< ExpressionRef > ()); }
#line 990 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 23:
#line 155 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< ExpressionRef > (), "||", yystack_[0].value.as< ExpressionRef > ()); }
#line 996 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 24:
#line 156 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1002 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 25:
#line 157 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "+=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1008 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 26:
#line 158 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "-=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1014 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 27:
#line 159 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "*=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1020 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 28:
#line 160 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "/=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1026 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 29:
#line 161 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), "%=", yystack_[0].value.as< ExpressionRef > ()); }
#line 1032 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 30:
#line 165 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = yystack_[0].value.as< SingleTermRef > (); }
#line 1038 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 31:
#line 166 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnaryTerm>("+", yystack_[0].value.as< SingleTermRef > ()); }
#line 1044 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 32:
#line 167 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnaryTerm>("-", yystack_[0].value.as< SingleTermRef > ()); }
#line 1050 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 33:
#line 168 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnaryTerm>("!", yystack_[0].value.as< SingleTermRef > ()); }
#line 1056 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 34:
#line 169 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnarySubstitutionTerm>("++*", yystack_[0].value.as< std::string > ()); }
#line 1062 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 35:
#line 170 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnarySubstitutionTerm>("--*", yystack_[0].value.as< std::string > ()); }
#line 1068 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 36:
#line 171 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnarySubstitutionTerm>("*++", yystack_[1].value.as< std::string > ()); }
#line 1074 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 37:
#line 172 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<UnarySubstitutionTerm>("*--", yystack_[1].value.as< std::string > ()); }
#line 1080 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 38:
#line 173 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<BinaryTerm>(yystack_[2].value.as< TermRef > (), "*", yystack_[0].value.as< TermRef > ()); }
#line 1086 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 39:
#line 174 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<BinaryTerm>(yystack_[2].value.as< TermRef > (), "/", yystack_[0].value.as< TermRef > ()); }
#line 1092 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 40:
#line 175 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< TermRef > () = createShared<BinaryTerm>(yystack_[2].value.as< TermRef > (), "%", yystack_[0].value.as< TermRef > ()); }
#line 1098 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 41:
#line 179 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< SingleTermRef > () = yystack_[0].value.as< ValueRef > (); }
#line 1104 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 42:
#line 180 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< SingleTermRef > () = createShared<FunctionCallSingleTerm>(yystack_[3].value.as< std::string > (), yystack_[1].value.as< ParameterVecRef > ()); }
#line 1110 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 43:
#line 181 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< SingleTermRef > () = createShared<EnclosedSingleTerm>(yystack_[1].value.as< ExpressionRef > ()); }
#line 1116 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 44:
#line 182 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< SingleTermRef > () = createShared<TupleElementSingleTerm>(yystack_[3].value.as< SingleTermRef > (), yystack_[1].value.as< ExpressionRef > ()); }
#line 1122 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 45:
#line 186 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = yystack_[0].value.as< ValueRef > (); }
#line 1128 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 46:
#line 187 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = yystack_[0].value.as< ValueRef > (); }
#line 1134 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 47:
#line 188 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = createShared<VariableValue>(yystack_[0].value.as< std::string > ()); }
#line 1140 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 48:
#line 192 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element::create<TypeMap::Bool>(yystack_[0].value.as< bool > ())); }
#line 1146 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 49:
#line 193 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element::create<TypeMap::Integer>(yystack_[0].value.as< int32_t > ())); }
#line 1152 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 50:
#line 194 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element::create<TypeMap::RealNumber>(yystack_[0].value.as< double > ())); }
#line 1158 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 51:
#line 195 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element::create<TypeMap::String>(yystack_[0].value.as< std::string > ())); }
#line 1164 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 52:
#line 199 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(createShared<std::vector<ParameterRef>>());
}
#line 1172 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 53:
#line 202 "SceneParser.yy" // lalr1.cc:859
    {
    ParameterVecRef elem = createShared<std::vector<ParameterRef>>();
    elem->push_back(yystack_[2].value.as< ParameterRef > ());
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(elem);
}
#line 1182 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 54:
#line 207 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(yystack_[1].value.as< ParameterVecRef > ());
}
#line 1190 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 55:
#line 213 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ArgumentDefinitionRef > () = createShared<ArgumentDefinition>(yystack_[0].value.as< std::string > ()); }
#line 1196 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 56:
#line 214 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ArgumentDefinitionRef > () = createShared<ArgumentDefinition>(yystack_[2].value.as< std::string > (), yystack_[0].value.as< ExpressionRef > ()); }
#line 1202 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 57:
#line 218 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentDefinitionVecRef > () = createShared<std::vector<ArgumentDefinitionRef>>();
}
#line 1210 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 58:
#line 221 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentDefinitionVecRef > () = createShared<std::vector<ArgumentDefinitionRef>>();
    yylhs.value.as< ArgumentDefinitionVecRef > ()->push_back(yystack_[0].value.as< ArgumentDefinitionRef > ());
}
#line 1219 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 59:
#line 225 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentDefinitionVecRef > () = yystack_[2].value.as< ArgumentDefinitionVecRef > ();
    yylhs.value.as< ArgumentDefinitionVecRef > ()->push_back(yystack_[0].value.as< ArgumentDefinitionRef > ());
}
#line 1228 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 60:
#line 232 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ParameterRef > () = createShared<Parameter>(nullptr, yystack_[0].value.as< ExpressionRef > ()); }
#line 1234 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 61:
#line 233 "SceneParser.yy" // lalr1.cc:859
    { yylhs.value.as< ParameterRef > () = createShared<Parameter>(yystack_[2].value.as< ExpressionRef > (), yystack_[0].value.as< ExpressionRef > ()); }
#line 1240 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 62:
#line 237 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterVecRef > () = createShared<std::vector<ParameterRef>>();
    yylhs.value.as< ParameterVecRef > ()->push_back(yystack_[2].value.as< ParameterRef > ());
    yylhs.value.as< ParameterVecRef > ()->push_back(yystack_[0].value.as< ParameterRef > ());
}
#line 1250 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 63:
#line 242 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterVecRef > () = yystack_[2].value.as< ParameterVecRef > ();
    yylhs.value.as< ParameterVecRef > ()->push_back(yystack_[0].value.as< ParameterRef > ());
}
#line 1259 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 64:
#line 249 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterVecRef > () = createShared<std::vector<ParameterRef>>();
}
#line 1267 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 65:
#line 252 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterVecRef > () = createShared<std::vector<ParameterRef>>();
    yylhs.value.as< ParameterVecRef > ()->push_back(yystack_[0].value.as< ParameterRef > ());
}
#line 1276 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 66:
#line 256 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ParameterVecRef > () = yystack_[2].value.as< ParameterVecRef > ();
    yylhs.value.as< ParameterVecRef > ()->push_back(yystack_[0].value.as< ParameterRef > ());
}
#line 1285 "SceneParser.tab.cc" // lalr1.cc:859
    break;


#line 1289 "SceneParser.tab.cc" // lalr1.cc:859
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


  const signed char SceneParser::yypact_ninf_ = -48;

  const signed char SceneParser::yytable_ninf_ = -1;

  const short int
  SceneParser::yypact_[] =
  {
     158,   -48,    50,   158,    39,    39,    39,   -37,   -15,   -48,
     -48,   -48,   -48,   293,    11,    19,   -14,   196,    98,   -48,
     317,    -7,    25,   -48,   -48,   -48,    43,   283,    22,    -1,
     142,    40,    25,    25,    25,   -48,   -48,   236,   236,   236,
     236,   236,   236,   236,   -48,   -48,   236,   236,    60,   -48,
     334,   -48,   -48,   236,   236,   247,   247,   236,   236,   236,
     236,   236,   236,   -48,   247,   247,   247,   236,   -48,   -48,
     236,   207,   -48,   236,   -48,   385,   -48,     0,   457,   457,
     457,   457,   457,   457,   406,   351,    28,   -48,     7,     7,
      20,    -7,    -7,     7,     7,   114,   114,   473,     9,   -48,
     -48,   -48,   441,   457,   -48,   -48,   -48,   -48,   236,   158,
     236,    41,   -48,     1,   -48,   -48,    29,   368,   236,   158,
      28,   158,   236,   457,   -48,   -48,   -48,   425,   158,   -48
  };

  const unsigned char
  SceneParser::yydefact_[] =
  {
       0,    12,     0,     0,     0,     0,     0,     0,     0,    48,
      49,    50,    51,    47,     0,     0,     0,     0,     0,     2,
       0,    13,    30,    41,    45,    46,     0,    60,     0,     0,
       0,    47,    31,    32,    33,    34,    35,    64,     0,     0,
       0,     0,     0,     0,    36,    37,     0,     0,     0,    10,
       0,     1,     3,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     4,     0,     0,     0,     0,    52,    43,
       0,     0,    54,     0,     5,    60,    65,     0,    24,    25,
      26,    27,    28,    29,     0,     0,    57,    11,    16,    17,
      47,    14,    15,    18,    19,    20,    21,    22,    23,    38,
      39,    40,     0,    61,    53,    62,    63,    42,     0,     0,
       0,    55,    58,     0,    44,    66,     6,     0,     0,     0,
       0,     0,     0,    56,     9,    59,     7,     0,     0,     8
  };

  const signed char
  SceneParser::yypgoto_[] =
  {
     -48,    69,   -17,    -2,   100,     5,   -48,   -48,   -48,   -47,
     -48,   -25,   -48,   -48
  };

  const signed char
  SceneParser::yydefgoto_[] =
  {
      -1,    18,    19,    20,    21,    22,    23,    24,    25,   112,
     113,    28,    29,    77
  };

  const unsigned char
  SceneParser::yytable_[] =
  {
      27,    52,    35,    72,   107,   119,    64,    65,    66,    32,
      33,    34,    76,    52,    46,    50,    53,    54,    55,    56,
      55,    56,    47,    37,    36,    48,    57,    58,    59,    60,
      61,    73,   108,   120,    67,    75,    78,    79,    80,    81,
      82,    83,     2,    37,    84,    85,   105,    68,   106,    44,
      45,    88,    89,     2,    71,    93,    94,    95,    96,    97,
      98,     4,     5,    86,   118,   102,     6,   111,   103,    75,
     121,    75,    30,   125,     9,    10,    11,    12,    31,     7,
       8,     0,    26,   115,     0,     9,    10,    11,    12,    13,
       0,     0,   116,     0,     0,     0,     0,     0,    51,     1,
       0,     2,   124,     3,   126,     0,    75,     0,   117,     4,
       5,   129,     0,     0,     6,     0,   123,     0,     0,     0,
     127,    53,    54,     0,     0,    55,    56,     7,     8,     0,
       0,    57,    58,     9,    10,    11,    12,    13,    14,     0,
      15,    16,    17,     1,     0,     2,     0,     3,    74,     0,
       0,     0,     0,     4,     5,    91,    92,     0,     6,     1,
       0,     2,     0,     3,    99,   100,   101,     0,     0,     4,
       5,     7,     8,     0,     6,     0,     0,     9,    10,    11,
      12,    13,    14,     0,    15,    16,    17,     7,     8,     0,
       0,     0,     0,     9,    10,    11,    12,    13,    14,     2,
      15,    16,    17,     0,     0,     0,     0,     4,     5,     0,
       2,   104,     6,     0,     0,     0,     0,     0,     4,     5,
       0,     0,     0,     6,     0,     7,     8,     0,     0,    49,
       0,     9,    10,    11,    12,    13,     7,     8,     0,     2,
       0,     0,     9,    10,    11,    12,    13,     4,     5,     0,
       2,     0,     6,     0,     0,     0,     0,     0,     4,     5,
       0,     0,     0,     6,     0,     7,     8,     0,     0,     0,
       0,     9,    10,    11,    12,    13,     7,     8,     0,     0,
       0,     0,     9,    10,    11,    12,    90,    69,     0,     0,
      53,    54,     0,     0,    55,    56,    37,     0,     0,     0,
      57,    58,    59,    60,    61,    62,     0,     0,     0,     0,
       0,     0,     0,     0,    70,     0,    38,    39,    40,    41,
      42,    43,    44,    45,    53,    54,     0,     0,    55,    56,
       0,     0,     0,     0,    57,    58,    59,    60,    61,    62,
       0,    53,    54,     0,     0,    55,    56,     0,     0,     0,
      63,    57,    58,    59,    60,    61,    62,     0,    53,    54,
       0,     0,    55,    56,     0,     0,     0,    87,    57,    58,
      59,    60,    61,    62,     0,    53,    54,     0,     0,    55,
      56,     0,     0,     0,   110,    57,    58,    59,    60,    61,
      62,     0,    53,    54,     0,     0,    55,    56,     0,     0,
       0,   122,    57,    58,    59,    60,    61,    62,     0,     0,
     109,     0,     0,    53,    54,     0,    70,    55,    56,     0,
       0,     0,     0,    57,    58,    59,    60,    61,    62,   128,
       0,     0,    53,    54,     0,     0,    55,    56,     0,     0,
       0,     0,    57,    58,    59,    60,    61,    62,    53,    54,
       0,   114,    55,    56,     0,     0,     0,     0,    57,    58,
      59,    60,    61,    62,    53,    54,     0,     0,    55,    56,
       0,     0,     0,     0,    57,    58,    59,    60,    61,    62,
      53,    54,     0,     0,    55,    56,     0,     0,     0,     0,
      57,    58,    59,    60
  };

  const short int
  SceneParser::yycheck_[] =
  {
       2,    18,    39,     4,     4,     4,    13,    14,    15,     4,
       5,     6,    37,    30,     3,    17,     7,     8,    11,    12,
      11,    12,     3,     3,    39,    39,    17,    18,    19,    20,
      21,    32,    32,    32,     9,    37,    38,    39,    40,    41,
      42,    43,     3,     3,    46,    47,    71,     4,    73,    29,
      30,    53,    54,     3,    32,    57,    58,    59,    60,    61,
      62,    11,    12,     3,    23,    67,    16,    39,    70,    71,
      41,    73,     3,   120,    35,    36,    37,    38,    39,    29,
      30,    -1,    32,   108,    -1,    35,    36,    37,    38,    39,
      -1,    -1,   109,    -1,    -1,    -1,    -1,    -1,     0,     1,
      -1,     3,   119,     5,   121,    -1,   108,    -1,   110,    11,
      12,   128,    -1,    -1,    16,    -1,   118,    -1,    -1,    -1,
     122,     7,     8,    -1,    -1,    11,    12,    29,    30,    -1,
      -1,    17,    18,    35,    36,    37,    38,    39,    40,    -1,
      42,    43,    44,     1,    -1,     3,    -1,     5,     6,    -1,
      -1,    -1,    -1,    11,    12,    55,    56,    -1,    16,     1,
      -1,     3,    -1,     5,    64,    65,    66,    -1,    -1,    11,
      12,    29,    30,    -1,    16,    -1,    -1,    35,    36,    37,
      38,    39,    40,    -1,    42,    43,    44,    29,    30,    -1,
      -1,    -1,    -1,    35,    36,    37,    38,    39,    40,     3,
      42,    43,    44,    -1,    -1,    -1,    -1,    11,    12,    -1,
       3,     4,    16,    -1,    -1,    -1,    -1,    -1,    11,    12,
      -1,    -1,    -1,    16,    -1,    29,    30,    -1,    -1,    33,
      -1,    35,    36,    37,    38,    39,    29,    30,    -1,     3,
      -1,    -1,    35,    36,    37,    38,    39,    11,    12,    -1,
       3,    -1,    16,    -1,    -1,    -1,    -1,    -1,    11,    12,
      -1,    -1,    -1,    16,    -1,    29,    30,    -1,    -1,    -1,
      -1,    35,    36,    37,    38,    39,    29,    30,    -1,    -1,
      -1,    -1,    35,    36,    37,    38,    39,     4,    -1,    -1,
       7,     8,    -1,    -1,    11,    12,     3,    -1,    -1,    -1,
      17,    18,    19,    20,    21,    22,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    31,    -1,    23,    24,    25,    26,
      27,    28,    29,    30,     7,     8,    -1,    -1,    11,    12,
      -1,    -1,    -1,    -1,    17,    18,    19,    20,    21,    22,
      -1,     7,     8,    -1,    -1,    11,    12,    -1,    -1,    -1,
      33,    17,    18,    19,    20,    21,    22,    -1,     7,     8,
      -1,    -1,    11,    12,    -1,    -1,    -1,    33,    17,    18,
      19,    20,    21,    22,    -1,     7,     8,    -1,    -1,    11,
      12,    -1,    -1,    -1,    33,    17,    18,    19,    20,    21,
      22,    -1,     7,     8,    -1,    -1,    11,    12,    -1,    -1,
      -1,    33,    17,    18,    19,    20,    21,    22,    -1,    -1,
       4,    -1,    -1,     7,     8,    -1,    31,    11,    12,    -1,
      -1,    -1,    -1,    17,    18,    19,    20,    21,    22,     4,
      -1,    -1,     7,     8,    -1,    -1,    11,    12,    -1,    -1,
      -1,    -1,    17,    18,    19,    20,    21,    22,     7,     8,
      -1,    10,    11,    12,    -1,    -1,    -1,    -1,    17,    18,
      19,    20,    21,    22,     7,     8,    -1,    -1,    11,    12,
      -1,    -1,    -1,    -1,    17,    18,    19,    20,    21,    22,
       7,     8,    -1,    -1,    11,    12,    -1,    -1,    -1,    -1,
      17,    18,    19,    20
  };

  const unsigned char
  SceneParser::yystos_[] =
  {
       0,     1,     3,     5,    11,    12,    16,    29,    30,    35,
      36,    37,    38,    39,    40,    42,    43,    44,    55,    56,
      57,    58,    59,    60,    61,    62,    32,    57,    65,    66,
      55,    39,    59,    59,    59,    39,    39,     3,    23,    24,
      25,    26,    27,    28,    29,    30,     3,     3,    39,    33,
      57,     0,    56,     7,     8,    11,    12,    17,    18,    19,
      20,    21,    22,    33,    13,    14,    15,     9,     4,     4,
      31,    32,     4,    32,     6,    57,    65,    67,    57,    57,
      57,    57,    57,    57,    57,    57,     3,    33,    57,    57,
      39,    58,    58,    57,    57,    57,    57,    57,    57,    58,
      58,    58,    57,    57,     4,    65,    65,     4,    32,     4,
      33,    39,    63,    64,    10,    65,    56,    57,    23,     4,
      32,    41,    33,    57,    56,    63,    56,    57,     4,    56
  };

  const unsigned char
  SceneParser::yyr1_[] =
  {
       0,    54,    55,    55,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    57,    57,    57,    57,    57,    57,    57,
      57,    57,    57,    57,    57,    57,    57,    57,    57,    57,
      58,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    59,    59,    59,    59,    60,    60,    60,    61,    61,
      61,    61,    62,    62,    62,    63,    63,    64,    64,    64,
      65,    65,    66,    66,    67,    67,    67
  };

  const unsigned char
  SceneParser::yyr2_[] =
  {
       0,     2,     1,     2,     2,     3,     5,     7,     9,     6,
       2,     3,     1,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     2,     2,     2,     2,     2,     2,     2,     3,     3,
       3,     1,     4,     3,     4,     1,     1,     1,     1,     1,
       1,     1,     3,     4,     3,     1,     3,     0,     1,     3,
       1,     3,     3,     3,     0,     1,     3
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const SceneParser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "\"(\"", "\")\"", "\"{\"",
  "\"}\"", "\"<\"", "\">\"", "\"[\"", "\"]\"", "\"+\"", "\"-\"", "\"*\"",
  "\"/\"", "\"%\"", "\"!\"", "\"<=\"", "\">=\"", "\"==\"", "\"!=\"",
  "\"&&\"", "\"||\"", "\"=\"", "\"+=\"", "\"-=\"", "\"*=\"", "\"/=\"",
  "\"%=\"", "\"++\"", "\"--\"", "\":\"", "\",\"", "\";\"", "CHAR", "BOOL",
  "INTEGER", "REALNUMBER", "STRING", "ID", "IF", "ELSE", "FOR", "FUNCTION",
  "RETURN", "PREC_SUBST", "PREC_LOGIC_OR", "PREC_LOGIC_AND", "PREC_EQ_REL",
  "PREC_INEQ_REL", "PREC_ADD", "PREC_MUL", "PREC_PRE_INC", "PREC_POST_INC",
  "$accept", "Statements", "Statement", "Expression", "Term", "SingleTerm",
  "Value", "ImmValue", "TupleValue", "ArgumentDefinition",
  "ArgumentDefinitions", "Parameter", "Elements", "Arguments", YY_NULLPTR
  };

#if YYDEBUG
  const unsigned short int
  SceneParser::yyrline_[] =
  {
       0,   115,   115,   121,   130,   131,   132,   133,   134,   135,
     136,   137,   138,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   179,   180,   181,   182,   186,   187,   188,   192,   193,
     194,   195,   199,   202,   207,   213,   214,   218,   221,   225,
     232,   233,   237,   242,   249,   252,   256
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
#line 1815 "SceneParser.tab.cc" // lalr1.cc:1167
#line 262 "SceneParser.yy" // lalr1.cc:1168


namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
