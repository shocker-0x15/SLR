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
      case 11: // API
        value.move< API > (that.value);
        break;

      case 34: // Argument
        value.move< ArgumentRef > (that.value);
        break;

      case 31: // Elements
      case 33: // Arguments
        value.move< ArgumentsRef > (that.value);
        break;

      case 26: // Expression
        value.move< ExpressionRef > (that.value);
        break;

      case 25: // Statement
        value.move< StatementRef > (that.value);
        break;

      case 24: // Statements
        value.move< StatementsRef > (that.value);
        break;

      case 27: // Term
      case 32: // Function
        value.move< TermRef > (that.value);
        break;

      case 28: // Value
      case 29: // ImmValue
      case 30: // TupleValue
        value.move< ValueRef > (that.value);
        break;

      case 12: // BOOL
        value.move< bool > (that.value);
        break;

      case 7: // ","
      case 8: // ":"
      case 9: // ";"
      case 10: // CHAR
        value.move< char > (that.value);
        break;

      case 14: // REALNUMBER
        value.move< double > (that.value);
        break;

      case 13: // INTEGER
        value.move< int32_t > (that.value);
        break;

      case 15: // STRING
      case 16: // ID
      case 18: // SUBSTITUTION
      case 19: // OPERATOR0
      case 20: // OPERATOR1
      case 22: // INC_DEC
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
      case 11: // API
        value.copy< API > (that.value);
        break;

      case 34: // Argument
        value.copy< ArgumentRef > (that.value);
        break;

      case 31: // Elements
      case 33: // Arguments
        value.copy< ArgumentsRef > (that.value);
        break;

      case 26: // Expression
        value.copy< ExpressionRef > (that.value);
        break;

      case 25: // Statement
        value.copy< StatementRef > (that.value);
        break;

      case 24: // Statements
        value.copy< StatementsRef > (that.value);
        break;

      case 27: // Term
      case 32: // Function
        value.copy< TermRef > (that.value);
        break;

      case 28: // Value
      case 29: // ImmValue
      case 30: // TupleValue
        value.copy< ValueRef > (that.value);
        break;

      case 12: // BOOL
        value.copy< bool > (that.value);
        break;

      case 7: // ","
      case 8: // ":"
      case 9: // ";"
      case 10: // CHAR
        value.copy< char > (that.value);
        break;

      case 14: // REALNUMBER
        value.copy< double > (that.value);
        break;

      case 13: // INTEGER
        value.copy< int32_t > (that.value);
        break;

      case 15: // STRING
      case 16: // ID
      case 18: // SUBSTITUTION
      case 19: // OPERATOR0
      case 20: // OPERATOR1
      case 22: // INC_DEC
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
            case 7: // ","

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 438 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 8: // ":"

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 445 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 9: // ";"

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 452 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 10: // CHAR

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 459 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 11: // API

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 466 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 12: // BOOL

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 473 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 13: // INTEGER

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 480 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 14: // REALNUMBER

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 487 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 15: // STRING

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 494 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 16: // ID

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 501 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 18: // SUBSTITUTION

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 508 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 19: // OPERATOR0

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 515 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 20: // OPERATOR1

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 522 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 22: // INC_DEC

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 529 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 24: // Statements

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 536 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 25: // Statement

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 543 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 26: // Expression

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 550 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 27: // Term

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 557 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 28: // Value

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 564 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 29: // ImmValue

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 571 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 30: // TupleValue

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 578 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 31: // Elements

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 585 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 32: // Function

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 592 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 33: // Arguments

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 599 "SceneParser.tab.cc" // lalr1.cc:636
        break;

      case 34: // Argument

#line 70 "SceneParser.yy" // lalr1.cc:636
        { /*yyoutput << $$;*/ }
#line 606 "SceneParser.tab.cc" // lalr1.cc:636
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

#line 725 "SceneParser.tab.cc" // lalr1.cc:745

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
      case 11: // API
        yylhs.value.build< API > ();
        break;

      case 34: // Argument
        yylhs.value.build< ArgumentRef > ();
        break;

      case 31: // Elements
      case 33: // Arguments
        yylhs.value.build< ArgumentsRef > ();
        break;

      case 26: // Expression
        yylhs.value.build< ExpressionRef > ();
        break;

      case 25: // Statement
        yylhs.value.build< StatementRef > ();
        break;

      case 24: // Statements
        yylhs.value.build< StatementsRef > ();
        break;

      case 27: // Term
      case 32: // Function
        yylhs.value.build< TermRef > ();
        break;

      case 28: // Value
      case 29: // ImmValue
      case 30: // TupleValue
        yylhs.value.build< ValueRef > ();
        break;

      case 12: // BOOL
        yylhs.value.build< bool > ();
        break;

      case 7: // ","
      case 8: // ":"
      case 9: // ";"
      case 10: // CHAR
        yylhs.value.build< char > ();
        break;

      case 14: // REALNUMBER
        yylhs.value.build< double > ();
        break;

      case 13: // INTEGER
        yylhs.value.build< int32_t > ();
        break;

      case 15: // STRING
      case 16: // ID
      case 18: // SUBSTITUTION
      case 19: // OPERATOR0
      case 20: // OPERATOR1
      case 22: // INC_DEC
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
#line 84 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementsRef > () = createShared<std::vector<StatementRef>>();
    if (yystack_[0].value.as< StatementRef > ())
        yylhs.value.as< StatementsRef > ()->push_back(yystack_[0].value.as< StatementRef > ());
    driver.statements = yylhs.value.as< StatementsRef > ();
}
#line 904 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 3:
#line 90 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementsRef > () = yystack_[1].value.as< StatementsRef > ();
    if (yystack_[0].value.as< StatementRef > ())
        yylhs.value.as< StatementsRef > ()->push_back(yystack_[0].value.as< StatementRef > ()); 
    driver.statements = yylhs.value.as< StatementsRef > ();
}
#line 915 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 4:
#line 99 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementRef > () = yystack_[1].value.as< ExpressionRef > ();
}
#line 923 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 5:
#line 102 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< StatementRef > () = createShared<ForStatement>(yystack_[8].value.as< ExpressionRef > (), yystack_[6].value.as< ExpressionRef > (), yystack_[4].value.as< ExpressionRef > (), yystack_[1].value.as< StatementsRef > ());
}
#line 931 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 6:
#line 105 "SceneParser.yy" // lalr1.cc:859
    {
    printf("Parsing aborted.\n");
    YYABORT;
}
#line 940 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 7:
#line 112 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ExpressionRef > () = yystack_[0].value.as< TermRef > ();
}
#line 948 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 8:
#line 115 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ExpressionRef > () = createShared<BinaryExpression>(yystack_[2].value.as< TermRef > (), yystack_[1].value.as< std::string > (), yystack_[0].value.as< TermRef > ());
}
#line 956 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 9:
#line 118 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ExpressionRef > () = createShared<SubstitutionExpression>(yystack_[2].value.as< std::string > (), yystack_[1].value.as< std::string > (), yystack_[0].value.as< ExpressionRef > ());
}
#line 964 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 10:
#line 124 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = yystack_[0].value.as< ValueRef > ();
}
#line 972 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 11:
#line 127 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = yystack_[0].value.as< TermRef > ();
}
#line 980 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 12:
#line 130 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = createShared<UnaryTerm>(yystack_[1].value.as< std::string > (), yystack_[0].value.as< TermRef > ());
}
#line 988 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 13:
#line 133 "SceneParser.yy" // lalr1.cc:859
    {
    
}
#line 996 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 14:
#line 136 "SceneParser.yy" // lalr1.cc:859
    {

}
#line 1004 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 15:
#line 139 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = createShared<BinaryTerm>(yystack_[2].value.as< TermRef > (), yystack_[1].value.as< std::string > (), yystack_[0].value.as< TermRef > ());
}
#line 1012 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 16:
#line 142 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = createShared<EnclosedTerm>(yystack_[1].value.as< ExpressionRef > ());
}
#line 1020 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 17:
#line 148 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = yystack_[0].value.as< ValueRef > ();
}
#line 1028 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 18:
#line 151 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = yystack_[0].value.as< ValueRef > ();
}
#line 1036 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 19:
#line 154 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<VariableValue>(yystack_[0].value.as< std::string > ());
}
#line 1044 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 20:
#line 160 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element(yystack_[0].value.as< bool > ()));
}
#line 1052 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 21:
#line 163 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element(yystack_[0].value.as< int32_t > ()));
}
#line 1060 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 22:
#line 166 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element(yystack_[0].value.as< double > ()));
}
#line 1068 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 23:
#line 169 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<ImmediateValue>(Element(yystack_[0].value.as< std::string > ()));
}
#line 1076 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 24:
#line 175 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(createShared<std::vector<ArgumentRef>>());
}
#line 1084 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 25:
#line 178 "SceneParser.yy" // lalr1.cc:859
    {
    ArgumentsRef elem = createShared<std::vector<ArgumentRef>>();
    elem->push_back(yystack_[2].value.as< ArgumentRef > ());
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(elem);
}
#line 1094 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 26:
#line 183 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ValueRef > () = createShared<TupleValue>(yystack_[1].value.as< ArgumentsRef > ());
}
#line 1102 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 27:
#line 189 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentsRef > () = createShared<std::vector<ArgumentRef>>();
    yylhs.value.as< ArgumentsRef > ()->push_back(yystack_[2].value.as< ArgumentRef > ());
    yylhs.value.as< ArgumentsRef > ()->push_back(yystack_[0].value.as< ArgumentRef > ());
}
#line 1112 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 28:
#line 194 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentsRef > () = yystack_[2].value.as< ArgumentsRef > ();
    yylhs.value.as< ArgumentsRef > ()->push_back(yystack_[0].value.as< ArgumentRef > ());
}
#line 1121 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 29:
#line 201 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< TermRef > () = createShared<FunctionTerm>(yystack_[3].value.as< API > (), yystack_[1].value.as< ArgumentsRef > ());
}
#line 1129 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 30:
#line 207 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentsRef > () = createShared<std::vector<ArgumentRef>>();
}
#line 1137 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 31:
#line 210 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentsRef > () = createShared<std::vector<ArgumentRef>>();
    yylhs.value.as< ArgumentsRef > ()->push_back(yystack_[0].value.as< ArgumentRef > ());
}
#line 1146 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 32:
#line 214 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentsRef > () = yystack_[2].value.as< ArgumentsRef > ();
    yylhs.value.as< ArgumentsRef > ()->push_back(yystack_[0].value.as< ArgumentRef > ());
}
#line 1155 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 33:
#line 221 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentRef > () = createShared<Argument>(nullptr, yystack_[0].value.as< ExpressionRef > ());
}
#line 1163 "SceneParser.tab.cc" // lalr1.cc:859
    break;

  case 34:
#line 224 "SceneParser.yy" // lalr1.cc:859
    {
    yylhs.value.as< ArgumentRef > () = createShared<Argument>(yystack_[2].value.as< ExpressionRef > (), yystack_[0].value.as< ExpressionRef > ());
}
#line 1171 "SceneParser.tab.cc" // lalr1.cc:859
    break;


#line 1175 "SceneParser.tab.cc" // lalr1.cc:859
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


  const signed char SceneParser::yypact_ninf_ = -28;

  const signed char SceneParser::yytable_ninf_ = -1;

  const short int
  SceneParser::yypact_[] =
  {
      80,   -28,    97,     9,   -28,   -28,   -28,   -28,   -14,    10,
     125,   -10,    46,   -28,     8,    -4,   -28,   -28,   -28,   -28,
      16,     1,     3,    14,   139,   139,   -28,   139,     2,   -28,
     -28,   -28,   -28,   -28,   125,   125,   -28,   -28,   139,   -28,
     139,   111,    18,     7,   -28,   -28,    20,    11,   -28,   -28,
     -28,   -28,   -28,   -28,   139,   139,   -28,    21,   139,    29,
      30,    80,    63,   -28
  };

  const unsigned char
  SceneParser::yydefact_[] =
  {
       0,     6,     0,     0,    20,    21,    22,    23,    19,     0,
       0,     0,     0,     2,     0,     7,    10,    17,    18,    11,
       0,    33,     0,     0,    30,     0,    14,     0,    19,    12,
      13,     1,     3,     4,     0,     0,    24,    16,     0,    26,
       0,     0,    33,     0,    31,     9,     0,     8,    15,    34,
      28,    25,    27,    29,     0,     0,    32,     0,     0,     0,
       0,     0,     0,     5
  };

  const signed char
  SceneParser::yypgoto_[] =
  {
     -28,   -27,   -11,    -2,    -7,   -28,   -28,   -28,   -28,   -28,
     -28,   -22
  };

  const signed char
  SceneParser::yydefgoto_[] =
  {
      -1,    12,    13,    14,    15,    16,    17,    18,    22,    19,
      43,    23
  };

  const unsigned char
  SceneParser::yytable_[] =
  {
      21,    32,    44,    29,    25,    37,    30,    39,    26,    38,
      40,    53,    24,    27,    54,    34,    35,    33,    50,    52,
      36,    41,    42,    45,    26,    46,    38,    47,    48,    55,
      58,    35,    56,    60,    62,    61,    49,     0,    42,    42,
       0,     0,     0,     0,     0,     0,    31,     1,     0,     2,
       0,    32,    42,    57,     0,     0,    59,     3,     4,     5,
       6,     7,     8,     9,     1,    10,     2,     0,    11,    63,
       0,     0,     0,     0,     3,     4,     5,     6,     7,     8,
       9,     1,    10,     2,     0,    11,     0,     0,     0,     0,
       0,     3,     4,     5,     6,     7,     8,     9,     0,    10,
       2,     0,    11,     0,    20,     0,     0,     0,     3,     4,
       5,     6,     7,     8,     2,    51,    10,     0,     0,    11,
       0,     0,     3,     4,     5,     6,     7,     8,     2,     0,
      10,     0,     0,    11,     0,     0,     3,     4,     5,     6,
       7,    28,     2,     0,    10,     0,     0,    11,     0,     0,
       3,     4,     5,     6,     7,     8,     0,     0,    10,     0,
       0,    11
  };

  const signed char
  SceneParser::yycheck_[] =
  {
       2,    12,    24,    10,    18,     4,    16,     4,    22,     8,
       7,     4,     3,     3,     7,    19,    20,     9,    40,    41,
       4,     7,    24,    25,    22,    27,     8,    34,    35,     9,
       9,    20,    54,     4,    61,     5,    38,    -1,    40,    41,
      -1,    -1,    -1,    -1,    -1,    -1,     0,     1,    -1,     3,
      -1,    62,    54,    55,    -1,    -1,    58,    11,    12,    13,
      14,    15,    16,    17,     1,    19,     3,    -1,    22,     6,
      -1,    -1,    -1,    -1,    11,    12,    13,    14,    15,    16,
      17,     1,    19,     3,    -1,    22,    -1,    -1,    -1,    -1,
      -1,    11,    12,    13,    14,    15,    16,    17,    -1,    19,
       3,    -1,    22,    -1,     7,    -1,    -1,    -1,    11,    12,
      13,    14,    15,    16,     3,     4,    19,    -1,    -1,    22,
      -1,    -1,    11,    12,    13,    14,    15,    16,     3,    -1,
      19,    -1,    -1,    22,    -1,    -1,    11,    12,    13,    14,
      15,    16,     3,    -1,    19,    -1,    -1,    22,    -1,    -1,
      11,    12,    13,    14,    15,    16,    -1,    -1,    19,    -1,
      -1,    22
  };

  const unsigned char
  SceneParser::yystos_[] =
  {
       0,     1,     3,    11,    12,    13,    14,    15,    16,    17,
      19,    22,    24,    25,    26,    27,    28,    29,    30,    32,
       7,    26,    31,    34,     3,    18,    22,     3,    16,    27,
      16,     0,    25,     9,    19,    20,     4,     4,     8,     4,
       7,     7,    26,    33,    34,    26,    26,    27,    27,    26,
      34,     4,    34,     4,     7,     9,    34,    26,     9,    26,
       4,     5,    24,     6
  };

  const unsigned char
  SceneParser::yyr1_[] =
  {
       0,    23,    24,    24,    25,    25,    25,    26,    26,    26,
      27,    27,    27,    27,    27,    27,    27,    28,    28,    28,
      29,    29,    29,    29,    30,    30,    30,    31,    31,    32,
      33,    33,    33,    34,    34
  };

  const unsigned char
  SceneParser::yyr2_[] =
  {
       0,     2,     1,     2,     2,    11,     1,     1,     3,     3,
       1,     1,     2,     2,     2,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     3,     4,     3,     3,     3,     4,
       0,     1,     3,     1,     3
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const SceneParser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "\"(\"", "\")\"", "\"{\"",
  "\"}\"", "\",\"", "\":\"", "\";\"", "CHAR", "API", "BOOL", "INTEGER",
  "REALNUMBER", "STRING", "ID", "FOR", "SUBSTITUTION", "OPERATOR0",
  "OPERATOR1", "NEG", "INC_DEC", "$accept", "Statements", "Statement",
  "Expression", "Term", "Value", "ImmValue", "TupleValue", "Elements",
  "Function", "Arguments", "Argument", YY_NULLPTR
  };

#if YYDEBUG
  const unsigned char
  SceneParser::yyrline_[] =
  {
       0,    84,    84,    90,    99,   102,   105,   112,   115,   118,
     124,   127,   130,   133,   136,   139,   142,   148,   151,   154,
     160,   163,   166,   169,   175,   178,   183,   189,   194,   201,
     207,   210,   214,   221,   224
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
#line 1603 "SceneParser.tab.cc" // lalr1.cc:1167
#line 229 "SceneParser.yy" // lalr1.cc:1168


namespace SLRSceneGraph {
    void SceneParser::error(const location_type& l, const std::string& m) {
        driver.error(l, m);
    }
}
