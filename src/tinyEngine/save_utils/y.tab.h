/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    TST = 258,
    TEND = 259,
    TNAME = 260,
    TCOUNT = 261,
    TSYNTS = 262,
    TSYNTSPRES = 263,
    TPOS = 264,
    TOB = 265,
    TEL = 266,
    TCYL = 267,
    TSIZE = 268,
    TOKPREEND = 269,
    TOKPREST = 270,
    TOKNORMAL = 271,
    TOKUNIFORM = 272,
    TOKNO_RANDOM = 273,
    TOKEXPLICIT_REGENERATION = 274,
    TOKREGENERATE_ON_STATE_CHANGE = 275,
    TOKREGENERATE_ON_GET = 276,
    TOKBS = 277,
    TOKBE = 278,
    TOKPS = 279,
    TOKPE = 280,
    TOKEQ = 281,
    TOKCOM = 282,
    TOKPARAMEND = 283,
    NUMBER = 284,
    NAME = 285
  };
#endif
/* Tokens.  */
#define TST 258
#define TEND 259
#define TNAME 260
#define TCOUNT 261
#define TSYNTS 262
#define TSYNTSPRES 263
#define TPOS 264
#define TOB 265
#define TEL 266
#define TCYL 267
#define TSIZE 268
#define TOKPREEND 269
#define TOKPREST 270
#define TOKNORMAL 271
#define TOKUNIFORM 272
#define TOKNO_RANDOM 273
#define TOKEXPLICIT_REGENERATION 274
#define TOKREGENERATE_ON_STATE_CHANGE 275
#define TOKREGENERATE_ON_GET 276
#define TOKBS 277
#define TOKBE 278
#define TOKPS 279
#define TOKPE 280
#define TOKEQ 281
#define TOKCOM 282
#define TOKPARAMEND 283
#define NUMBER 284
#define NAME 285

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 34 "presets.y" /* yacc.c:1909  */

        float number;
        char *string;

#line 119 "y.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
