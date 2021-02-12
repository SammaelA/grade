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
    TPOS = 262,
    TOB = 263,
    TEL = 264,
    TCYL = 265,
    TSIZE = 266,
    TOKPREEND = 267,
    TOKPREST = 268,
    TOKNORMAL = 269,
    TOKUNIFORM = 270,
    TOKNO_RANDOM = 271,
    TOKEXPLICIT_REGENERATION = 272,
    TOKREGENERATE_ON_STATE_CHANGE = 273,
    TOKREGENERATE_ON_GET = 274,
    TOKBS = 275,
    TOKBE = 276,
    TOKPS = 277,
    TOKPE = 278,
    TOKEQ = 279,
    TOKCOM = 280,
    TOKPARAMEND = 281,
    NUMBER = 282,
    NAME = 283
  };
#endif
/* Tokens.  */
#define TST 258
#define TEND 259
#define TNAME 260
#define TCOUNT 261
#define TPOS 262
#define TOB 263
#define TEL 264
#define TCYL 265
#define TSIZE 266
#define TOKPREEND 267
#define TOKPREST 268
#define TOKNORMAL 269
#define TOKUNIFORM 270
#define TOKNO_RANDOM 271
#define TOKEXPLICIT_REGENERATION 272
#define TOKREGENERATE_ON_STATE_CHANGE 273
#define TOKREGENERATE_ON_GET 274
#define TOKBS 275
#define TOKBE 276
#define TOKPS 277
#define TOKPE 278
#define TOKEQ 279
#define TOKCOM 280
#define TOKPARAMEND 281
#define NUMBER 282
#define NAME 283

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 34 "presets.y" /* yacc.c:1909  */

        float number;
        char *string;

#line 115 "y.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
