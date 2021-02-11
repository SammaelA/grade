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
    TOKPREEND = 258,
    TOKPREST = 259,
    TOKNORMAL = 260,
    TOKUNIFORM = 261,
    TOKNO_RANDOM = 262,
    TOKEXPLICIT_REGENERATION = 263,
    TOKREGENERATE_ON_STATE_CHANGE = 264,
    TOKREGENERATE_ON_GET = 265,
    TOKBS = 266,
    TOKBE = 267,
    TOKPS = 268,
    TOKPE = 269,
    TOKEQ = 270,
    TOKCOM = 271,
    TOKPARAMEND = 272,
    NUMBER = 273,
    NAME = 274
  };
#endif
/* Tokens.  */
#define TOKPREEND 258
#define TOKPREST 259
#define TOKNORMAL 260
#define TOKUNIFORM 261
#define TOKNO_RANDOM 262
#define TOKEXPLICIT_REGENERATION 263
#define TOKREGENERATE_ON_STATE_CHANGE 264
#define TOKREGENERATE_ON_GET 265
#define TOKBS 266
#define TOKBE 267
#define TOKPS 268
#define TOKPE 269
#define TOKEQ 270
#define TOKCOM 271
#define TOKPARAMEND 272
#define NUMBER 273
#define NAME 274

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 33 "presets.y" /* yacc.c:1909  */

        float number;
        char *string;

#line 97 "y.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
