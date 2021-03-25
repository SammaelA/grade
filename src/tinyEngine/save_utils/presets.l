%{
#include <stdio.h>
#include "y.tab.h"
%}
%%
grove_end              return TEND;
grove                  return TST;
name                  return TNAME;
trees_count return TCOUNT;
pos return TPOS;
count return TCOUNT;
synthetic_precision return TSYNTSPRES;
synthetic_count return TSYNTS;
box return TOB;
ellipsoid return TEL;
cylinder return TCYL;
Box return TOB;
Ellipsoid return TEL;
Cylinder return TCYL;
size return TSIZE;
preset_end              return TOKPREEND;
preset                  return TOKPREST;
normal                  return TOKNORMAL;
uniform                 return TOKUNIFORM;
NORMAL                  return TOKNORMAL;
UNIFORM                 return TOKUNIFORM;
NO_RANDOM return TOKNO_RANDOM;
EXPLICIT_REGENERATION return TOKEXPLICIT_REGENERATION;
REGENERATE_ON_STATE_CHANGE return TOKREGENERATE_ON_STATE_CHANGE;
REGENERATE_ON_GET return TOKREGENERATE_ON_GET;
\{                      return TOKBS;
\}                      return TOKBE;
\(                      return TOKPS;
\)                      return TOKPE;
\=                      return TOKEQ;
\,                      return TOKCOM;
\;                      return TOKPARAMEND;
[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)    yylval.number=atof(yytext); return NUMBER;
[a-zA-Z_-][a-zA-Z0-9_-]*    yylval.string=strdup(yytext); return NAME;
\n                      /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;
%%