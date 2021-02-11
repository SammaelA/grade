%{
#include <stdio.h>
#include "y.tab.h"
%}
%%

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