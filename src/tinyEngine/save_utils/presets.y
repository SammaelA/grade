%require "3.0"

%{
#include <stdio.h>
#include <string.h>
#include "parser.h"
extern void start_preset(void);
extern void end_preset(void);
extern void start_param(void);
extern void end_param(void);

void yyerror(const char *str)
{
	fprintf(stderr,"error: %s\n",str);
}

int yywrap()
{
	return 1;
}

int parse_file(char *file)
{
    freopen(file,"r",stdin);
	yyparse();
    freopen("/dev/stdin", "r", stdin);
}

%}

%token TOKPREEND TOKPREST TOKNORMAL TOKUNIFORM TOKNO_RANDOM TOKEXPLICIT_REGENERATION TOKREGENERATE_ON_STATE_CHANGE TOKREGENERATE_ON_GET TOKBS TOKBE TOKPS TOKPE TOKEQ TOKCOM TOKPARAMEND
%union 
{
        float number;
        char *string;
}

%token <number> NUMBER
%token <string> NAME
%%

commands:
	| commands command

	;


command:
	preset_start preset_name preset_body preset_end
    {

    }
	;

preset_start:
	TOKPREST 
	{
        start_preset();
	}
	;

preset_name:
	NAME
	{
        cur_preset_name = strdup($1);
	}
	;
preset_body:
	TOKBS presets_parameters TOKBE
	;
preset_end:
	TOKPREEND
	{
        end_preset();
	}
	;
presets_parameters:
    | presets_parameters presets_param
	;
presets_param:
    param_name TOKEQ param1 param_end |
    param_name TOKEQ param2 param_end |
    param_name TOKEQ param3 param_end |
    param_name TOKEQ param4 param_end
    {

    }
    ;
param_name:
    NAME
    {
        cur_parameter_name = $1;
    }
    ;
param_end:
    TOKPARAMEND
    {
        end_param();
    }
    ;
param1:
    NUMBER
    {
        start_param();
        param_base = $1;
        param_type = 1;
    };
param2:
    param1 TOKCOM TOKBS list TOKBE
    {
        param_type = 2;
    };
param3:
    param2 TOKCOM distribution TOKPS NUMBER TOKCOM NUMBER TOKCOM regenerate_type TOKPE 
    {
        distr_a = $5;
        distr_b = $7;
        param_type = 3;
    };
param4:
    param3 TOKCOM TOKPS NUMBER TOKCOM NUMBER TOKPE 
    {
        param_from = $4;
        param_to = $6;
        param_type = 4;
    };
distribution:
    norm | unif
    ;
norm:
    TOKNORMAL
    {
        dist_type = 1;
    };
unif:
    TOKUNIFORM
    {
        dist_type = 2;
    };
regenerate_type:
    r1 | r2 | r3 | r4
    ;
r1: TOKNO_RANDOM
    {
        regen_type = 1;
    };
r2: TOKEXPLICIT_REGENERATION
    {
        regen_type = 2;
    };
r3: TOKREGENERATE_ON_STATE_CHANGE
    { 
        regen_type = 3;
    };
r4: TOKREGENERATE_ON_GET
    {
        regen_type = 4;
    };
list:
     | list list_elem
    {

    };
list_elem:
    NUMBER
    {
        param_list[params] = $1;
        params++;
    };