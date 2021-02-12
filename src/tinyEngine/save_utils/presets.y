%require "3.0"

%{
#include <stdio.h>
#include <string.h>
#include "parser.h"
#include "ggd_parser.h"
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

%token TST TEND TNAME TCOUNT TPOS TOB TEL TCYL TSIZE TOKPREEND TOKPREST TOKNORMAL TOKUNIFORM TOKNO_RANDOM TOKEXPLICIT_REGENERATION TOKREGENERATE_ON_STATE_CHANGE TOKREGENERATE_ON_GET TOKBS TOKBE TOKPS TOKPE TOKEQ TOKCOM TOKPARAMEND
%union 
{
        float number;
        char *string;
}

%token <number> NUMBER
%token <string> NAME
%%
commands:
    prsts | groves
    ;
prsts:
	| prsts prst

	;


prst:
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
    groves:
	grove | groves grove

	;


grove: ggd_head TOKBS count pos size presets obstacles TOKBE TEND
{
    dat.ggds_c++;
};
ggd_head: TST NAME

{
    CG.ttds_c = 0;
    CG.obsts_c = 0;
    CG.name = strdup($2);
};
count: TCOUNT TOKEQ NUMBER
{
    CG.count = $3;
};
pos: TPOS TOKEQ TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CG.pos[0] = $4;
    CG.pos[1] = $5;
    CG.pos[2] = $6;

};
size: TSIZE TOKEQ TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CG.size[0] = $4;
    CG.size[1] = $5;
    CG.size[2] = $6;
};
presets: preset | presets preset
{

};
preset: TOKPREST TOKEQ TOKBS NAME NAME NAME TOKBE
{
    CT.id = CG.ttds_c;
    CT.name = strdup($4);
    CT.wood = strdup($5);
    CT.leaf = strdup($6);
    CG.ttds_c++;
};
obstacles: | obstacles obstacle
{

};
obstacle: ob_type TOKEQ TOKBS opos oa ob oc TOKBE
{
    CG.obsts_c++;
};
ob_type: ob_b | ob_e | ob_c
{

};
ob_b: TOB
{
    CO.type = 1;
};
ob_e: TEL
{
    CO.type = 2;
};
ob_c: TCYL
{
    CO.type = 3;
};
opos: TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CO.pos[0] = $2;
    CO.pos[1] = $3;
    CO.pos[2] = $4;
};
oa: TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CO.a[0] = $2;
    CO.a[1] = $3;
    CO.a[2] = $4;
};
ob: TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CO.b[0] = $2;
    CO.b[1] = $3;
    CO.b[2] = $4;
};
oc: TOKPS NUMBER  NUMBER  NUMBER TOKPE
{
    CO.c[0] = $2;
    CO.c[1] = $3;
    CO.c[2] = $4;
};