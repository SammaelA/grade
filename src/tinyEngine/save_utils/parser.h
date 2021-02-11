#pragma once
#include <stdio.h>
typedef struct _Param
{
char *name;
double param_base;
double param_from;
double param_to;
double distr_a;
double distr_b;
double *param_list;
int params;
int param_type;
int dist_type;
int regen_type;
} Param;
typedef struct _Preset
{
    char *name;
    Param params[50];
    int params_n;
}Preset;
typedef struct _ParseData
{
    Preset presets[100];
    int presets_n;
}ParseData;
extern char *cur_preset_name;
extern char *cur_parameter_name;
extern double param_base;
extern double param_from;
extern double param_to;
extern double distr_a;
extern double distr_b;
extern double *param_list;
extern int params;
extern int param_type;
extern int dist_type;
extern int regen_type;
extern ParseData pd;