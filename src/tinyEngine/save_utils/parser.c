#include "parser.h"
#include <malloc.h>
#include <string.h>
extern int parse_file(void);
ParseData pd;
char *cur_preset_name = NULL;
 char *cur_parameter_name = NULL;
 double param_base;
 double param_from;
 double param_to;
 double distr_a;
 double distr_b;
 double *param_list = NULL;
 int params = 0;
 int param_type;
 int dist_type;
 int regen_type;
 void start_preset()
 {
    cur_preset_name = (char *)calloc(100,sizeof(char));
    pd.presets[pd.presets_n].params_n = 0;
 }
 void end_preset()
 {
    pd.presets[pd.presets_n].name = cur_preset_name;
    pd.presets_n++;
 }
 void start_param()
 {
    param_list = (double *)calloc(10,sizeof(double));
    params = 0;
 }
 void end_param()
 {
     int pn = pd.presets[pd.presets_n].params_n;
     pd.presets[pd.presets_n].params[pn].dist_type = dist_type;
     pd.presets[pd.presets_n].params[pn].distr_a = distr_a;
     pd.presets[pd.presets_n].params[pn].distr_b = distr_b;
     pd.presets[pd.presets_n].params[pn].name = strdup(cur_parameter_name);
     pd.presets[pd.presets_n].params[pn].param_base = param_base;
     pd.presets[pd.presets_n].params[pn].param_from = param_from;
     pd.presets[pd.presets_n].params[pn].param_to = param_to;
     pd.presets[pd.presets_n].params[pn].param_type = param_type;
     pd.presets[pd.presets_n].params[pn].params = params;
     pd.presets[pd.presets_n].params[pn].regen_type = regen_type;
     pd.presets[pd.presets_n].params[pn].param_list = param_list;
     pd.presets[pd.presets_n].params_n++;
 }