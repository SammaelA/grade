#pragma once

typedef struct _obst
{
    int type;
    float pos[3];
    float a[3];
    float b[3];
    float c[3];
}obst;
typedef struct _ttd
{
    int id;
    char *name;
    char *wood;
    char *leaf;
} ttd;
typedef struct _ggd
{
    char *name;
    int count;
    int synts_count;
    int synts_precision;
    float pos[3];
    float size[3];
    ttd ttds[10];
    obst obsts[100];
    int ttds_c;
    int obsts_c;
} ggd;
typedef struct _ggddat
{
    ggd ggds[50];
    int ggds_c;
} ggddat;

extern ggddat dat;
#define CG dat.ggds[dat.ggds_c]
#define CT CG.ttds[CG.ttds_c]
#define CO CG.obsts[CG.obsts_c]