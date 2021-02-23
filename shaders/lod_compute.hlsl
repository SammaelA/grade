#version 430
#define THREADS 1
layout( local_size_x = THREADS ) in;

struct SliceVertexData 
{
    vec4 position;
	vec4 tcs;
};
layout(std140, binding=2) buffer Slices 
{
    SliceVertexData sliceVertexes[];
};
struct LOD_info
{
    uint offset;
    uint pad;
    vec2 min_max;
};
struct InstanceData
{
    vec4 center_self;
    vec4 center_par;
    mat4 projection_camera;
};
struct IndexData
{
    uint index;
    //uint pad;
    float mn;
    float mx;
};
layout(std140, binding=0) buffer _LODs
{
    readonly LOD_info LODs[];
};
layout(std140, binding=1) buffer _models_intervals
{
    readonly uvec4 models_intervals[];
};
layout(std140, binding=3) buffer _instances
{
    readonly InstanceData instances[];
};
layout(std140, binding=4) buffer _instance_indexes//same size as instances
{
    writeonly IndexData instance_indexes[];
};
layout(std140, binding=5) buffer _models_instance_count//same size as models_intervals
{
    volatile uint models_instance_count[];
};

uniform uint lods_count;
uniform vec3 camera_pos;
uniform float trans;
void main()
{
    uvec2 interval = models_intervals[gl_WorkGroupID.x].xy;
    int k = -1;
    for (int i=0;i<lods_count - 1;i++)
    {
        if (LODs[i].offset<=interval.x && LODs[i + 1].offset>=interval.y)
        {
            k = i;
            break;
        }
    }
    if (interval.x >= LODs[lods_count - 1].offset)
        k = int(lods_count) - 1;
    if (k == -1)
    {
        return;
    }
    vec2 min_max = LODs[k].min_max;

    uint cnt = (interval.y - interval.x)/ THREADS + 1;
    uint st, en;
    st = interval.x + cnt * gl_LocalInvocationIndex;
    en = interval.x + cnt * (gl_LocalInvocationIndex + 1);
    en = en > interval.y ? interval.y : en;
    st = interval.x;
    en = interval.y;
    uint inst_num = 0;
    for (uint i = st; i < en; i++)
    {
        float mx_dist = min_max.y - length(instances[i].center_par.xyz - camera_pos);
	    float mn_dist = length(instances[i].center_self.xyz - camera_pos) - min_max.x;

        if ((mx_dist >= -trans) && (mn_dist >= -trans))//test if we need this to draw this instance
        {
            vec2 a_mult = vec2(0,0);
            if (mn_dist < trans)
            {
                a_mult = vec2(0.5*(mn_dist/trans + 1), length(instances[i].center_par.xyz - camera_pos));
            }
            else if (mx_dist < trans)
            {
                a_mult = vec2(0.5*(mx_dist/trans + 1), -length(instances[i].center_self.xyz - camera_pos));
            }
            else
            {
                a_mult = vec2(10,0);
            }
            instance_indexes[st + inst_num].index = i;
            instance_indexes[st + inst_num].mn = a_mult.x;
            instance_indexes[st + inst_num].mx = a_mult.y;
            inst_num++;
        }
        else
        {
            //instance_indexes[i].index = 0;
            //instance_indexes[i].pad = 0;
            //instance_indexes[i].mn = 0.0;
            //instance_indexes[i].mx = 0.0;
        }
    }
    /*if (gl_LocalInvocationIndex == 0)
    {
        //models_instance_count[1] = 3;
        models_instance_count[gl_WorkGroupID.x] = inst_num;
    }
    else*/
        atomicAdd(models_instance_count[gl_WorkGroupID.x], inst_num);
}
/*
uniform uint sv_size;

void main()
{
    uint st,en;
    uint cnt = sv_size / THREADS + 1;

    st = cnt * gl_LocalInvocationIndex;
    en = cnt * (gl_LocalInvocationIndex + 1);
    en = en > sv_size ? sv_size : en;
    //sliceVertexes[gl_LocalInvocationIndex].position.y = gl_LocalInvocationIndex;
    for (uint i = st;i<en;i++)
        sliceVertexes[i].position.y += 0.02;
}*/