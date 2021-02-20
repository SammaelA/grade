#version 430
#define THREADS 128
layout( local_size_x = THREADS ) in;


layout(std140, binding=5) buffer _models_instance_count//same size as models_intervals
{
    writeonly uint models_instance_count[];
};

uniform uint count;

void main()
{
    uint cnt = count / THREADS + 1;
    uint st, en;
    st = cnt * gl_LocalInvocationIndex;
    en = cnt * (gl_LocalInvocationIndex + 1);

    en = en > count ? count : en;

    for (uint i = st; i < en; i++)
    {
        models_instance_count[i] = 0;
    }
}
