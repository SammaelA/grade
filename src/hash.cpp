#include "hash.h"
#include "tinyEngine/utility.h"

void Hash::weight_and_normalize()
{
    double s = 0;
    if (weights.size() != start_points.size())
    {
        logerr("hash weight and start_points size mismatch");
    }
    start_points.push_back(data.size());
    for (int i=0;i<weights.size();i++)
    {
        for (int j=start_points[i];j<start_points[i+1];j++)
        {
            data[j] *= weights[i];
            s += SQR(data[j]);
        }
    }
    s = sqrt(s);
    for (int i=0;i<data.size();i++)
    {
        data[i] /= s;
    }

    weighted = true;
    normalized = true;
}
void Hash::weight()
{
    if (weights.size() != start_points.size())
    {
        logerr("hash weight and start_points size mismatch");
    }
    start_points.push_back(data.size());
    for (int i=0;i<weights.size();i++)
    {
        for (int j=start_points[i];j<start_points[i+1];j++)
        {
            data[j] *= weights[i];
        }
    }
    weighted = true;
}
void Hash::normalize()
{
    double s = 0;

    for (int i=0;i<data.size();i++)
    {
        s += SQR(data[i]);
    }

    s = sqrt(s);
    
    for (int i=0;i<data.size();i++)
    {
        data[i] /= s;
    }
    normalized = true;
}
float Hash::L2_dist(Hash &h1, Hash &h2)
{
    if (!h1.normalized || !h1.weighted)
        h1.weight_and_normalize();
    if (!h2.normalized || !h2.weighted)
        h2.weight_and_normalize();
    
    if (h1.data.size() != h2.data.size())
    {
        logerr("size of feature vectors mismatch %d %d",h1.data.size(), h2.data.size());
        return 1000;
    }
    else
    {
        double a = 0;
        for (int i = 0;i<h1.data.size();i++)
        {
            a += SQR(h1.data[i] - h2.data[i]);
        }
        return sqrt(a);
    }
}
float Hash::L1_dist(Hash &h1, Hash &h2)
{
    if (h1.data.size() != h2.data.size())
    {
        logerr("size of feature vectors mismatch %d %d",h1.data.size(), h2.data.size());
        return 1000;
    }
    else
    {
        double a = 0;
        for (int i = 0;i<h1.data.size();i++)
        {
            a += abs(h1.data[i] - h2.data[i]);
        }
        return a;
    }
}