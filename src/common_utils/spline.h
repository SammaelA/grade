#pragma once
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
struct SplineSet
{
    double a;
    double b;
    double c;
    double d;
    double x;
    double get(double y)
    {
        return a + (y - x)*(b + (y - x)*(c + d*(y -x)));
    }
};
std::vector<SplineSet> spline(std::vector<double> &x, std::vector<double> &y);