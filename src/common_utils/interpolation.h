#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
namespace interpolation
{
  template <typename T>
  struct Spline
  {
    // cubic spline y(x) = a + x*(b+x*(c+d*x))
    T a;
    T b;
    T c;
    T d;
    T x;
    T get(T y)
    {
      return a + (y - x) * (b + (y - x) * (c + d * (y - x)));
    }
  };

  template <typename T>
  std::vector<Spline<T>> spline(std::vector<T> &x, std::vector<T> &y)
  {
    int n = x.size() - 1;
    std::vector<T> a;
    a.insert(a.begin(), y.begin(), y.end());
    std::vector<T> b(n);
    std::vector<T> d(n);
    std::vector<T> h(n);

    for (int i = 0; i < n; ++i)
      h[i] = (x[i + 1] - x[i]);

    std::vector<T> alpha(n);
    alpha[0] = 0;
    for (int i = 1; i < n; ++i)
      alpha[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

    std::vector<T> c(n + 1);
    std::vector<T> l(n + 1);
    std::vector<T> mu(n + 1);
    std::vector<T> z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; ++i)
    {
      l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l[i];
      z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n - 1; j >= 0; --j)
    {
      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
      d[j] = (c[j + 1] - c[j]) / 3 / h[j];
    }

    std::vector<Spline<T>> output_set(n);
    for (int i = 0; i < n; ++i)
    {
      output_set[i].a = a[i];
      output_set[i].b = b[i];
      output_set[i].c = c[i];
      output_set[i].d = d[i];
      output_set[i].x = x[i];
    }
    return output_set;
  }
}