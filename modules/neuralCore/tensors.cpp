#include "tensors.h"
#include <cassert>
#include <algorithm>
#include <cstdio>
#include <string>
#include <cstring>

namespace nn
{
  TensorView::TensorView(ValueType *data, const std::vector<IndexType> &sizes) : TensorView(data, get_scheme(sizes))
  {
  }

  TensorView::TensorView(ValueType *data, const Scheme &_scheme) : Dim(get_dimensions(_scheme)),
                                                                   _data(data),
                                                                   scheme(_scheme),
                                                                   total_size(get_total_size(_scheme, get_dimensions(_scheme)))
  {
  }

  int get_dimensions(const Scheme &scheme)
  {
    int dims = 0;
    while (dims < MAX_DIM && scheme[dims].size > 0)
      dims++;
    return dims;
  }

  Scheme get_scheme(const std::vector<IndexType> &sizes)
  {
    assert(sizes.size() <= MAX_DIM);
    Scheme s;
    for (int i = 0; i < sizes.size(); i++)
      s[i].size = sizes[i];
    return s;
  }

  IndexType get_total_size(const Scheme &scheme, int Dim)
  {
    if (Dim == 0)
      return 0;
    // TODO: support non-compact tensors
    IndexType total_sz = 1;
    for (int i = 0; i < Dim; i++)
      total_sz *= scheme[i].size;
    return total_sz;
  }

  IndexType get_total_size(const std::vector<IndexType> &sizes)
  {
    IndexType total_sz = 1;
    for (auto sz : sizes)
      total_sz *= sz;
    return total_sz;
  }

  FullIndex linear_to_full_index(const TensorView &t, IndexType index)
  {
    FullIndex full_index;
    for (int i = 0; i < t.Dim; i++)
    {
      full_index[i] = index % t.scheme[i].size;
      index /= t.scheme[i].size;
    }
    return full_index;
  }

  void print_scheme(int Dim, const Scheme &scheme)
  {
    printf("[ ");
    for (int i = 0; i < Dim; i++)
      printf("%u ", scheme[i].size);
    printf("]");
  }

  int compact_dims(const Scheme &scheme)
  {
    // TODO: support non-compact tensors
    return get_dimensions(scheme);
  }

  void print(const TensorView &t)
  {
    FullIndex prev_i({0u});
    std::vector<std::string> delims = {" ", "\n", "\n========\n", "\n\n#####4#####\n\n",
                                       "\n\n#####5#####\n\n", "\n\n#####6#####\n\n", "\n\n#####7#####\n\n"};
    printf("%d-dimentional tensor ", t.Dim);
    print_scheme(t.Dim, t.scheme);
    printf("\n");

    IndexType offset = 0;
    TV_ITERATE(t, 1, offset,
               {
                 const ValueType &val = t.get(offset);
                 FullIndex index = linear_to_full_index(t, offset);
                 int delim = -1;
                 for (int i = 0; i < t.Dim; i++)
                 {
                   if (index[i] != prev_i[i])
                     delim = i;
                 }
                 if (delim >= 0)
                   printf("%s", delims[std::min((int)(delims.size() - 1), delim)].c_str());
                 prev_i = index;

                 if constexpr (std::is_floating_point<ValueType>::value)
                   printf("%8.4f ", (float)val);
                 else
                   printf("%d ", (int)val);
               });
    printf("\n");
  }

  // Creates new tensor view with the same data, but different interpretation
  TensorView reshape(const TensorView &t, const Scheme &scheme)
  {
    int non_compact = t.Dim - compact_dims(t.scheme);

    IndexType new_size = 1;
    int new_dims = get_dimensions(scheme);
    for (int i = 0; i < new_dims; i++)
      new_size *= scheme[i].size;

    // total size shouldn't change
    assert(new_size == t.total_size);

    // should have the same number of non-compact dims
    assert(non_compact == new_dims - compact_dims(scheme));

    // it should preserve non-compact dims layout, i.e.
    //[20,10,<5>] --> [2,10,10,<5>] is ok
    //[10, <20> ] --> [2, ]
    for (int i = 0; i < non_compact; i++)
      assert(t.scheme[t.Dim - i - 1].size == scheme[new_dims - i - 1].size);

    return TensorView(t._data, scheme);
  }
  TensorView reshape(const TensorView &t, const std::vector<IndexType> &sizes)
  {
    return reshape(t, get_scheme(sizes));
  }

  // Returns (Dim-1)-dimentional tensor with given index
  TensorView slice(const TensorView &t, IndexType index)
  {
    assert(t.Dim > 1);
    assert(index < t.scheme[t.Dim - 1].size);

    Scheme new_scheme;
    for (int i = 0; i < t.Dim - 1; i++)
      new_scheme[i] = t.scheme[i];

    // TODO: support non-compact tensors
    IndexType offset = index * t.step(t.Dim - 1);
    return TensorView(t._data + offset, new_scheme);
  }

  // Returns (Dim-1)-dimentional tensor with given range [first, second)
  TensorView slice(const TensorView &t, std::pair<IndexType, IndexType> range)
  {
    assert(t.Dim == compact_dims(t.scheme));
    assert(range.first < range.second);
    assert(range.first < t.scheme[t.Dim - 1].size);
    assert(range.second <= t.scheme[t.Dim - 1].size);

    Scheme new_scheme = t.scheme;
    new_scheme[t.Dim - 1].size = range.second - range.first;
    IndexType offset = range.first * t.step(t.Dim - 1);
    return TensorView(t._data + offset, new_scheme);
  }

  void fill(TensorView &t, ValueType value)
  {
    if (compact_dims(t.scheme) == t.Dim)
      std::fill_n(t._data, t.total_size, value);
    else
    {
      IndexType offset = 0;
      TV_ITERATE(t, 1, offset, { t.get(offset) = value; });
    }
  }

  void copy(const TensorView &t, TensorView out)
  {
    assert(t.total_size == out.total_size);
    assert(compact_dims(t.scheme) == t.Dim);
    assert(compact_dims(out.scheme) == out.Dim);

    memcpy(out._data, t._data, sizeof(ValueType) * t.total_size);
  }

  // Multiplies tensor with Dim>=2 by a vertex
  // if Dim=2 it is standart matrix-vector multiplication
  // if Dim>2 tensor is treated as a (Dim-2)-dimentional array of matrices, each
  // of them is multiplied by a vector.
  // out is a (Dim-1)-dimentional tensor
  void vec_mul(const TensorView &t, const TensorView &v, TensorView out)
  {
    assert(compact_dims(t.scheme) >= 2);
    assert(compact_dims(out.scheme) >= 1);
    assert(compact_dims(v.scheme) >= 1);
    assert(v.scheme[0].size == t.scheme[0].size);
    assert(out.Dim == t.Dim - 1);
    assert(out.scheme[0].size == t.scheme[1].size);
    for (int i = 2; i < t.Dim; i++)
      assert(t.scheme[i].size == out.scheme[i - 1].size);

    IndexType step_size = v.total_size;
    IndexType offset1 = 0;
    IndexType offset2 = 0;
    TV_ITERATE_2(t, t.scheme[0].size, offset1, out, 1, offset2,
                 {
                   out.get(offset2) = 0;
                   for (IndexType i = 0; i < step_size; i++)
                     out.get(offset2) += v.get(i) * t.get(offset1 + i);
                 });
  }

  // Tensor is treated as a (Dim-2)-dimentional array of compact matrices,
  // each of them is transposed
  void transpose(const TensorView &t, TensorView out)
  {
    assert(compact_dims(t.scheme) >= 2);
    assert(compact_dims(out.scheme) >= 2);
    assert(out.Dim == t.Dim);
    assert(t.scheme[0].size == out.scheme[1].size);
    assert(t.scheme[1].size == out.scheme[0].size);
    for (int i = 2; i < t.Dim; i++)
      assert(t.scheme[i].size == out.scheme[i].size);
    assert(t.scheme[0].size == out.scheme[1].size);
    assert(t.scheme[1].size == out.scheme[0].size);

    for (int i = 2; i < t.Dim; i++)
      assert(t.scheme[i].size == out.scheme[i].size);

    IndexType step_size = t.scheme[0].size * t.scheme[1].size;
    IndexType offset1 = 0;
    IndexType offset2 = 0;
    TV_ITERATE_2(t, step_size, offset1, out, step_size, offset2,
                 {
                   for (IndexType i = 0; i < t.scheme[1].size; i++)
                     for (IndexType j = 0; j < t.scheme[0].size; j++)
                       out.get(offset2 + j * t.scheme[1].size + i) = t.get(offset1 + i * t.scheme[0].size + j);
                 });
  }

  // t1 += t2
  void add(TensorView t1, const TensorView &t2)
  {
    add(t1, t2, t1);
  }

  // out = t1+t2
  // if Dim(t1) = Dim(t2) perform element-wise operation
  // else if Dim(t1) > Dim(t2), t1 is treated as an array of Dim(t2)
  // tensors and t2 is added to each of them
  // out should have exactly the same size as t1 by every dimension
  // t2 should be compact
  void add(const TensorView &t1, const TensorView &t2, TensorView out)
  {
    assert(t1.Dim == out.Dim);
    for (int i = 0; i < t1.Dim; i++)
      assert(t1.scheme[i].size == out.scheme[i].size);

    assert(t1.Dim >= t2.Dim);
    assert(compact_dims(t2.scheme) == t2.Dim);
    for (int i = 0; i < t2.Dim; i++)
      assert(t1.scheme[i].size == t2.scheme[i].size);

    IndexType step_size = t2.total_size;
    IndexType offset1 = 0;
    IndexType offset2 = 0;
    TV_ITERATE_2(t1, step_size, offset1, out, step_size, offset2,
                 {
                   for (IndexType i = 0; i < step_size; i++)
                     out.get(offset2 + i) = t1.get(offset1 + i) + t2.get(i);
                 });
  }

  void sum(const TensorView &t, TensorView out, const std::vector<int> &dimensions)
  {
  }

  // performs operation ⊗ between two vectors
  // returns matrix
  // [x, y]⊗[a, b] = [x*a, x*b]
  //                 [y*a, y*b]
  void vector_outer_product(const TensorView &v1, const TensorView &v2, TensorView out)
  {
    assert(v1.Dim == 1);
    assert(compact_dims(v1.scheme) == 1);
    assert(v2.Dim == 1);
    assert(compact_dims(v2.scheme) == 1);
    assert(out.Dim == 2);
    assert(compact_dims(out.scheme) == 2);
    assert(v1.total_size * v2.total_size == out.total_size);

    for (IndexType i = 0; i < v1.total_size; i++)
      for (IndexType j = 0; j < v2.total_size; j++)
        out.get(i * v2.total_size + j) = v1.get(i) * v2.get(j);
  }
}