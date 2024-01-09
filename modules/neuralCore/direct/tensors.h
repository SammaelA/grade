#pragma once
#include <vector>
#include <cstdint>
#include <array>

namespace nnd
{
  using AxisIdType = uint8_t;
  using IndexType = uint32_t;
  using Shape = std::vector<IndexType>;
  using ValueType = float;
  constexpr int MAX_DIM = 4;
  struct DimInfo
  {
    IndexType size = 0;
  };
  using Scheme = std::array<DimInfo, MAX_DIM>;
  using FullIndex = std::array<IndexType, MAX_DIM>;

  struct TensorView
  {
    int Dim = 0;
    ValueType *_data = nullptr;
    Scheme scheme;
    IndexType total_size = 0;

    TensorView() = default; // empty view
    TensorView(ValueType *data, const Scheme &_scheme);
    TensorView(ValueType *data, const std::vector<IndexType> &sizes);
    ~TensorView() = default;
    TensorView(const TensorView &other) = default;
    TensorView(TensorView &&other) = default;
    TensorView &operator=(const TensorView &other) = default;
    TensorView &operator=(TensorView &&other) = default;

    inline IndexType size(int dimension) const
    {
      return scheme[dimension].size;
    }

    inline IndexType step(int dimension) const
    {
      IndexType s = 1;
      for (int i = 0; i < dimension; i++)
        s *= scheme[i].size;
      return s;
    }

    inline const ValueType &get(IndexType i) const
    {
      return _data[i];
    }
    inline ValueType &get(IndexType i)
    {
      return _data[i];
    }

    // TODO: support non-compact tensors
    inline const ValueType &get(IndexType column, IndexType row) const
    {
      return _data[step(1) * row + column];
    }
    inline ValueType &get(IndexType column, IndexType row)
    {
      return _data[step(1) * row + column];
    }

    // TODO: support non-compact tensors
    inline const ValueType &get(IndexType column, IndexType row, IndexType layer) const
    {
      return _data[step(2) * layer + step(1) * row + column];
    }
    inline ValueType &get(IndexType column, IndexType row, IndexType layer)
    {
      return _data[step(2) * layer + step(1) * row + column];
    }

    // TODO: support non-compact tensors
    inline const ValueType &get(IndexType column, IndexType row, IndexType layer, IndexType group) const
    {
      return _data[step(3) * group + step(2) * layer + step(1) * row + column];
    }
    inline ValueType &get(IndexType column, IndexType row, IndexType layer, IndexType group)
    {
      return _data[step(3) * group + step(2) * layer + step(1) * row + column];
    }
  };

  int get_dimensions(const Scheme &scheme);
  Scheme get_scheme(const std::vector<IndexType> &sizes);
  IndexType get_total_size(const Scheme &scheme, int Dim);
  IndexType get_total_size(const std::vector<IndexType> &sizes);
  FullIndex linear_to_full_index(const TensorView &t, IndexType index);
  int compact_dims(const Scheme &scheme);

  void print(const TensorView &view);
  TensorView reshape(const TensorView &t, const std::vector<IndexType> &sizes);
  TensorView reshape(const TensorView &t, const Scheme &scheme);
  TensorView slice(const TensorView &t, IndexType index);
  TensorView slice(const TensorView &t, std::pair<IndexType, IndexType> range);

  void fill(TensorView &t, ValueType value);
  void copy(const TensorView &t, TensorView out);
  void vec_mul(const TensorView &t, const TensorView &v, TensorView out);
  void transpose(const TensorView &t, TensorView out);
  void sum(const TensorView &t, TensorView out, const std::vector<int> &dimensions);
  void add(TensorView t1, const TensorView &t2);
  void add(const TensorView &t1, const TensorView &t2, TensorView out);
  void vector_outer_product(const TensorView &v1, const TensorView &v2, TensorView out);

#define FOR_EACH(t, out, Func)                        \
  {                                                   \
    assert(t.Dim == out.Dim);                         \
    for (int i = 0; i < t.Dim; i++)                   \
      assert(t.scheme[i].size == out.scheme[i].size); \
    for (IndexType i = 0; i < t.total_size; i++)      \
      out.get(i) = Func(t.get(i));                    \
  }

#define FOR_EACH_INPLACE(t, Func) FOR_EACH(t, t, Func)

// TODO: support non-compact tensors
#define TV_ITERATE(t, step_size, offset, F)                    \
  {                                                            \
    IndexType steps_count = t.total_size / step_size;          \
    for (IndexType __step = 0; __step < steps_count; __step++) \
    {                                                          \
      offset = __step * step_size;                             \
      F                                                        \
    }                                                          \
  }

// TODO: support non-compact tensors
#define TV_ITERATE_2(t1, step_size1, offset1, t2, step_size2, offset2, F) \
  {                                                                       \
    IndexType steps_count = t1.total_size / step_size;                    \
    for (IndexType __step = 0; __step < steps_count; __step++)            \
    {                                                                     \
      offset1 = __step * step_size1;                                      \
      offset2 = __step * step_size2;                                      \
      F                                                                   \
    }                                                                     \
  }
}