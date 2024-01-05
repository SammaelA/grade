#include "tensor_compiler.h"
#include <cassert>
#include <map>
#include <string>

namespace nn
{
  TensorToken::TensorToken()
  {
    Dim = 0;
    id = tp->add_var(*this);
  }
  TensorToken::TensorToken(float val)
  {
    Dim = 0;
    sizes[0] = 1;
    id = tp->add_var(*this);
    tp->ftt(id, val);
  }
  TensorToken::TensorToken(const std::vector<unsigned> &shape)
  {
    assert(shape.size() <= TensorCompiler::MAX_DIM);
    Dim = shape.size();
    for (int i = 0; i < Dim; i++)
      sizes[i] = shape[i];
    id = tp->add_var(*this);
  }
  TensorToken::TensorToken(const unsigned _sizes[TensorCompiler::MAX_DIM]) : TensorToken(_sizes[0], _sizes[1], _sizes[2], _sizes[3])
  {
  }
  TensorToken::TensorToken(int sz_0, int sz_1, int sz_2, int sz_3) : TensorToken((unsigned)sz_0, (unsigned)sz_1, (unsigned)sz_2, (unsigned)sz_3) {}
  TensorToken::TensorToken(unsigned sz_0, unsigned sz_1, unsigned sz_2, unsigned sz_3)
  {
    static_assert(TensorCompiler::MAX_DIM == 4);
    Dim = sz_0 == 0 ? 0 : (sz_1 == 0 ? 1 : (sz_2 == 0 ? 2 : (sz_3 == 0 ? 3 : 4)));
    sizes[0] = sz_0;
    sizes[1] = sz_1;
    sizes[2] = sz_2;
    sizes[3] = sz_3;
    id = tp->add_var(*this);
  }
  TensorToken::TensorToken(const TensorToken &other)
  {
    Dim = other.Dim;
    for (int i = 0; i < Dim; i++)
      sizes[i] = other.sizes[i];
    id = tp->add_var(*this);
    tp->add_command(TensorProgram::MOV, other.id, 0, id);
  }
  TensorToken &TensorToken::operator=(const TensorToken &other)
  {
    bool same_size = (Dim == other.Dim);
    if (same_size)
    {
      for (int i = 0; i < Dim; i++)
        same_size = same_size && (sizes[i] == other.sizes[i]);
    }
    Dim = other.Dim;
    for (int i = 0; i < TensorCompiler::MAX_DIM; i++)
      sizes[i] = other.sizes[i];
    if (!same_size) // TODO: do something more clear for end user
    {
      id = tp->add_var(*this);
      printf("TensorToken warning: reassigning tensor with different size. It will create new id for the same variable. Mind your step!\n");
    }
    tp->add_command(TensorProgram::MOV, other.id, 0, id);
    return *this;
  }

  unsigned TensorToken::total_size() const
  {
    unsigned size = 1;
    for (int i = 0; i < Dim; i++)
      size *= sizes[i];
    return size;
  }

  void TensorToken::check_dimensions_for_arithmetics(const TensorToken &other) const
  {
    if (other.total_size() == 1)
      return;
    if (Dim < other.Dim)
      printf("TensorToken: check failed %u < %u\n", Dim, other.Dim);
    assert(Dim >= other.Dim);
    for (int i = 0; i < other.Dim; i++)
    {
      if (sizes[i] != other.sizes[i])
        printf("TensorToken: check failed %u != %u\n", sizes[i], other.sizes[i]);
      assert(sizes[i] == other.sizes[i]);
    }
  }

  TensorToken &TensorToken::operator+=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    tp->add_command(TensorProgram::ADD, id, other.id, id);
    return *this;
  }

  TensorToken TensorToken::operator+(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(sizes);
    tp->add_command(TensorProgram::ADD, id, other.id, res.id);
    return res;
  }

  TensorToken &TensorToken::operator*=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    tp->add_command(TensorProgram::MUL, id, other.id, id);
    return *this;
  }

  TensorToken TensorToken::operator*(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(sizes);
    tp->add_command(TensorProgram::MUL, id, other.id, res.id);
    return res;
  }

  TensorToken &TensorToken::operator-=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    tp->add_command(TensorProgram::SUB, id, other.id, id);
    return *this;
  }

  TensorToken TensorToken::operator-(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(sizes);
    tp->add_command(TensorProgram::SUB, id, other.id, res.id);
    return res;
  }

  TensorToken &TensorToken::operator/=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    tp->add_command(TensorProgram::DIV, id, other.id, id);
    return *this;
  }

  TensorToken TensorToken::operator/(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(sizes);
    tp->add_command(TensorProgram::DIV, id, other.id, res.id);
    return res;
  }

  TensorToken TensorToken::sum(int Dims) const
  {
    if (Dim == 0) // sum of scalar is this scalar itself
      return *this;
    if (Dims == -1)
      Dims = Dim;
    assert(Dims > 0);
    assert(Dims <= Dim);
    unsigned res_Dim = Dim - Dims; // remaining dimensions
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    for (int i = 0; i < res_Dim; i++)
      res_sizes[i] = sizes[i + Dims];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::SUM, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::outer_sum() const
  {
    if (Dim == 0) // sum of scalar is this scalar itself
      return *this;
    unsigned res_Dim = 1;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = total_size() / sizes[Dim - 1];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::O_SUM, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::transpose() const
  {
    assert(Dim > 1);
    unsigned res_Dim = Dim;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = sizes[1];
    res_sizes[1] = sizes[0];
    for (int i = 2; i < res_Dim; i++)
      res_sizes[i] = sizes[i];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::TRANSP, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::get(unsigned n) const
  {
    assert(Dim > 0);
    assert(n < sizes[Dim - 1]);

    unsigned res_size = 1;
    unsigned res_Dim = Dim - 1;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    for (int i = 0; i < res_Dim; i++)
    {
      res_sizes[i] = sizes[i];
      res_size *= res_sizes[i];
    }
    TensorToken res(res_sizes);

    tp->add_command(TensorProgram::COPY, id, 0, res.id, n * res_size, 0, res_size);
    return res;
  }

  TensorToken TensorToken::get(std::pair<unsigned, unsigned> range) const
  {
    unsigned from = range.first;
    unsigned to = range.second;
    assert(Dim > 0);
    assert(from < to);
    assert(to <= sizes[Dim - 1]);

    unsigned res_size = 1;
    unsigned res_Dim = Dim;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    for (int i = 0; i < Dim - 1; i++)
    {
      res_sizes[i] = sizes[i];
      res_size *= res_sizes[i];
    }
    res_sizes[Dim - 1] = to - from;
    TensorToken res(res_sizes);

    tp->add_command(TensorProgram::COPY, id, 0, res.id, from * res_size, 0, (to - from) * res_size);
    return res;
  }

  void TensorToken::set(unsigned n, const TensorToken &t)
  {
    assert(t.Dim == Dim - 1);
    assert(n <= sizes[Dim - 1]);
    unsigned sub_size = 1;
    for (int i = 0; i < Dim - 1; i++)
    {
      assert(sizes[i] == t.sizes[i]);
      sub_size *= sizes[i];
    }
    tp->add_command(TensorProgram::COPY, t.id, 0, id, 0, n * sub_size, sub_size);
  }

  void TensorToken::set_flatten(std::pair<unsigned, unsigned> range, const TensorToken &t)
  {
    unsigned from = range.first;
    unsigned to = range.second;
    assert(t.total_size() == to - from);
    assert(from < to);
    assert(to <= total_size());

    tp->add_command(TensorProgram::COPY, t.id, 0, id, from, 0, (to - from));
  }

  void TensorToken::copy_to(std::pair<unsigned, unsigned> to_range, const TensorToken &t, std::pair<unsigned, unsigned> from_range)
  {
    assert(to_range.second <= total_size());
    assert(from_range.second <= t.total_size());
    assert(to_range.second - to_range.first == from_range.second - from_range.first);
    tp->add_command(TensorProgram::COPY, t.id, 0, id, from_range.first, to_range.first, to_range.second - to_range.first);
  }

  void TensorToken::set(std::pair<unsigned, unsigned> range, const TensorToken &t)
  {
    unsigned from = range.first;
    unsigned to = range.second;
    assert(t.Dim == Dim);
    assert(from < to);
    assert(to <= sizes[Dim - 1]);
    assert(to - from <= t.sizes[Dim - 1]);
    unsigned sub_size = 1;
    for (int i = 0; i < Dim - 1; i++)
    {
      assert(sizes[i] == t.sizes[i]);
      sub_size *= sizes[i];
    }
    tp->add_command(TensorProgram::COPY, t.id, 0, id, 0, from * sub_size, (to - from) * sub_size);
  }

  TensorToken TensorToken::reshape(std::vector<unsigned> new_shape) const
  {
    unsigned new_size = 1;
    for (unsigned s : new_shape)
      new_size *= s;
    unsigned size = 1;
    for (int i = 0; i < Dim; i++)
      size *= sizes[i];

    assert(new_size == size);
    assert(new_shape.size() <= TensorCompiler::MAX_DIM);

    unsigned res_Dim = new_shape.size();
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    for (int i = 0; i < new_shape.size(); i++)
      res_sizes[i] = new_shape[i];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MOV, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::flatten() const
  {
    unsigned size = 1;
    for (int i = 0; i < Dim; i++)
      size *= sizes[i];
    return reshape({size});
  }

  void TensorToken::fill(float val)
  {
    tp->add_command(TensorProgram::FILL, 0, 0, id, *((unsigned *)(&val)));
  }

  TensorToken TensorToken::vector_outer_product(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim >= 1);
    assert(A.Dim < TensorCompiler::MAX_DIM);
    assert(B.Dim == A.Dim);
    for (int i = 1; i < A.Dim; i++)
      assert(A.sizes[i] == B.sizes[i]);

    unsigned res_Dim = A.Dim + 1;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = B.sizes[0];
    res_sizes[1] = A.sizes[0];
    for (int i = 2; i < res_Dim; i++)
      res_sizes[i] = A.sizes[i - 1];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::OUTER_P, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::vector_outer_product_sum(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim >= 1);
    assert(A.Dim < TensorCompiler::MAX_DIM);
    assert(B.Dim == A.Dim);
    for (int i = 1; i < A.Dim; i++)
      assert(A.sizes[i] == B.sizes[i]);

    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = B.sizes[0];
    res_sizes[1] = A.sizes[0];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::OUTER_PS, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::mat_mul_t(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim == 2);
    assert(B.Dim == 2);
    assert(A.sizes[0] == B.sizes[0]);

    unsigned res_Dim = 2;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = A.sizes[1];
    res_sizes[1] = B.sizes[1];
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MATMUL_T, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::mat_vec_mul(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim == 2);
    assert(B.Dim >= 1);
    assert(A.sizes[0] == B.sizes[0]);

    unsigned res_Dim = B.Dim;
    unsigned res_sizes[TensorCompiler::MAX_DIM] = {0, 0, 0, 0};
    res_sizes[0] = A.sizes[1];
    for (int i = 1; i < B.Dim; i++)
      res_sizes[i] = B.sizes[i];
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MATMUL_T, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::pow(const TensorToken &A, const TensorToken &B)
  {
    assert(B.total_size() == 1);
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::POW, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::exp(const TensorToken &A)
  {
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::EXP, A.id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::sin(const TensorToken &A)
  {
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::SIN, A.id, 0, res.id);
    return res;
  }
  
  TensorToken TensorToken::cos(const TensorToken &A)
  {
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::COS, A.id, 0, res.id);
    return res;
  }
  
  TensorToken TensorToken::log(const TensorToken &A)
  {
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::LOG, A.id, 0, res.id);
    return res;
  }
  
}