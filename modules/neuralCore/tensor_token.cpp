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
    assert(shape.size() <= TensorProgram::MAX_DIM);
    Dim = shape.size();
    for (int i = 0; i < Dim; i++)
      sizes[i] = shape[i];
    id = tp->add_var(*this);
  }
  TensorToken::TensorToken(const unsigned _sizes[TensorProgram::MAX_DIM]) : 
               TensorToken(_sizes[0], _sizes[1], _sizes[2], _sizes[3], 
                           _sizes[4], _sizes[5], _sizes[6], _sizes[7])
  {
  }
  TensorToken::TensorToken(int sz_0, int sz_1, int sz_2, int sz_3, int sz_4, int sz_5, int sz_6, int sz_7) : 
               TensorToken((unsigned)sz_0, (unsigned)sz_1, (unsigned)sz_2, (unsigned)sz_3,
                           (unsigned)sz_4, (unsigned)sz_5, (unsigned)sz_6, (unsigned)sz_7) {}
  TensorToken::TensorToken(unsigned sz_0, unsigned sz_1, unsigned sz_2, unsigned sz_3, unsigned sz_4, 
                           unsigned sz_5, unsigned sz_6, unsigned sz_7)
  {
    static_assert(TensorProgram::MAX_DIM == 8);

    sizes[0] = sz_0;
    sizes[1] = sz_1;
    sizes[2] = sz_2;
    sizes[3] = sz_3;
    sizes[4] = sz_4;
    sizes[5] = sz_5;
    sizes[6] = sz_6;
    sizes[7] = sz_7;
    Dim = 8;
    for (int i=0;i<8;i++)
    {
      if (sizes[i] == 0)
      {
        Dim = i;
        break;
      }
    }
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
    for (int i = 0; i < TensorProgram::MAX_DIM; i++)
      sizes[i] = other.sizes[i];
    
    //reassigning tensor token. It means that the same variable will have another tensor_id
    if (!same_size)
    {
      assert(!(tp->vars[id].is_input || tp->vars[id].is_output)); //it wil mess up your input/output declarations!
      id = tp->add_var(*this);
      //printf("TensorToken warning: reassigning tensor with different size. It will create new id for the same variable. Mind your step!\n");
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
    if (total_size() == 1 || other.total_size() == 1)
      return;
    if (Dim < other.Dim)
      printf("TensorToken: check failed %u < %u\n", Dim, other.Dim);
    assert(Dim >= other.Dim);
    for (int i = 0; i < other.Dim; i++)
    {
      if (sizes[i] != other.sizes[i])
        printf("TensorToken: Dim %d check failed %u != %u\n", i, sizes[i], other.sizes[i]);
      assert(sizes[i] == other.sizes[i]);
    }
  }

  void TensorToken::g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B, const TensorToken &C,
                          unsigned steps, unsigned step_size, unsigned group_size)
  {
    assert(TensorProgram::cmd_properties[cmd].cls == TensorProgram::CmdClass::ARITHMETICS);
    if (A.total_size() == 1)
      tp->add_command(cmd, A.id, B.id, C.id, 1, B.total_size(), 1, 1);
    else
      tp->add_command(cmd, A.id, B.id, C.id, steps, step_size, group_size, 0);
  }

  TensorToken TensorToken::g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B,
                                 unsigned steps, unsigned step_size, unsigned group_size)
  {
    TensorToken res(A.sizes);
    g_2op(cmd, A, B, res, steps, step_size, group_size);
    return res;
  }

  TensorToken TensorToken::g_2op(TensorProgram::CommandType cmd, const TensorToken &A, const TensorToken &B,
                                 unsigned start_dim)
  {
    assert(start_dim + B.Dim <= A.Dim);
    for (int i=0;i<B.Dim;i++)
      assert(A.sizes[start_dim + i] == B.sizes[i]);
    unsigned group_size = 1;
    for (int i=0;i<start_dim;i++)
      group_size *= A.sizes[i];
    TensorToken res(A.sizes);
    g_2op(cmd, A, B, res, A.total_size()/(B.total_size()*group_size), B.total_size(), group_size);
    return res;
  }

  TensorToken &TensorToken::operator+=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    g_2op(TensorProgram::ADD, *this, other, *this, this->total_size()/other.total_size(), other.total_size(), 1);
    return *this;
  }

  TensorToken TensorToken::operator+(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(total_size() == 1 ? other.sizes : sizes);
    g_2op(TensorProgram::ADD, *this, other, res, this->total_size()/other.total_size(), other.total_size(), 1);
    return res;
  }

  TensorToken &TensorToken::operator*=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    g_2op(TensorProgram::MUL, *this, other, *this, this->total_size()/other.total_size(), other.total_size(), 1);
    return *this;
  }

  TensorToken TensorToken::operator*(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(total_size() == 1 ? other.sizes : sizes);
    g_2op(TensorProgram::MUL, *this, other, res, this->total_size()/other.total_size(), other.total_size(), 1);
    return res;
  }

  TensorToken &TensorToken::operator-=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    g_2op(TensorProgram::SUB, *this, other, *this, this->total_size()/other.total_size(), other.total_size(), 1);
    return *this;
  }

  TensorToken TensorToken::operator-(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(total_size() == 1 ? other.sizes : sizes);
    g_2op(TensorProgram::SUB, *this, other, res, this->total_size()/other.total_size(), other.total_size(), 1);
    return res;
  }

  TensorToken &TensorToken::operator/=(const TensorToken &other)
  {
    check_dimensions_for_arithmetics(other);
    g_2op(TensorProgram::DIV, *this, other, *this, this->total_size()/other.total_size(), other.total_size(), 1);
    return *this;
  }

  TensorToken TensorToken::operator/(const TensorToken &other) const
  {
    check_dimensions_for_arithmetics(other);
    TensorToken res(total_size() == 1 ? other.sizes : sizes);
    g_2op(TensorProgram::DIV, *this, other, res, this->total_size()/other.total_size(), other.total_size(), 1);
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
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
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
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
    res_sizes[0] = total_size() / sizes[Dim - 1];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::O_SUM, id, 0, res.id);
    return res;
  }


  TensorToken TensorToken::minimum(int Dims) const
  {
    if (Dim == 0) // min of scalar is this scalar itself
      return *this;
    if (Dims == -1)
      Dims = Dim;
    assert(Dims > 0);
    assert(Dims <= Dim);
    unsigned res_Dim = Dim - Dims; // remaining dimensions
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < res_Dim; i++)
      res_sizes[i] = sizes[i + Dims];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MINIMUM, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::maximum(int Dims) const
  {
    if (Dim == 0) // max of scalar is this scalar itself
      return *this;
    if (Dims == -1)
      Dims = Dim;
    assert(Dims > 0);
    assert(Dims <= Dim);
    unsigned res_Dim = Dim - Dims; // remaining dimensions
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < res_Dim; i++)
      res_sizes[i] = sizes[i + Dims];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MAXIMUM, id, 0, res.id);
    return res;
  }

  TensorToken TensorToken::transpose(unsigned transp_dim) const
  {
    assert(Dim > 1);
    assert(transp_dim+1 < Dim);
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < TensorProgram::MAX_DIM; i++)
      res_sizes[i] = sizes[i];
    res_sizes[transp_dim] = sizes[transp_dim+1];
    res_sizes[transp_dim+1] = sizes[transp_dim];
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::TRANSP, id, 0, res.id, transp_dim);
    return res;
  }

  TensorToken TensorToken::get(unsigned n) const
  {
    assert(Dim > 0);
    assert(n < sizes[Dim - 1]);

    unsigned res_size = 1;
    unsigned res_Dim = Dim - 1;
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
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
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
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
    assert(new_shape.size() <= TensorProgram::MAX_DIM);

    unsigned res_Dim = new_shape.size();
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
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

  void TensorToken::issue_command(TensorProgram::CommandType type, const TensorToken &A, const TensorToken &B, const TensorToken &C, 
                                  unsigned arg0, unsigned arg1, unsigned arg2, unsigned arg3, unsigned arg4)
  {
    tp->add_command(type, A.id, B.id, C.id, arg0, arg1, arg2, arg3, arg4);
  }

  void TensorToken::fill(float val)
  {
    tp->add_command(TensorProgram::FILL, 0, 0, id, *((unsigned *)(&val)));
  }

  TensorToken TensorToken::add_padding(unsigned left_pad, unsigned right_pad, int pad_Dim) const
  {
    assert(Dim > pad_Dim);
    unsigned pad_mult = 1;
    for (int i=0;i<pad_Dim;i++)
      pad_mult *= sizes[i];
    
    unsigned res_sizes[TensorProgram::MAX_DIM];
    for (int i = 0; i < TensorProgram::MAX_DIM; i++)
      res_sizes[i] = sizes[i];
    res_sizes[pad_Dim] = sizes[pad_Dim] + left_pad + right_pad;
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::PAD, id, 0, res.id, total_size()/(pad_mult*sizes[pad_Dim]), 
                    pad_mult*sizes[pad_Dim], pad_mult*left_pad, pad_mult*right_pad);
    return res;
  }

  TensorToken TensorToken::flip(unsigned axis) const
  {
    assert(Dim > axis);
    unsigned res_sizes[TensorProgram::MAX_DIM];
    for (int i = 0; i < TensorProgram::MAX_DIM; i++)
      res_sizes[i] = sizes[i];
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::FLIP, id, 0, res.id, axis);
    return res;
  }

  void TensorToken::random(unsigned seed)
  {
    tp->add_command(TensorProgram::URAND, 0, 0, id, seed);
  }

  TensorToken TensorToken::vector_outer_product(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim >= 1);
    assert(A.Dim < TensorProgram::MAX_DIM);
    assert(B.Dim == A.Dim);
    for (int i = 1; i < A.Dim; i++)
      assert(A.sizes[i] == B.sizes[i]);

    unsigned res_Dim = A.Dim + 1;
    unsigned res_sizes[TensorProgram::MAX_DIM] = {0, 0, 0, 0, 0, 0, 0, 0};
    res_sizes[0] = B.sizes[0];
    res_sizes[1] = A.sizes[0];
    for (int i = 2; i < res_Dim; i++)
      res_sizes[i] = A.sizes[i - 1];

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::OUTER_P, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::mat_mul_t(const TensorToken &A, const TensorToken &B)
  {
    assert(A.Dim == 2);
    assert(B.Dim == 1 || B.Dim == 2);
    assert(A.sizes[0] == B.sizes[0]);

    unsigned res_Dim = B.Dim;
    unsigned res_sizes[TensorProgram::MAX_DIM] = {B.sizes[1], A.sizes[1], 0, 0};
    if (B.Dim == 1)
    {
      res_sizes[0] = A.sizes[1];
      res_sizes[1] = 0;
    }
    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::MATMUL_T, A.id, B.id, res.id);
    return res;
  }

  TensorToken TensorToken::pow(const TensorToken &A, const TensorToken &B)
  {
    A.check_dimensions_for_arithmetics(B);
    TensorToken res(A.sizes);
    g_2op(TensorProgram::POW, A, B, res, A.total_size()/B.total_size(), B.total_size(), 1);
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

  TensorToken TensorToken::sqrt(const TensorToken &A)
  {
    TensorToken res(A.sizes);
    tp->add_command(TensorProgram::SQRT, A.id, 0, res.id);
    return res;
  }

  /*
  if kernel.Dim = 2 (HxW) then A is treated as an array of matrices (NxiHxiW) and conv2D returns an array (NxoHxoW) 
  if kernel.Dim = 3 (KxHxW) then A is treated as an array of K-channel images (NxKxiHxiW) and conv2D returns an array of matrices (NxoHxoW) 
  if kernel.Dim = 4 (LxKxHxW) then A is treated as an array of K-channel images (NxKxiHxiW) and conv2D returns an array of L-channel images (NxLxoHxoW)
  in all cases N=0 is allowed, which reduces Dim of result by 1 (i.e. (iHxiW) -> (oHxoW))
  padding is not applied, borders are ignored 
  */
  TensorToken TensorToken::conv2D(const TensorToken &A, const TensorToken &kernel, unsigned stride)
  {
    assert(stride > 0);
    assert(kernel.Dim >= 2 && kernel.Dim <= 4);
    if (kernel.Dim == 2)
      assert(A.Dim >= 2);
    else
    {
      assert(A.Dim >= 3);
      assert(A.sizes[2] == kernel.sizes[2]); //number of channels
    }
    unsigned oW = (A.sizes[0] - kernel.sizes[0])/stride + 1;
    unsigned oH = (A.sizes[1] - kernel.sizes[1])/stride + 1;

    unsigned res_sizes[TensorProgram::MAX_DIM] = {0,0,0,0};
    res_sizes[0] = oW;
    res_sizes[1] = oH;
    if (kernel.Dim == 2)
    {
      for (int i=2;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i] = A.sizes[i];
    }
    else if (kernel.Dim == 3)
    {
      for (int i=3;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i-1] = A.sizes[i];
    }
    else //if (kernel.Dim == 4)
    {
      res_sizes[2] = kernel.sizes[3];
      for (int i=3;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i] = A.sizes[i];
    }

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::CONV_2D, A.id, kernel.id, res.id, stride);
    return res;
  }
  
  /*
  if kernel.Dim = 3 (DxHxW) then A is treated as an array of 3D tensors (NxiDxiHxiW) and conv3D returns an array (NxoDxoHxoW) 
  if kernel.Dim = 4 (KxDxHxW) then A is treated as an array of K-channel voxel grids (NxKxiDxiHxiW) and conv3D returns an array of 3D tensors (NxxoDoHxoW) 
  if kernel.Dim = 5 (LxKxDxHxW) then A is treated as an array of K-channel voxel grids (NxKxiDxiHxiW) and conv3D returns an array of L-channel voxel grids (NxLxoDxoHxoW)
  in all cases N=0 is allowed, which reduces Dim of result by 1 (i.e. (iDxiHxiW) -> (oDxoHxoW))
  padding is not applied, borders are ignored 
  */
  TensorToken TensorToken::conv3D(const TensorToken &A, const TensorToken &kernel, unsigned stride)
  {
    assert(stride > 0);
    assert(kernel.Dim >= 3 && kernel.Dim <= 5);
    if (kernel.Dim == 3)
      assert(A.Dim >= 3);
    else
    {
      assert(A.Dim >= 4);
      assert(A.sizes[3] == kernel.sizes[3]); //number of channels
    }
    unsigned oW = (A.sizes[0] - kernel.sizes[0])/stride + 1;
    unsigned oH = (A.sizes[1] - kernel.sizes[1])/stride + 1;
    unsigned oD = (A.sizes[2] - kernel.sizes[2])/stride + 1;

    unsigned res_sizes[TensorProgram::MAX_DIM] = {0,0,0,0,0,0,0,0};
    res_sizes[0] = oW;
    res_sizes[1] = oH;
    res_sizes[2] = oD;
    if (kernel.Dim == 3)
    {
      for (int i=3;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i] = A.sizes[i];
    }
    else if (kernel.Dim == 4)
    {
      for (int i=4;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i-1] = A.sizes[i];
    }
    else //if (kernel.Dim == 5)
    {
      res_sizes[3] = kernel.sizes[4];
      for (int i=4;i<TensorProgram::MAX_DIM;i++)
        res_sizes[i] = A.sizes[i];
    }

    TensorToken res(res_sizes);
    tp->add_command(TensorProgram::CONV_3D, A.id, kernel.id, res.id, stride);
    return res;
  }

}