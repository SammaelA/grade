#include <vector>

typedef float my_float;

class SimpleMeshData
{
  static const int elem_size = 8;
  std::vector<my_float> data;
public:
  std::vector<my_float> get_tringle_data(int i)
  {
    std::vector<my_float> res;
    for (int ind = i * elem_size * 3; ind < data.size() && ind < (i + 1) * elem_size * 3; ++ind)
    {
      res.push_back(data[ind]);
    }
    return res;
  }
  void set_tringle_data(int i, std::vector<my_float> tr)
  {

    for (int ind = i * elem_size * 3, cnt = 0; ind < data.size() && cnt < elem_size * 3 && cnt < tr.size(); ++ind, ++cnt)
    {
      data[ind] = tr[cnt];
    }
  }
  void add_data(std::vector<my_float> tr)
  {
    for (auto i : tr)
    {
      data.push_back(i);
    }
  }
  SimpleMeshData operator|(const SimpleMeshData &elem)
  {
    //or
    return *this;
  }
  SimpleMeshData operator-(const SimpleMeshData &elem)
  {
    //subtract
    return *this;
  }
  SimpleMeshData operator&(const SimpleMeshData &elem)
  {
    //and
    return *this;
  }
};