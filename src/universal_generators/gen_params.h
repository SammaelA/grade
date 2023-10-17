#include <list>

typedef float my_float;

class Params
{
  std::list<my_float> data;
public:
  Params() {}
  Params(const Params &&p)
  {
    for (auto i : p.data)
    {
      data.push_back(i);
    }
  }

  void add(my_float item)
  {
    data.push_back(item);
  }
  
  my_float get()
  {
    my_float res = 0;
    if (data.end() != data.begin())
    {
      res = *data.end();
      data.pop_back();
    }
    return res;
  }

  my_float read()
  {
    my_float res = 0;
    if (data.end() != data.begin())
    {
      res = *data.end();
    }
    return res;
  }
  void read_blk(/*blk file*/)
  {

  }
  void write_blk(/*blk file*/)
  {

  }
};