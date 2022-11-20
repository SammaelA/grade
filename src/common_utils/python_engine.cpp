#include "python_engine.h"
#include <Python.h>

//Hacky but simple
static bool initialized = false;
namespace python_engine
{
  void init()
  {
    if (!initialized)
    {
      initialized = true;
      Py_Initialize();
    }
  }
}