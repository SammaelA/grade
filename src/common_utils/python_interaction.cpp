#include "python_interaction.h"
#include "common_utils/utility.h"
#include <vector>
#include <regex>
#include <fcntl.h>
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/transform.hpp>

bool python_started = false;
std::vector<std::string> split(const std::string& input, const std::string& regex) 
{
    // passing -1 as the submatch index parameter performs splitting
    std::regex re(regex);
    std::sregex_token_iterator
        first{input.begin(), input.end(), re, -1},
        last;
    return {first, last};
}

char *PythonHelper::python_func_get_str(char *val) 
{
    char *ret = NULL;

    // Загрузка объекта get_value из func.py
    pObjct = PyDict_GetItemString(pDict, (const char *) "get_value");
    if (!pObjct) {
        return ret;
    }

    do {
        // Проверка pObjct на годность.
        if (!PyCallable_Check(pObjct)) {
            break;
        }

        pVal = PyObject_CallFunction(pObjct, (char *) "(s)", val);
        if (pVal != NULL) {
            PyObject* pResultRepr = PyObject_Repr(pVal);

            // Если полученную строку не скопировать, то после очистки ресурсов Python её не будет.
            // Для начала pResultRepr нужно привести к массиву байтов.
            printf("%d %d %d\n",PyBytes_Check(pResultRepr), PyBytes_CheckExact(pResultRepr), PyUnicode_Check(pResultRepr));
            ret = strdup(PyBytes_AS_STRING(PyUnicode_AsEncodedString(pResultRepr, "utf-8", "ERROR")));

            Py_XDECREF(pResultRepr);
            Py_XDECREF(pVal);
        } else {
            PyErr_Print();
        }
    } while (0);

    return ret;
}
void PythonHelper::python_func(std::string function_name, std::string input)
{

    // Загрузка объекта get_value из func.py
    pObjct = PyDict_GetItemString(pDict, function_name.c_str());
    if (!pObjct) 
        return;

    // Проверка pObjct на годность.
    if (!PyCallable_Check(pObjct)) 
        return;

    pVal = PyObject_CallFunction(pObjct, (char *) "(s)", input.c_str());
}
/**
 * Получение значения переменной содержащей значение типа int
 */
int PythonHelper::python_func_get_val(char *val) 
{
    int ret = 0;

    // Получить объект с именем val
    pVal = PyDict_GetItemString(pDict, (const char *) val);
    if (!pVal) {
        return ret;
    }

    // Проверка переменной на long
    if (PyLong_Check(pVal)) {
        ret = _PyLong_AsInt(pVal);
    } else {
        PyErr_Print();
    }

    return ret;
}
void PythonHelper::init(std::string _scripts_dir)
{
    if (!python_started)
    {
        Py_Initialize();
        python_started = true;
    }
    scripts_dir = _scripts_dir;
}
void PythonHelper::finish()
{
    Py_Finalize();
}
void PythonHelper::run_script(std::string script_file_name)
{
    run_script(script_file_name, "get_hashes");
}
void PythonHelper::run_script(std::string script_file_name, std::string args)
{
        do 
        {
        // Загрузка модуля sys
        sys = PyImport_ImportModule("sys");
        sys_path = PyObject_GetAttrString(sys, "path");
        // Путь до наших исходников Python
        folder_path = PyUnicode_FromString(scripts_dir.c_str());
        PyList_Append(sys_path, folder_path);

        auto *sys_args = PyList_New(0);
        std::vector<std::string> splitted_args = split(args, "\\s+");
        for (auto &s : splitted_args)
        {
            auto *obj = PyUnicode_FromString(s.c_str());
            PyList_Append(sys_args, obj);
            Py_XDECREF(obj);
        }
        PyObject_SetAttrString(sys, "argv",sys_args);
        Py_XDECREF(sys_args);

        // Загрузка func.py
        pName = PyUnicode_FromString(script_file_name.c_str());
        if (!pName) 
        {
            break;
        }
        static bool block_print = true;
        if (block_print)
        {
            PyRun_SimpleString("import os");
            PyRun_SimpleString("import sys");
            PyRun_SimpleString("sys.stdout = open(os.devnull, 'w')");
            PyRun_SimpleString("sys.stderr = open(os.devnull, 'w')");
            block_print = false;
        }
        script_file_name = scripts_dir+"/" +script_file_name+".py";
        //logerr("sfn %s", script_file_name.c_str());
        FILE *f = fopen(script_file_name.c_str(), "rb");
        PyRun_SimpleFile(f, script_file_name.c_str());
        pModule = PyImport_AddModule("__main__");
        // Загрузить объект модуля
        
        if (!pModule) 
        {
            break;
        }

        // Словарь объектов содержащихся в модуле
        pDict = PyModule_GetDict(pModule);
        if (!pDict) 
        {
            break;
        }
        Py_XDECREF(sys_path);
        Py_XDECREF(sys);
        Py_XDECREF(folder_path);
        Py_XDECREF(pName);
        return;
    } while (0);

    // Печать ошибки
    PyErr_Print();
}

void PythonHelper::finish_script()
{
    //PyObject * poAttrList = PyObject_Dir(pModule);

    //PyObject * poAttrIter = PyObject_GetIter(poAttrList);

    //PyObject * poAttrName;

    //Py_XDECREF(pDict);
}

    void PythonHelper::get_int(std::string &name, int *res)
    {
        pVal = PyDict_GetItemString(pDict, name.c_str());
        if (!pVal) 
            return;

        // Проверка переменной на long
        if (PyLong_Check(pVal)) 
        {
            *res = _PyLong_AsInt(pVal);
        } 
        else 
        {
            PyErr_Print();
        }
    }
    void PythonHelper::get_bytes(std::string &name, int *sz, char **data)
    {
        pVal = PyDict_GetItemString(pDict, name.c_str());
        if (!pVal)
            return;
            // Проверка переменной на long
            if (PyBytes_Check(pVal)) 
            {
                *sz = PyBytes_GET_SIZE(pVal);
                //printf("bytes size %d\n", *sz);
                *data = PyBytes_AS_STRING(pVal);
            } 
            else 
            {
                PyErr_Print();
            }
    }

void PythonHelper::get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, int *sz, char **data)
{
    if (size_x)
        *size_x = 0;
    if (size_y)
        *size_y = 0;
    if (sz)
        *sz = 0;
    std::string bytes_name = name + "_" + "bytes";
    std::string x_name = name + "_" + "size_x";
    std::string y_name = name + "_" + "size_y";

    get_int(x_name, size_x);
    get_int(y_name, size_y);
    get_bytes(bytes_name, sz, data);
}
bool PythonHelper::get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, double **data)
{
    int sz = 0;
    char *raw_data = NULL;
    get_numpy_2d_array_double(name, size_x, size_y, &sz, &raw_data);
    int sx = *size_x;
    int sy = *size_y;

    if (!(sx&&sy&&sz&&data))
    {
        logerr("launched script does not contain bytes representation of array %s", name.c_str());
        return false;
    }

    if (sz != sx*sy*sizeof(double))
    {
        logerr("array %s is not a double array", name.c_str());
        return false;
    }

    *data = new double[sx*sy];
    memcpy(*data, raw_data, sz);
    //logerr("data created %d %d", sx, sy);
    return true;
}
bool PythonHelper::get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, float **data)
{
    int sz = 0;
    char *raw_data = NULL;
    get_numpy_2d_array_double(name, size_x, size_y, &sz, &raw_data);
    int sx = *size_x;
    int sy = *size_y;

    if (!(sx&&sy&&sz&&data))
    {
        logerr("launched script does not contain bytes representation of array %s", name.c_str());
        return false;
    }

    if (sz != sx*sy*sizeof(double))
    {
        logerr("array %s is not a double array", name.c_str());
        return false;
    }

    *data = new float[sx*sy];
    for (int i=0;i<sx*sy;i++)
    {
        double d = 0;
        memcpy(&d, raw_data + sizeof(double)*i, sizeof(double));
        (*data)[i] = d;
    }
    return true;
}

void show_errors()
{
  // get the error details
  std::cout << "An error occurred:" << std::endl;
  PyObject *pExcType, *pExcValue, *pExcTraceback;
  PyErr_Fetch(&pExcType, &pExcValue, &pExcTraceback);
  if (pExcType != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcType);
    std::cout << "- EXC type: " << PyUnicode_AsUTF8(pRepr) << std::endl;
    Py_DecRef(pRepr);
    Py_DecRef(pExcType);
  }
  if (pExcValue != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcValue);
    std::cout << "- EXC value: " << PyUnicode_AsUTF8(pRepr) << std::endl;
    Py_DecRef(pRepr);
    Py_DecRef(pExcValue);
  }
  if (pExcTraceback != NULL)
  {
    PyObject *pRepr = PyObject_Repr(pExcValue);
    std::cout << "- EXC traceback: " << PyUnicode_AsUTF8(pRepr) << std::endl;
    Py_DecRef(pRepr);
    Py_DecRef(pExcTraceback);
  }
}

void PythonHelper::test()
{
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("sys.path.append(\"/home/sammael/grade/scripts\")");
    PyRun_SimpleString("print(sys.path)");
    PyObject *pName, *pModule, *pFunc, *pArgs, *pValue, *pIndex;
    float *data_initial = nullptr;
    float *data_to_transport = nullptr;
    int data_floats = 0;

    pName = PyUnicode_FromString((char*)"emb_test");
    pModule = PyImport_Import(pName);
    if (!pModule)
      show_errors();
    
    PyObject *mitsubaContext;
    {
      PyObject *initFunc, *initArgs, *basePath;
      basePath = PyUnicode_FromString("resources/mitsuba_data/");
      initArgs = PyTuple_Pack(1, basePath);
      initFunc = PyObject_GetAttrString(pModule, (char*)"init");
      mitsubaContext = PyObject_CallObject(initFunc, initArgs);
      Py_DECREF(initFunc);
      Py_DECREF(initArgs);
      Py_DECREF(basePath);
    }

    if (!mitsubaContext)
      show_errors();
    {
      PyObject *func, *args, *params, *params_bytes, *params_name;
      params_name = PyUnicode_FromString("vertex_positions");
      args = PyTuple_Pack(2, mitsubaContext, params_name);
      func = PyObject_GetAttrString(pModule, (char*)"get_params");
      params = PyObject_CallObject(func, args);
      if (!params)
        show_errors();
      auto t = Py_TYPE(params);
      logerr("type %s", t->tp_name);
      params_bytes = PyObject_Bytes(params);
      if (!params_bytes)
        show_errors();
      int sz = PyBytes_Size(params_bytes);
      data_floats = sz/sizeof(float);
      char *data = PyBytes_AsString(params_bytes);
      data_initial = new float[data_floats];
      data_to_transport = new float[data_floats];
      memcpy(data_initial, data, sz);
      logerr("C++: recieve %d floats from Python", data_floats);
      for (int i=0;i<10;i++)
        logerr("%f", (float)(data_initial[i]));
      Py_DECREF(args);
      Py_DECREF(func);
      Py_DECREF(params); 
      Py_DECREF(params_bytes);
      Py_DECREF(params_name); 
    }
    {
      float rotate = 0.0;
      glm::vec3 translate = glm::vec3(0,0.3,0);
      glm::mat4 tr_mat = glm::translate(glm::rotate(rotate, glm::vec3(0,1,0)), translate);
      for (int j = 0;j<data_floats;j+=3)
      {
        glm::vec4 res = tr_mat*glm::vec4(data_initial[j], data_initial[j+1], data_initial[j+2], 1);
        data_to_transport[j] = res.x;
        data_to_transport[j+1] = res.y;
        data_to_transport[j+2] = res.z;
      }

      logerr("C++: send %d floats to Python", data_floats);
      for (int i=0;i<10;i++)
        logerr("%f", (float)(data_to_transport[i]));

      PyObject *func, *args, *params_n, *params_bytes, *params, *params_name;
      params_name = PyUnicode_FromString("vertex_positions");
      params_n = PyLong_FromLong(data_floats);
      params_bytes = PyBytes_FromStringAndSize((const char *)data_to_transport, sizeof(float)*data_floats);
      args = PyTuple_Pack(4, mitsubaContext, params_name, params_bytes, params_n);
      func = PyObject_GetAttrString(pModule, (char*)"set_params");
      params = PyObject_CallObject(func, args);
      
      Py_DECREF(args);
      Py_DECREF(func);
      Py_DECREF(params); 
      Py_DECREF(params_n); 
      Py_DECREF(params_bytes); 
      Py_DECREF(params_name);
    }
    pFunc = PyObject_GetAttrString(pModule, (char*)"render");
    if (!pFunc)
      show_errors();
    for (int i=0;i<1;i++)
    {
      
      pIndex = PyLong_FromLong(i);
      pArgs = PyTuple_Pack(2, pIndex, mitsubaContext);
      pValue = PyObject_CallObject(pFunc, pArgs);
      if (!pValue)
        show_errors();
      double result = PyFloat_AsDouble(pValue);
      Py_DECREF(pValue);
      Py_DECREF(pIndex);
      Py_DECREF(pArgs);
      logerr("%d val %f",i, (float)result);
    }
    Py_DECREF(pFunc);
    Py_DECREF(mitsubaContext);
    Py_DECREF(pModule);

    if (data_initial)
      delete[] data_initial;
    if (data_to_transport)
      delete[] data_to_transport;
    Py_Finalize();
}