#pragma once

#include <string>
#include <Python.h>

class PythonHelper
{
public:
    PythonHelper() {};
    void init(std::string _scripts_dir = "./scripts");
    void finish();
    void run_script(std::string script_file_name, std::string args);
    void run_script(std::string script_file_name);
    bool get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, double **data);
    bool get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, float **data);
    void finish_script();
    void test();
private:
    PyObject *python_init();
    void python_clear();
    void python_func(std::string function_name, std::string input);
    char *python_func_get_str(char *val);
    int python_func_get_val(char *val);
    void get_int(std::string &name, int *val);
    void get_bytes(std::string &name, int *sz, char **data);
    void get_numpy_2d_array_double(std::string name, int *size_x, int *size_y, int *sz, char **data);

    std::string scripts_dir;

    PyObject *pName = NULL, *pModule = NULL;
    PyObject *pDict = NULL, *pObjct = NULL, *pVal = NULL;
    PyObject* sys = NULL;
    PyObject* sys_path = NULL;
    PyObject* folder_path = NULL;
};