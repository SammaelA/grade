#include "python_interaction.h"
#include "../utility.h"
#include <vector>
#include <regex>
#include <fcntl.h>

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
}

void PythonHelper::test()
{
    puts("Test func:");
    init("./scripts/deep_hashing");
    run_script("get_hashes");
    int x,y;
    double *res = nullptr;
    bool status = get_numpy_2d_array_double("arr", &x, &y, &res);
    logerr("%d %d xy",x ,y);
    if (status && res)
    {
        for (int i = 0;i<y;i++)
        {
            for (int j=0;j<x;j++)
            {
                debug("%f ",(float)(res[x*i+j]));
            }
            debugnl();
        }
    }
    /*char *str2 = python_func_get_str("Hello from Python!");
    int ret = 0;
    char *str = nullptr;
    // Получить объект с именем val
    pVal = PyDict_GetItemString(pDict, (const char *)"d");

    // Проверка переменной на long
    if (PyBytes_Check(pVal)) {
        int sz = PyBytes_GET_SIZE(pVal);
        printf("bytes size %d\n", sz);
        str = PyBytes_AS_STRING(pVal);
    } else {
        PyErr_Print();
    }

    int i=0;
    if (str)
    {
        while (i < 32)
        {
            printf("%d ",(int)str[i]);
            i++;
        }
    }
    printf("\n");
    puts("Strings:");
    printf("\tString: %s\n", str2);

    puts("Attrs:");
    printf("\ta: %d\n", python_func_get_val("a"));
    printf("\tb: %d\n", python_func_get_val("b"));
    printf("\tc: %d\n", python_func_get_val("c"));*/
    finish_script();
    finish();
}