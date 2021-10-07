#include "python_interaction.h"
#include "../utility.h"

/*
PyObject *PythonHelper::python_init() 
{
    // Инициализировать интерпретатор Python
    Py_Initialize();

    do {
        // Загрузка модуля sys
        sys = PyImport_ImportModule("sys");
        sys_path = PyObject_GetAttrString(sys, "path");
        // Путь до наших исходников Python
        folder_path = PyUnicode_FromString((const char*) ".");
        PyList_Append(sys_path, folder_path);

        // Загрузка func.py
        pName = PyUnicode_FromString("func");
        if (!pName) {
            break;
        }

        // Загрузить объект модуля
        pModule = PyImport_Import(pName);
        if (!pModule) {
            break;
        }

        // Словарь объектов содержащихся в модуле
        pDict = PyModule_GetDict(pModule);
        if (!pDict) {
            break;
        }

        return pDict;
    } while (0);

    // Печать ошибки
    PyErr_Print();
}

void PythonHelper::python_clear() 
{
    // Вернуть ресурсы системе
    Py_XDECREF(pDict);
    Py_XDECREF(pModule);
    Py_XDECREF(pName);

    Py_XDECREF(folder_path);
    Py_XDECREF(sys_path);
    Py_XDECREF(sys);

    // Выгрузка интерпретатора Python
    Py_Finalize();
}*/

/**
 * Передача строки в качестве аргумента и получение строки назад
 */
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
            //ret = strdup(PyUnicode_AS_DATA(pResultRepr));
            ret = strdup(PyBytes_AS_STRING(PyUnicode_AsEncodedString(pResultRepr, "utf-8", "ERROR")));

            Py_XDECREF(pResultRepr);
            Py_XDECREF(pVal);
        } else {
            PyErr_Print();
        }
    } while (0);

    return ret;
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
    Py_Initialize();
    scripts_dir = _scripts_dir;
}
void PythonHelper::finish()
{
    Py_Finalize();
}
void PythonHelper::run_script(std::string script_file_name, std::string args)
{
    //TODO
}
void PythonHelper::run_script(std::string script_file_name)
{
        do 
        {
        // Загрузка модуля sys
        sys = PyImport_ImportModule("sys");
        sys_path = PyObject_GetAttrString(sys, "path");
        // Путь до наших исходников Python
        folder_path = PyUnicode_FromString(scripts_dir.c_str());
        PyList_Append(sys_path, folder_path);

        // Загрузка func.py
        pName = PyUnicode_FromString(script_file_name.c_str());
        if (!pName) 
        {
            break;
        }

        // Загрузить объект модуля
        pModule = PyImport_Import(pName);
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

        return;
    } while (0);

    // Печать ошибки
    PyErr_Print();
}

void PythonHelper::finish_script()
{
    Py_XDECREF(pDict);
    Py_XDECREF(pModule);
    Py_XDECREF(pName);

    Py_XDECREF(folder_path);
    Py_XDECREF(sys_path);
    Py_XDECREF(sys);
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
                printf("bytes size %d\n", *sz);
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