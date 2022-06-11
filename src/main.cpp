#include "parser.h"
#include "cities_generator/cities_gen_main.h"
#include "tinyEngine/TinyEngine.h"
View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
int main(int argc, char *argv[])
{
    if (argc == 1)
        return cities_generator_main();
    else
        return parser::parser_main(argc, argv);
    //return old_main(argc,argv);
}