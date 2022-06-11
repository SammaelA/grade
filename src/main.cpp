#include "parser.h"
//#include "cities_generator/cities_gen_main.h"

int main(int argc, char *argv[])
{
    //if (argc == 1)
    //    return cities_generator_main();
    //else
        return parser::parser_main(argc, argv);
    //return old_main(argc,argv);
}