#include <CppUTest/CommandLineTestRunner.h>

int verbose_flag = 0;

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}
