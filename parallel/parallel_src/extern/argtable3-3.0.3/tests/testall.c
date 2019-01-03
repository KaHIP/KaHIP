#include <stdio.h>

#include "CuTest.h"

CuSuite* get_arglit_testsuite();
CuSuite* get_argstr_testsuite();
CuSuite* get_argint_testsuite();
CuSuite* get_argdate_testsuite();
CuSuite* get_argdbl_testsuite();
CuSuite* get_argfile_testsuite();
CuSuite* get_argrex_testsuite();

void RunAllTests(void)
{
	CuString *output = CuStringNew();
	CuSuite* suite = CuSuiteNew();

	CuSuiteAddSuite(suite, get_arglit_testsuite());
	CuSuiteAddSuite(suite, get_argstr_testsuite());
	CuSuiteAddSuite(suite, get_argint_testsuite());
	CuSuiteAddSuite(suite, get_argdate_testsuite());
	CuSuiteAddSuite(suite, get_argdbl_testsuite());
	CuSuiteAddSuite(suite, get_argfile_testsuite());
	CuSuiteAddSuite(suite, get_argrex_testsuite());

	CuSuiteRun(suite);
	CuSuiteSummary(suite, output);
	CuSuiteDetails(suite, output);
	printf("%s\n", output->buffer);
}

int main(void)
{
	RunAllTests();
    return 0;
}
