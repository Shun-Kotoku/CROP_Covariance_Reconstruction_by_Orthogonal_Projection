#ifdef _WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include "Memory.h"
#include "Initialization.h"

void runMemoryCalculation(int argc, char** argv) {
	Data data;
	initializeVariables(data);
	calcNonlinearMemory(data);
	cout << "\n";
}


int main(int argc, char** argv)
{
	runMemoryCalculation(argc, argv);
	cout << "Done...\n";
	return 0;
}
