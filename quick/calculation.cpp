#include <cstdlib>
#include <stdio.h>
using namespace std;
// comment
int main(int argc, char** argv) {
	if(argc<3){
		printf("Usage: %s p n \n", argv[0]);
      		return 0;
    	}

	const double p = atof(argv[1]);
	const double n = atof(argv[2]);

	double s = 1./(1. - p + p/n);
	//printf("n = %f\n", n);
	//printf("p = %f\n", p);	
	printf("s = %f\n", s);
}
