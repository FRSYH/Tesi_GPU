#include<stdio.h>
#include<stdlib.h>


int mod_long_GPU(long long n, long long p) {
	long long v = n, x = 0;

	if (v >= p) {
		v = n % p;
	}
	else {
		if (v < 0) {
			x = n / p;
			v = n - (x*p);
			v += p;
		}
	}
	int r = v;
	return r;
}

int sub_mod_GPU(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa-bb),p);
}

int mul_mod_GPU(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa*bb),p);
}


int main(void){

int module = 49069;
int n = 49068;

int r = mul_mod_GPU(n, n, module);
int r2 = mul_mod_GPU(r, n, module);
printf("r = %d\n", r);
printf("r2 = %d\n", r2);

int m = mod_long_GPU( (long long) n*n, (long long) module);
printf("m = %d\n", m);


}