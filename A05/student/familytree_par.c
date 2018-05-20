#include "familytree.h"
#include <omp.h>
#include <math.h>

#define T 10

__attribute__((optimize(3), __target__("avx,ssse3,sse4a")))

void traverseDo(tree *node, int depth);
int compute_IQ_local(int data);
int is_prime_local(int n);

void traverse(tree *node, int numThreads){
	omp_set_num_threads(numThreads);
	omp_set_nested(1);
	// omp_set_max_active_levels(5);
	#pragma omp parallel
	{
		#pragma omp single
		{
			traverseDo(node, 0);
		}
	}
}

void traverseDo(tree *node, int depth){
	if(node != NULL){
		#pragma omp task final(depth >= T)
		{
			traverse(node->right, depth+1);
		}
		#pragma omp task final(depth >= T)
		{
			traverse(node->left, depth+1);
		}
		
		node->IQ = compute_IQ_local(node->data);
		genius[node->id] = node->IQ;

	}
}

int compute_IQ_local(int data)
{
	// nothing meaningful, just a function that is kind of computationally expensive
	// such that parallelising makes sense :)
	int i; int sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(i = 0; i < data; i++){
		if(is_prime_local(i))
			sum += 1;
	}
	int modSum = sum % 100;
	return 70 + (modSum)*(modSum)*(modSum)*(modSum)/1000000;
}

int is_prime_local(int n)
{
	int i, flag = 1;
	int sqn = sqrt(n);
	if (n%2==0)
	{
		return 0;
	}
	for(i=3; i<=sqn; i+=2)
	{
		if(n%i==0)
		{
			flag=0;
			break;
		}
	}
	return flag;
}
