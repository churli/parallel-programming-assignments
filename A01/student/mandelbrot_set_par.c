#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <pthread.h>

#include "mandelbrot_set.h"

struct pthread_args
{
	int id;
	int iStart;
	int iEnd;
	int jStart;
	int jEnd;
	int x_resolution;
	int y_resolution;
	int max_iter;
	double view_x0;
	double view_x1;
	double view_y0;
	double view_y1;
	double x_stepsize;
	double y_stepsize;
	int palette_shift;
	unsigned char* img;
};

void* mandelbrot_kernel(void* args)
{
	double y;
	double x;

	complex double Z;
	complex double C;

	int k;

	struct pthread_args* a = (struct pthread_args*) args;
	// // Debug - do nothing if thread id is odd
	// if (a->id % 2)
	// 	return NULL;
	// //
	unsigned char (*img)[a->x_resolution][3]
		= (unsigned char (*)[a->x_resolution][3]) a->img;
	for (int i = a->iStart; i < a->iEnd; i++)
	{
		for (int j = a->jStart; j < a->jEnd; j++)
		{
			y = a->view_y1 - i * a->y_stepsize;
			x = a->view_x0 + j * a->x_stepsize;

			Z = 0 + 0 * I;
			C = x + y * I;

			k = 0;

			do
			{
				Z = Z * Z + C;
				k++;
			} while (cabs(Z) < 2 && k < a->max_iter);

			if (k == a->max_iter)
			{
				memcpy(img[i][j], "\0\0\0", 3);
			}
			else
			{
				int index = (k + a->palette_shift)
				            % (sizeof(colors) / sizeof(colors[0]));
				memcpy(img[i][j], colors[index], 3);
			}
		}
	}
}

void mandelbrot_draw(int x_resolution, int y_resolution, int max_iter,
	                double view_x0, double view_x1, double view_y0, double view_y1,
	                double x_stepsize, double y_stepsize,
	                int palette_shift, unsigned char (*img)[x_resolution][3],
						 int num_threads) {
	// TODO:
	// implement your solution in this file.

	// Create threads
	pthread_t* threads = (pthread_t*) malloc(num_threads * sizeof(pthread_t));
	struct pthread_args* args = (struct pthread_args*) malloc(num_threads * sizeof(struct pthread_args));
	for (int i=0; i<num_threads; ++i)
	{
		int x_subres = x_resolution / num_threads;
		int y_subres = y_resolution / num_threads;
		// Set all args
		args[i].id = i;
		// TODO: find a smarter way to compute blocks
		args[i].iStart = i*y_subres;
		args[i].jStart = 0;
		if (i == num_threads-1)
		{
			args[i].iEnd = y_resolution;
		}
		else
		{
			args[i].iEnd = (i+1)*y_subres;
		}
		args[i].jEnd = x_resolution;
		//
		args[i].x_resolution = x_resolution;
		args[i].y_resolution = y_resolution;
		args[i].max_iter = max_iter;
		args[i].view_x0 = view_x0;
		args[i].view_x1 = view_x1;
		args[i].view_y0 = view_y0;
		args[i].view_y1 = view_y1;
		args[i].x_stepsize = x_stepsize;
		args[i].y_stepsize = y_stepsize;
		args[i].palette_shift = palette_shift;
		args[i].img = (unsigned char *) img;
		// Actually create the thread
		pthread_create(threads+i, NULL, mandelbrot_kernel, args+i);
	}
	for (int i=0; i<num_threads; ++i)
	{
		pthread_join(threads[i], NULL);
	}
	free(threads);
	free(args);
}
