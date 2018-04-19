#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <pthread.h>
#include <errno.h>

#include "mandelbrot_set.h"

struct static_data
{
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

struct pthread_args
{
	int id;
	int* taskCounter;
	pthread_mutex_t* taskCounterMutex;
	int maxTasks; // Max number of tasks
	int rowStride; // How many rows to cover with 1 task computation
	struct static_data* staticData;
};

int getNextTask(int* counter, pthread_mutex_t* mutex, int stride)
{
	int nextTask;
	pthread_mutex_lock(mutex);
	nextTask = *counter;
	*counter += stride;
	pthread_mutex_unlock(mutex);
	return nextTask;
}

void* mandelbrot_kernel(void* args)
{
	double y;
	double x;

	complex double Z;
	complex double C;

	int k;

	// // Debug - do nothing if thread id is odd
	// if (a->id % 2)
	// 	return NULL;
	// //

	struct pthread_args* a = (struct pthread_args*) args;
	struct static_data* sd = a->staticData;
	// pthread_mutex_lock(a->lock);
	unsigned char (*img)[sd->x_resolution][3]
		= (unsigned char (*)[sd->x_resolution][3]) sd->img;
	int* taskCounter = a->taskCounter;
	pthread_mutex_t* mutex = a->taskCounterMutex;
	int taskStride = a->rowStride;
	int currentTask = getNextTask(taskCounter, mutex, taskStride);
	while (currentTask < a->maxTasks)
	{
		int endTask = currentTask + taskStride;
		for (int i = currentTask; i < endTask; ++i)
		{
			for (int j = 0; j < sd->x_resolution; j++)
			{
				y = sd->view_y1 - i * sd->y_stepsize;
				x = sd->view_x0 + j * sd->x_stepsize;

				Z = 0 + 0 * I;
				C = x + y * I;

				k = 0;

				do
				{
					Z = Z * Z + C;
					k++;
				} while (cabs(Z) < 2 && k < sd->max_iter);

				if (k == sd->max_iter)
				{
					memcpy(img[i][j], "\0\0\0", 3);
				}
				else
				{
					int index = (k + sd->palette_shift)
					            % (sizeof(colors) / sizeof(colors[0]));
					memcpy(img[i][j], colors[index], 3);
				}
			}
		}
		currentTask = getNextTask(taskCounter, mutex, taskStride);
	}
}

void mandelbrot_draw(int x_resolution, int y_resolution, int max_iter,
	                double view_x0, double view_x1, double view_y0, double view_y1,
	                double x_stepsize, double y_stepsize,
	                int palette_shift, unsigned char (*img)[x_resolution][3],
						 int num_threads) {
	// TODO:
	// implement your solution in this file.

	// Create thread-specific data
	pthread_t* threads = (pthread_t*) calloc(num_threads, sizeof(pthread_t));
	struct static_data* staticData = (struct static_data*) calloc(1, sizeof(struct static_data));
	struct pthread_args* args = (struct pthread_args*) calloc(num_threads, sizeof(struct pthread_args));
	int taskCounter = 0;
	pthread_mutex_t* taskCounterMutex = (pthread_mutex_t*) calloc(1, sizeof(pthread_mutex_t));
	pthread_mutex_init(taskCounterMutex, NULL);
	// Initializing static data
	staticData->x_resolution = x_resolution;
	staticData->y_resolution = y_resolution;
	staticData->max_iter = max_iter;
	staticData->view_x0 = view_x0;
	staticData->view_x1 = view_x1;
	staticData->view_y0 = view_y0;
	staticData->view_y1 = view_y1;
	staticData->x_stepsize = x_stepsize;
	staticData->y_stepsize = y_stepsize;
	staticData->palette_shift = palette_shift;
	staticData->img = (unsigned char *) img;
	// Initialize thread args
	for (int i = 0; i<num_threads; ++i)
	{
		args[i].id = i;
		args[i].taskCounter = &taskCounter;
		args[i].taskCounterMutex = taskCounterMutex;
		args[i].maxTasks = staticData->y_resolution;
		args[i].rowStride = 1;
		args[i].staticData = staticData;
		// Actually create the thread
		pthread_create(threads+i, NULL, mandelbrot_kernel, args+i);
	}
	for (int i = 0; i<num_threads; ++i)
	{
		pthread_join(threads[i], NULL);
	}
	// Cleanup
	free(threads);
	free(args);
	free(taskCounterMutex);
}
