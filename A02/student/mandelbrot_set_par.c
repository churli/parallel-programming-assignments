#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <pthread.h>

#include "mandelbrot_set.h"

void mandelbrot_draw(int x_resolution, int y_resolution, int max_iter,
			                double view_x0, double view_x1, double view_y0, double view_y1,
						                double x_stepsize, double y_stepsize,
									                int palette_shift, unsigned char (*image)[x_resolution][3],
																	 int num_threads) {
	// TODO: implement your solution in this file.
}
