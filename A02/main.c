#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <omp.h>
#include <time.h>

#include "mandelbrot_set.h"

int main(int argc, char **argv)
{
	struct timespec begin, end;
	int num_threads = 1;
	int max_iter = 1000; // sizeof(colors);
	int x_resolution = 1000;
	int y_resolution = 1000;
	int palette_shift = 0;
	double view_x0 = -2.00;
	double view_x1 = +1.25;
	double view_y0 = -1.25;
	double view_y1 = +1.25;
	char file_name[256] = "mandelbrot.ppm";
	int no_output = 0;

	double x_stepsize;
	double y_stepsize;

	char batchMode = 0;
    int numTrials = 1;

	int c;
	while ((c = getopt(argc, argv, "A:t:p:i:r:v:n:f:")) != -1)
	{
		switch (c)
		{
		case 'A':
                batchMode = 1;
                if (sscanf(optarg, "%d", &numTrials) != 1)
                {
                    goto error;
                }
                break;

		case 't':
			if (sscanf(optarg, "%d", &num_threads) != 1)
				goto error;
			break;

		case 'p':
			if (sscanf(optarg, "%d", &palette_shift) != 1)
				goto error;
			break;

		case 'i':
			if (sscanf(optarg, "%d", &max_iter) != 1)
				goto error;
			break;

		case 'r':
			if (sscanf(optarg, "%dx%d", &x_resolution, &y_resolution) != 2)
				goto error;
			break;

		case 'v':
			if (sscanf(optarg, "[%lf,%lf]x[%lf,%lf]", &view_x0, &view_x1,
			           &view_y0, &view_y1) != 4)
				goto error;
			break;
		case 'n':
			if (sscanf(optarg, "%d", &no_output) != 1)
				goto error;
			break;
		case 'f':
			strncpy(file_name, optarg, sizeof(file_name));
			break;

		case '?':
error: printf(
			    "Usage:\n"
			    "-A \t num of trials to average\n"
			    "-t \t number of threads used in computation\n"
			    "-i \t maximum number of iterations per pixel\n"
			    "-r \t image resolution to be computed\n"
			    "-v \t view frame that should be computed\n"
			    "-p \t shift palette to change colors\n"
			    "-f \t output file name\n"
				 "-n \t no output(default: 0)\n"
			    "\n"
			    "Example:\n"
			    "%s -t 1 -r 720x480 -i 5000 -v [-2.0,1.0]x[-1.0,1.0] -f mandelbrot.ppm\n",
			    argv[0]);
			exit(EXIT_FAILURE);
			break;
		}
	}

	x_stepsize = (view_x1 - view_x0) / x_resolution;
	y_stepsize = (view_y1 - view_y0) / y_resolution;

	if (batchMode)
        stderr = fopen("/dev/null","w");
    
    fprintf(stderr, "Following settings are used for computation:\n"
           "Threads: %d\n"
           "Max. iterations: %d\n"
           "Resolution: %dx%d\n"
           "View frame: [%lf,%lf]x[%lf,%lf]\n"
           "Stepsize x = %lf y = %lf\n", num_threads, max_iter, x_resolution,
           y_resolution, view_x0, view_x1, view_y0, view_y1, x_stepsize,
           y_stepsize);
    
    if (!no_output)
    {
        fprintf(stderr, "Output file: %s\n", file_name);
    }
    else
    {
        fprintf(stderr, "No output will be writen\n");
    }
    
    FILE *file = NULL;
    if (!no_output)
    {
        if ((file = fopen(file_name, "w")) == NULL)
        {
            perror("fopen");
            exit(EXIT_FAILURE);
        }
    }
    
    unsigned char (*image)[x_resolution][3];
    image = (unsigned char (*)[x_resolution][3]) (malloc(x_resolution * y_resolution * sizeof(char[3])));
    
    if (image == NULL)
    {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

	// compute mandelbrot	
	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (int trial=0; trial<numTrials; ++trial)
    {
		mandelbrot_draw(x_resolution, y_resolution, max_iter, view_x0, view_x1, view_y0, view_y1, x_stepsize, y_stepsize, palette_shift, image, num_threads);
	}
	clock_gettime(CLOCK_MONOTONIC, &end);
    
    if (!no_output)
    {
        if (fprintf(file, "P6\n%d %d %d\n", x_resolution, y_resolution, 255) < 0)
        {
            perror("fprintf");
            exit(EXIT_FAILURE);
        }
        size_t bytes_written = fwrite(image, 1,
                                      x_resolution * y_resolution * sizeof(char[3]), file);
        if (bytes_written < x_resolution * y_resolution * sizeof(char[3]))
        {
            perror("fwrite");
            exit(EXIT_FAILURE);
        }
        fclose(file);
    }
    
    fprintf(stderr, "\n\n");
    double execTimeGlobalSec = ((double) end.tv_sec + 1.0e-9 * end.tv_nsec) -
                         ((double) begin.tv_sec + 1.0e-9 * begin.tv_nsec);
    double execTimeAvgSec = execTimeGlobalSec/numTrials;
    if (batchMode)
        printf("xresolution,%d,yresolution,%d,threads,%d,averagingTrials,%d,timeGlobal[s],%.5f,timeAvg[s],%.5f\n",
               x_resolution,
               y_resolution,
               num_threads,
               numTrials,
               execTimeGlobalSec,
               execTimeAvgSec);
    else
        printf("Time: %.5f seconds\n", execTimeGlobalSec);
    
    free(image);
    return 0;
}
