#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"

#include <future>

using namespace boost::gil;

template <typename Out> struct halfdiff_cast_channels; // forward declaration

template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

	typedef typename channel_type<DstView>::type dst_channel_t;

	auto H = src.height();
	auto W = src.width();
	std::cout << "Image is of size: " << W << "x" << H << std::endl;
	
	int stride = 6; // Rows per thread
	std::future<void>** futureSet = 
		(std::future<void>**) calloc((H/stride)+1, sizeof(std::future<void>*));

	int n = 0; // Thread num
    for (int y = 0; y < H; y+=stride, ++n)
    {
    	std::cout << n << std::endl;
    	*(futureSet[n]) = std::async(std::launch::async, 
            [&, y]{
		        int lim = y+stride;
            	for (int j=y; j<lim && j<H; ++j)
            	{
			        typename SrcView::x_iterator src_it = src.row_begin(j);
			        typename DstView::x_iterator dst_it = dst.row_begin(j);

			        for (int x = 1; x < W - 1; ++x)
			        {
			    		static_transform(
			    			src_it[x - 1], src_it[x + 1], dst_it[x], 
			    			halfdiff_cast_channels<dst_channel_t>()
			    			);
			        }
			    }
    		}
    	);
    }
    futureSet[n] = NULL; // Ensure stop in next iteration

    for (int i=0; i<H && futureSet[i]!=NULL; ++i)
    {
    	std::cout << "Waiting for thread: " << i << std::endl;
    	futureSet[i]->wait();
    }

    free(futureSet);

}

#endif // !_X_GRADIENT_PAR_
