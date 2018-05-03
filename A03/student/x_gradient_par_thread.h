#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"

#include <thread>
#include <future>

using namespace boost::gil;

template <typename Out> struct halfdiff_cast_channels; // forward declaration

template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

	typedef typename channel_type<DstView>::type dst_channel_t;

	auto H = src.height();
	auto W = src.width();
	std::cout << "Image is of size: " << W << "x" << H << std::endl;
	
	std::thread** threadSet = 
		(std::thread**) calloc(H, sizeof(std::thread*));

	int stride = 6; // Rows per thread
	int n = 0; // Thread num
    for (int y = 0; y < H; y+=stride, ++n)
    {
		threadSet[n] = new std::thread(
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

    for (int i=0; i<H && threadSet[i]!=NULL; ++i)
    {
		threadSet[i]->join();
		delete threadSet[i];
    }

    free(threadSet);

}

#endif // !_X_GRADIENT_PAR_
