#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"

using namespace boost::gil;

template <typename Out> struct halfdiff_cast_channels; // forward declaration

template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

	// TODO put your solution in here.

}

#endif // !_X_GRADIENT_PAR_
