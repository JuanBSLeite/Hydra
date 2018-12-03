/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2016 - 2018 Antonio Augusto Alves Junior
 *
 *   This file is part of Hydra Data Analysis Framework.
 *
 *   Hydra is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Hydra is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Hydra.  If not, see <http://www.gnu.org/licenses/>.
 *
 *---------------------------------------------------------------------------*/

/*
 * FFTW.h
 *
 *  Created on: 13/11/2018
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef FFTW_H_
#define FFTW_H_

/**
 *
 */
#include <hydra/detail/FFTPolicy.h>
#include<hydra/detail/fftw/WrappersFFTW.h>
#include<hydra/detail/fftw/BaseFFTW.h>
#include<hydra/detail/fftw/ComplexToRealFFTW.h>
#include<hydra/detail/fftw/RealToComplexFFTW.h>
#include<hydra/detail/fftw/ComplexToComplexFFTW.h>

namespace hydra {

template<typename T>
struct detail::FFTPolicy<T, detail::FFTW>
{
	typedef ComplexToComplexFFTW<T> C2C;
	typedef    RealToComplexFFTW<T> R2C;
	typedef    ComplexToRealFFTW<T> C2R;
};

	namespace fft {

		typedef detail::FFTPolicy<double, detail::FFTW> fftw_f64_t;
		typedef detail::FFTPolicy< float, detail::FFTW> fftw_f32_t;

		static const fftw_f32_t fftw_f32=fftw_f32_t();

		static const fftw_f64_t fftw_f64=fftw_f64_t();


	}  // namespace fft

}  // namespace hydra

#endif /* FFTW_H_ */
