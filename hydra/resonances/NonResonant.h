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
 * NR.h
 *
 *  Created on: 18/12/2017
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef NR_H_
#define NR_H_


#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <cassert>
#include <utility>



/**
 * \ingroup common_functions
 *
 * \class NR
 *
 * Implementation describing the ARGUS background shape.
 *
 * \tparam ArgIndex : index of the argument when evaluating on multidimensional data. Default is 0.
 */

namespace hydra{

class NonResonant: public hydra::BaseFunctor<NonResonant, hydra::complex<double>, 2>
{

	using hydra::BaseFunctor<NonResonant, hydra::complex<double>, 2>::_par;

public:

	NonResonant() = delete;

	NonResonant(hydra::Parameter const& c_re, hydra::Parameter const& c_im):
			hydra::BaseFunctor<NonResonant, hydra::complex<double>, 2>{c_re, c_im}
	{}


	 __hydra_dual__
	NonResonant( NonResonant const& other):
	hydra::BaseFunctor<NonResonant, hydra::complex<double>, 2>(other)
	{}

	 __hydra_dual__
	NonResonant& operator=( NonResonant const& other)
	{
		if(this==&other) return *this;

		hydra::BaseFunctor<NonResonant, hydra::complex<double>, 2>::operator=(other);

		return *this;
	}

	 __hydra_dual__  inline
	hydra::complex<double> Evaluate(unsigned int n, hydra::Vector4R* p)  const {

		return hydra::complex<double>(_par[0], _par[1]);
	}

};

}

#endif /* NR_H_ */
