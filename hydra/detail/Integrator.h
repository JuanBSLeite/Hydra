/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2016 Antonio Augusto Alves Junior
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
 * Integrator.h
 *
 *  Created on: 31/08/2016
 *      Author: Antonio Augusto Alves Junior
 */

/**
 * \file
 * \ingroup numerical_integration
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <hydra/detail/Config.h>
#include <hydra/Types.h>
#include <thrust/pair.h>

namespace hydra {

template<typename ALGORITHM, size_t N>
struct Integrator{

	//tag
	typedef void hydra_integrator_tag;

	static const size_t dimension=N;

	Integrator()
	{};

	template<typename FUNCTOR>
	__host__ inline
	thrust::pair<GReal_t, GReal_t> GetIntegral( FUNCTOR const& fFunctor ){

	   static_cast<ALGORITHM*>(this)->Integrate(fFunctor );

	   return thrust::make_pair( static_cast<ALGORITHM*>(this)->GetResult(),
	    		static_cast<ALGORITHM*>(this)->GetAbsError() );
	}


	__host__ inline
	thrust::pair<const GReal_t*, const GReal_t*> GetLimits( )
	{
		return thrust::make_pair( static_cast<ALGORITHM*>(this)->GetLowerLimit(),
				static_cast<ALGORITHM*>(this)->GetUpperLimit() );

	}




};


}  // namespace hydra



#endif /* INTEGRATOR_H_ */
