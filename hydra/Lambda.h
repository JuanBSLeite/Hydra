/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2016 - 2019 Antonio Augusto Alves Junior
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
 * Lambda.h
 *
 *  Created on: Jan 22, 2019
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef LAMBDA_H_
#define LAMBDA_H_


#include <hydra/detail/Config.h>
#include <hydra/Types.h>
#include <hydra/detail/Print.h>
#include <hydra/Integrator.h>
#include <hydra/Parameter.h>
#include <hydra/detail/utility/Utility_Tuple.h>
#include <hydra/detail/FunctorTraits.h>
#include <hydra/detail/Parameters.h>


#include <hydra/detail/external/thrust/iterator/detail/tuple_of_iterator_references.h>
#include <hydra/detail/external/thrust/iterator/zip_iterator.h>
#include <hydra/detail/external/thrust/tuple.h>
#include <hydra/detail/external/thrust/detail/type_traits.h>
#include <array>
#include <initializer_list>
#include <memory>

namespace hydra {

template<typename Closure, typename ReturnType, size_t NPARAM>
class  Lambda: public Closure, public detail::Parameters<NPARAM>
{

public:

	//tag
    typedef   void hydra_functor_tag;
	typedef   ReturnType return_type;
	typedef   std::true_type is_functor;




}; // namespace hydra


#endif /* LAMBDA_H_ */
