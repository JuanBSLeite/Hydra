/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2016 - 2017 Antonio Augusto Alves Junior
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
 * DeltaDMassBackground.h
 *
 *  Created on: Jul 31, 2018
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef DELTADMASSBACKGROUND_H_
#define DELTADMASSBACKGROUND_H_


namespace hydra {


template<unsigned int ArgIndex=0>
class DeltaDMassBackground: public BaseFunctor<DeltaDMassBackground<ArgIndex>, double, 4>
{
	using BaseFunctor<DeltaDMassBackground<ArgIndex>, double, 4>::_par;

public:

	DeltaDMassBackground() = delete;

	DeltaDMassBackground(Parameter const& threshold, Parameter const& A, Parameter const& B, Parameter const& C):
		BaseFunctor<DeltaDMassBackground<ArgIndex>, double, 4>({threshold, A, B, C})
		{}

	__hydra_host__ __hydra_device__
	DeltaDMassBackground(DeltaDMassBackground<ArgIndex> const& other ):
	BaseFunctor<DeltaDMassBackground<ArgIndex>, double,4>(other)
	{}

	__hydra_host__ __hydra_device__
	DeltaDMassBackground<ArgIndex>&
	operator=(DeltaDMassBackground<ArgIndex> const& other ){
		if(this==&other) return  *this;
		BaseFunctor<DeltaDMassBackground<ArgIndex>,double, 4>::operator=(other);
		return  *this;
	}

	template<typename T>
	__hydra_host__ __hydra_device__
	inline double Evaluate(unsigned int, T* x)  const	{

		double arg   = (x[ArgIndex] - _par[0]);
		double ratio = (arg/_par[0]);
		double val   = arg>0 ? (1- ::exp(-arg/_par[3]))*::pow(ratio, _par[1]) + _par[2]*(ratio-1) : 0;

		return  CHECK_VALUE( (val>0 ? val : 0), "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f ", _par[0], _par[1], _par[2], _par[3]);

	}

	template<typename T>
	__hydra_host__ __hydra_device__
	inline double Evaluate(T x)  const {

		double arg   = get<ArgIndex>(x) - _par[0];
		double ratio = (arg/_par[0]);
		double val   = arg>0 ? (1- ::exp(-arg/_par[3]))*::pow(ratio, _par[1]) + _par[2]*(ratio-1) : 0;

		return  CHECK_VALUE( (val>0 ? val : 0), "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f ", _par[0], _par[1], _par[2], _par[3]);

	}

};


}// namespace hydra

#endif /* DELTADMASSBACKGROUND_H_ */
