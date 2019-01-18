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
 * FunctionWrapper_cpp14.inl
 *
 *  Created on: 06/12/2018
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef FUNCTIONWRAPPER_CPP14_INL_
#define FUNCTIONWRAPPER_CPP14_INL_


#include <hydra/detail/Config.h>
#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include <hydra/detail/utility/Generic.h>
#include <type_traits>
#include <hydra/Function.h>


namespace hydra {

namespace function_wrapper
{

struct SingleArg{};

}


template<typename Lambda, typename ReturnType, size_t N>
class LambdaWrapper:public BaseFunctor<LambdaWrapper<Lambda, ReturnType, N>, ReturnType,N >
{

public:

	LambdaWrapper()=delete;

	/**
	 * Constructor for parametrized lambdas
	 * @param lambda
	 * @param parameters
	 */
	LambdaWrapper(L const& lambda,	std::array<Parameter, N> const& parameters):
				BaseFunctor<LambdaWrapper<Lambda, ReturnType, N>, ReturnType,N >(parameters),
		fLambda(lambda)
	{}

	/**
	 * Copy constructor
	 */
	__hydra_host__ __hydra_device__
	inline	LambdaWrapper(LambdaWrapper<Lambda, ReturnType, N> const& other ):
	BaseFunctor<LambdaWrapper<Lambda, ReturnType, N>, ReturnType,N>(other),
	fLambda( other.GetLambda())
	{	}

	/**
	 * Assignment operator
	 */
	__hydra_host__ __hydra_device__
	inline LambdaWrapper<Lambda, ReturnType, N>
	operator=(LambdaWrapper<Lambda, ReturnType, N> const& other )
	{
		if(this==&other) return *this;

		BaseFunctor<LambdaWrapper<Lambda, ReturnType, N>, ReturnType,N>::operator=(other);
		fLambda=other.GetLambda();

		return *this;
	}


	/**
	 * Get the underlying lambda
	 */
	__hydra_host__ __hydra_device__
	inline const Lambda& GetLambda() const {return fLambda; }


	template< typename T, size_t M=N >
	__hydra_host__ __hydra_device__
	inline typename std::enable_if< (M>0), ReturnType >::type
	Evaluate(T&& a)   const {

		return fLambda(this->GetNumberOfParameters(), this->GetParameters(), std::forward<T>(a) );
	}



	template< typename ...T, size_t M=N >
	__hydra_host__ __hydra_device__
	inline typename std::enable_if< (M==0), ReturnType >::type
	Evaluate(T&& a)   const {

		return fLambda( std::forward<T>(a) );
	}


private:
	L fLambda;
};


/**
 * @ingroup functor
 * @brief Function template for wrap a C++14 lambda into a hydra lambda.
 * @param f single argument C++14 lambda *
 * @param pars parameters.
 * @return LambdaWrapper object.
 */
template<typename L, typename ...T>
auto wrap_lambda(L const& f,  T const& ...pars)
-> LambdaWrapper<L, typename std::result_of<L(hydra::function_wrapper::SingleArg)>::type, sizeof...(T)>
{
	typedef typename std::result_of<L(SingleArg,
			std::array<Parameter, sizeof...(T)>)>::type result_type;

	std::array<Parameter, sizeof...(T)> parameters{ pars...};


	return LambdaWrapper<L, result_type,0>(f, parameters);
}

/**
 * @ingroup functor
 * @brief Function template for wrap a C++14 lambda into a hydra lambda.
 * @param f single argument C++14 lambda
 * @return LambdaWrapper object
 */
template<typename L>
auto wrap_lambda(L const& f)
-> LambdaWrapper<L, typename std::result_of<L(hydra::function_wrapper::SingleArg)>::type,0>
{
	typedef hydra::function_wrapper::SingleArg SingleArg;

	typedef typename std::result_of<L(SingleArg)>::type result_type;

	return LambdaWrapper<L, result_type,0>(f);
}



}



#endif /* FUNCTIONWRAPPER_CPP14_INL_ */
