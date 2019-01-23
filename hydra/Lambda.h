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
    typedef   Closure super_type;

    Lambda()=delete;

    Lambda(Closure const& closure,  typename std::enable_if<NPARAM==0, bool>::type hiden_flag=true):
    super_type(closure),
    detail::Parameters<NPARAM>(),
    fCacheIndex(-1),
    fCached(0),
    fNorm(1.0)
    {}

    Lambda(Closure const& closure, std::array<Parameter,NPARAM> const& init_parameters,
    		 typename std::enable_if<NPARAM==0, bool>::type  hiden_flag=true):
    super_type(closure),
    detail::Parameters<NPARAM>( init_parameters ),
    fCacheIndex(-1),
    fCached(0),
    fNorm(1.0)
    { }


	__hydra_host__ __hydra_device__
	Lambda(Lambda<Closure,ReturnType, NPARAM> const& other):
	super_type( other),
	detail::Parameters<NPARAM>( other),
	fCacheIndex( other.GetCacheIndex() ),
	fCached( other.IsCached() ),
	fNorm(other.GetNorm())
	{ }

	__hydra_host__ __hydra_device__
	inline Lambda<Closure, ReturnType, NPARAM>&
	operator=(Lambda<Closure, ReturnType, NPARAM> const & other )
	{
		if(this == &other) return *this;

		detail::Parameters<NPARAM>::operator=( other );
		this->fCacheIndex     = other.GetCacheIndex();
		this->fCached         = other.IsCached();
		this->fNorm = other.GetNorm();


		return *this;
	}

	__hydra_host__ __hydra_device__
	inline int GetCacheIndex() const { return this->fCacheIndex; }

	__hydra_host__ __hydra_device__
	inline void SetCacheIndex(int index) {fCacheIndex = index;}

	__hydra_host__ __hydra_device__
	inline bool IsCached() const
	{ return this->fCached;}

	__hydra_host__ __hydra_device__
	inline void SetCached(bool cached=true)
	{ fCached = cached; }


	void PrintRegisteredParameters()
	{

		HYDRA_CALLER ;
		HYDRA_MSG <<HYDRA_ENDL;
		HYDRA_MSG << "Registered parameters begin:" << HYDRA_ENDL;
		this->PrintParameters();

		HYDRA_MSG <<"Normalization " << fNorm << HYDRA_ENDL;
		HYDRA_MSG <<"Registered parameters end." << HYDRA_ENDL;
		HYDRA_MSG <<HYDRA_ENDL;
		return;
	}


	__hydra_host__ __hydra_device__
	inline GReal_t GetNorm() const {
		return fNorm;
	}

	__hydra_host__ __hydra_device__
	inline void SetNorm(GReal_t norm) {
		fNorm = norm;
	}


	template<typename T, size_t M = NPARAM>
	__hydra_host__ __hydra_device__
	inline typename HYDRA_EXTERNAL_NS::thrust::detail::enable_if<
	        detail::is_tuple_type<T>::value && M == 0,
	        return_type>::type
	operator()( T&&  x )  const
	{
		return  super_type::operator()( std::forward<T>(x));
	}

	template<typename T, size_t M = NPARAM>
	__hydra_host__ __hydra_device__
	inline typename HYDRA_EXTERNAL_NS::thrust::detail::enable_if<
	detail::is_tuple_type<T>::value && M != 0,
	return_type>::type
	operator()( T&&  x )  const
	{
		return  super_type::operator()(this->GetNumberOfParameters(),
				this->GetParameters(), std::forward<T>(x));
	}

	template<typename T, size_t M = NPARAM>
	__hydra_host__ __hydra_device__
	inline typename HYDRA_EXTERNAL_NS::thrust::detail::enable_if<
	(!detail::is_tuple_type<T>::value) && M == 0,
	return_type>::type
	operator()( T&&  x )  const
	{
		return  super_type::operator()(  HYDRA_EXTERNAL_NS::thrust::tie(std::forward<T>(x)));
	}

	template<typename T, size_t M = NPARAM>
	__hydra_host__ __hydra_device__
	inline typename HYDRA_EXTERNAL_NS::thrust::detail::enable_if<
	(!detail::is_tuple_type<T>::value) && M != 0,
	return_type>::type
	operator()( T&&  x )  const
	{
		return  super_type::operator()(this->GetNumberOfParameters(),
				this->GetParameters(), HYDRA_EXTERNAL_NS::thrust::tie(std::forward<T>(x)));
	}



private:

   int fCacheIndex;
	bool fCached;
   GReal_t fNorm;


};

template<typename ReturnType, typename Closure, typename ...T>
Lambda<Closure, ReturnType, sizeof...(T)>
wrap_lambda(Closure const& f,  T const& ...pars) {

	return Lambda< Closure, ReturnType,sizeof...(T)>(f, std::array<Parameter, sizeof...(T)>{ pars...});
}


template<typename ReturnType, typename Closure>
Lambda<Closure, ReturnType, 0>
wrap_lambda(Closure const& f) {

	return Lambda< Closure, ReturnType,0>(f);
}


}// namespace hydra


#endif /* LAMBDA_H_ */
