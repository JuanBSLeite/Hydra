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
 * GounarisSakurai.h
 *
 *  Created on: 18/12/2017
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef GounarisSakurai_H_
#define GounarisSakurai_H_


#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>

#include <hydra/functions/GounarisSakuraiLineShape.h>
#include <hydra/functions/CosHelicityAngle.h>
#include <hydra/functions/ZemachFunctions.h>

#include <tuple>
#include <limits>
#include <stdexcept>
#include <cassert>
#include <utility>

/**
 * \ingroup common_functions
 *
 * \class GounarisSakurai
 *
 * Implementation describing the ARGUS background shape.
 *
 * \tparam ArgIndex : index of the argument when evaluating on multidimensional data. Default is 0.
 */

namespace hydra{

template<unsigned int CHANNEL, hydra::Wave L>
class GounarisSakurai: public hydra::BaseFunctor<GounarisSakurai<CHANNEL,L>, hydra::complex<double>, 4>
{
	using hydra::BaseFunctor<GounarisSakurai<CHANNEL,L>, hydra::complex<double>, 4>::_par;

	constexpr static unsigned int _I1 = CHANNEL-1;
	constexpr static unsigned int _I2 = (CHANNEL!=3)*CHANNEL;
	constexpr static unsigned int _I3 = 3-( (CHANNEL-1) + (CHANNEL!=3)*CHANNEL );


public:

	GounarisSakurai() = delete;

	GounarisSakurai(hydra::Parameter const& c_re, hydra::Parameter const& c_im,
			hydra::Parameter const& mass, hydra::Parameter const& width,
			double mother_mass,	double daugther1_mass,
			double daugther2_mass, double daugther3_mass,
			double radi):
			hydra::BaseFunctor<GounarisSakurai<CHANNEL,L>, hydra::complex<double>, 4>{c_re, c_im, mass,width},
			fLineShape(mass,width, mother_mass, daugther1_mass, daugther2_mass, daugther3_mass, radi)
	{}


    __hydra_dual__
	GounarisSakurai( GounarisSakurai< CHANNEL,L> const& other):
	hydra::BaseFunctor<GounarisSakurai<CHANNEL ,L>, hydra::complex<double>, 4>(other),
	fLineShape(other.GetLineShape())
	{}

    __hydra_dual__  inline
	GounarisSakurai< CHANNEL ,L>&
	operator=( GounarisSakurai< CHANNEL ,L> const& other)
	{
		if(this==&other) return *this;

		hydra::BaseFunctor<GounarisSakurai<CHANNEL ,L>, hydra::complex<double>, 4>::operator=(other);
		fLineShape=other.GetLineShape();

		return *this;
	}

    __hydra_dual__  inline
	hydra::GounarisSakuraiLineShape<L> const& GetLineShape() const {	return fLineShape; }

    __hydra_dual__  inline
	hydra::complex<double> Evaluate(unsigned int n, hydra::Vector4R* p)  const {


		hydra::Vector4R p1 = p[_I1];
		hydra::Vector4R p2 = p[_I2];
		hydra::Vector4R p3 = p[_I3];


		fLineShape.SetParameter(0, _par[2]);
		fLineShape.SetParameter(1, _par[3]);
		

		double theta = fCosDecayAngle( (p1+p2+p3), (p1+p2), p1 );
		double angular = fAngularDist(theta);
		auto r = hydra::complex<double>(_par[0], _par[1])*fLineShape((p1+p2).mass())*angular;

		return r;

	}

private:

	mutable hydra::GounarisSakuraiLineShape<L> fLineShape;
	hydra::CosHelicityAngle fCosDecayAngle;
	hydra::ZemachFunction<L> fAngularDist;


};
}

#endif /* GounarisSakurai_H_ */
