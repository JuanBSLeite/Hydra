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
 * Flatte.h
 *
 *  Created on: 18/12/2017
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef Flatte_H_
#define Flatte_H_


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

#include <hydra/functions/FlatteLineShape.h>
#include <hydra/functions/CosHelicityAngle.h>
#include <hydra/functions/ZemachFunctions.h>

/**
 * \ingroup common_functions
 *
 * \class Flatte
 *
 * Implementation describing the ARGUS background shape.
 *
 * \tparam ArgIndex : index of the argument when evaluating on multidimensional data. Default is 0.
 */

namespace hydra{

template<unsigned int CHANNEL, hydra::Wave L>
class Flatte: public hydra::BaseFunctor<Flatte<CHANNEL,L>, hydra::complex<double>, 5>
{
	using hydra::BaseFunctor<Flatte<CHANNEL,L>, hydra::complex<double>, 5>::_par;

	constexpr static unsigned int _I1 = CHANNEL-1;
	constexpr static unsigned int _I2 = (CHANNEL!=3)*CHANNEL;
	constexpr static unsigned int _I3 = 3-( (CHANNEL-1) + (CHANNEL!=3)*CHANNEL );


public:

	Flatte() = delete;

	Flatte(hydra::Parameter const& c_re, hydra::Parameter const& c_im,
			hydra::Parameter const& mass, hydra::Parameter const& GPP,hydra::Parameter const& GKK,
			double mother_mass,	double daugther1_mass,
			double daugther2_mass, double daugther3_mass,
			double radi):
			hydra::BaseFunctor<Flatte<CHANNEL,L>, hydra::complex<double>, 5>{c_re, c_im, mass, GPP, GKK},
			fLineShape(mass, GPP, GKK, mother_mass, daugther1_mass, daugther2_mass, daugther3_mass, radi)
	{}


    __hydra_dual__
	Flatte( Flatte< CHANNEL,L> const& other):
	hydra::BaseFunctor<Flatte<CHANNEL ,L>, hydra::complex<double>, 5>(other),
	fLineShape(other.GetLineShape())
	{}

    __hydra_dual__  inline
	Flatte< CHANNEL ,L>&
	operator=( Flatte< CHANNEL ,L> const& other)
	{
		if(this==&other) return *this;

		hydra::BaseFunctor<Flatte<CHANNEL ,L>, hydra::complex<double>, 5>::operator=(other);
		fLineShape=other.GetLineShape();

		return *this;
	}

    __hydra_dual__  inline
	hydra::FlatteLineShape<L> const& GetLineShape() const {	return fLineShape; }

    __hydra_dual__  inline
	hydra::complex<double> Evaluate(unsigned int n, hydra::Vector4R* p)  const {


		hydra::Vector4R p1 = p[_I1];
		hydra::Vector4R p2 = p[_I2];
		hydra::Vector4R p3 = p[_I3];


		fLineShape.SetParameter(0, _par[2]);
		fLineShape.SetParameter(1, _par[3]);
		fLineShape.SetParameter(2, _par[4]);

		double theta = fCosDecayAngle( (p1+p2+p3), (p1+p2), p1 );
		double angular = fAngularDist(theta);
		auto r = hydra::complex<double>(_par[0], _par[1])*fLineShape((p1+p2).mass());//*angular;

		return r;

	}

private:

	mutable hydra::FlatteLineShape<L> fLineShape;
	hydra::CosHelicityAngle fCosDecayAngle;
	hydra::ZemachFunction<L> fAngularDist;


};
}

#endif /* Flatte_H_ */
