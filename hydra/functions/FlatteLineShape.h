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
 * FlatteLineShape.h
 *
 *  Created on: 26/12/2017
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef FlatteLineShape_H_
#define FlatteLineShape_H_

#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <hydra/Complex.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/BlattWeisskopfFunctions.h>


#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <cmath>


namespace hydra {


/**
 * \ingroup common_functions
 *
 * \class FlatteLineShape
 *
 * Breit-Wigner line shape for 3 body resonant decays \f$ A -> r c , r-> a b\f$ ,
 * where A is a "long-lived" particle and \f$ a, b\f$ and \f$c\f$ are the final states.
 * The lineshape is defined by the expression:
 *
 * \f[
 *  R(m_{a,b}|m_0,\Lambda_0) = B'_{L_A}(d, p_0, p)(\frac{p}{m_A})^{L_A} \times \\
 *  		BW(m_{a,b}|m_0,\Lambda_0) \times B'_{L_r}(d, q_0, q)(\frac{q}{q_r})^{L_r}
 * \f]
 *
 * where Breit-Wigner amplitude is given by:
 *
 *\f[ BW(m_{ab}|m_0,\Lambda_0)= \frac{1}{m_0^2 - m_{ab}^2 - im_0\Lambda(m_{ab})} \f]
 *
 *and
 *\f[  \Lambda(m_{ab}) = \Lambda_0(\frac{q}{q_0})^{2L_{r}+1}\frac{m_0}{m}B'_{L_r}(d, q_0, q)\f]
 *
 *@tparam ResonanceWave hydra::Wave resonance decay vertex wave
 *@tparam MotherWave hydra::Wave mother particle decay vertex wave
 */
template<Wave ResonanceWave, Wave MotherWave=SWave, unsigned int ArgIndex=0>
class FlatteLineShape : public BaseFunctor<FlatteLineShape< ResonanceWave,MotherWave,ArgIndex>, hydra::complex<double>, 3>
{
	using BaseFunctor<FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>, hydra::complex<double>, 3>::_par;

public:

	FlatteLineShape()=delete;

	/**
	 *
	 * @param mass resonance mass.
	 * @param width resonance width.
	 * @param mother_mass resonance mother mass.
	 * @param daugther1_mass resonance daughter particle 1 mass
	 * @param daugther2_mass resonance daughter particle 2 mass
	 * @param daugther3_mass daughter particle 2 mass
	 * @param radi decay vertex radio.
	 */
	FlatteLineShape(Parameter const& mass, Parameter const& GPP, Parameter const& GKK,
			double mother_mass,
			double daugther1_mass, double daugther2_mass, double bachelor_mass,
			double radi):
		BaseFunctor<FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>, hydra::complex<double>, 3>{mass,GPP,GKK},
		fDaughter1Mass(daugther1_mass),
		fDaughter2Mass(daugther2_mass),
		fBachelorMass(bachelor_mass),
		fMotherMass(mother_mass),
		fRadi(radi)
	{}

	__hydra_host__  __hydra_device__
	FlatteLineShape(FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>  const& other):
		BaseFunctor<FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>, hydra::complex<double>, 3>(other),
		fDaughter1Mass(other.GetDaughter1Mass()),
		fDaughter2Mass(other.GetDaughter2Mass()),
		fBachelorMass(other.GetBachelorMass()),
		fMotherMass(other.GetMotherMass()),
		fRadi(other.GetRadi())
		{}

	__hydra_host__  __hydra_device__
	FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>&
	operator=(FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>  const& other)
	{
		if(this==&other) return  *this;

		BaseFunctor<FlatteLineShape<ResonanceWave,MotherWave,ArgIndex>,
			hydra::complex<double>, 3>::operator=(other);

		fDaughter1Mass= other.GetDaughter1Mass();
		fDaughter2Mass= other.GetDaughter2Mass();
		fBachelorMass= other.GetBachelorMass();
		fMotherMass= other.GetMotherMass();
		fRadi= other.GetRadi();

		 return  *this;
	}

	__hydra_host__  __hydra_device__ inline
	double GetDaughter1Mass() const {
		return fDaughter1Mass;
	}

	__hydra_host__  __hydra_device__ inline
	void SetDaughter1Mass(double daughter1Mass) {
		fDaughter1Mass = daughter1Mass;
	}

	__hydra_host__  __hydra_device__ inline
	double GetDaughter2Mass() const {
		return fDaughter2Mass;
	}

	__hydra_host__  __hydra_device__ inline
	void SetDaughter2Mass(double daughter2Mass) {
		fDaughter2Mass = daughter2Mass;
	}

	__hydra_host__  __hydra_device__ inline
	double GetBachelorMass() const {
		return fBachelorMass;
	}

	__hydra_host__  __hydra_device__ inline
	void SetBachelorMass(double daughter3Mass) {
		fBachelorMass = daughter3Mass;
	}

	__hydra_host__  __hydra_device__ inline
	double GetMotherMass() const {
		return fMotherMass;
	}

	__hydra_host__  __hydra_device__ inline
	void SetMotherMass(double motherMass) {
		fMotherMass = motherMass;
	}

	__hydra_host__  __hydra_device__ inline
	double GetRadi() const {
		return fRadi;
	}

	__hydra_host__  __hydra_device__ inline
	void SetRadi(double radi) {
		fRadi = radi;
	}

	template<typename T>
	__hydra_host__ __hydra_device__ inline
	hydra::complex<double> Evaluate(unsigned int, const T*x)  const	{

		const double m = x[ArgIndex] ;

		const double resonance_mass  = _par[0];
		const double resonance_GPP = _par[1];
		const double resonance_GKK = _par[2];

		return  m > (fDaughter1Mass+fDaughter2Mass) && m<(fMotherMass-fBachelorMass) ?
				LineShape(m, resonance_mass, resonance_GPP, resonance_GKK): hydra::complex<double>(0.0, 0.0) ;

	}

	template<typename T>
	__hydra_host__ __hydra_device__ inline
	hydra::complex<double> Evaluate(T& x)  const {

		double m =  get<ArgIndex>(x);

		const double resonance_mass  = _par[0];
		const double resonance_GPP = _par[1];
		const double resonance_GKK = _par[2];

		return  m > (fDaughter1Mass+fDaughter2Mass) && m<(fMotherMass-fBachelorMass) ?
				LineShape(m, resonance_mass, resonance_GPP, resonance_GKK): hydra::complex<double>(0.0, 0.0) ;
	}



private:

	 __hydra_host__ __hydra_device__   inline
	 hydra::complex<double> LineShape(const double s, const double resonance_mass, const double resonance_GPP, const double resonance_GKK ) const {

    double resmass            = resonance_mass;
    double g1                 = resonance_GPP;
    double g2                 = resonance_GKK * g1;
    
    double pipmass = 0.13957018;
    double pi0mass = 0.1349766;
    double kpmass  = 0.493677;
    double k0mass  = 0.497614;

    double twopimasssq  = 4 * pipmass * pipmass;
    double twopi0masssq = 4 * pi0mass * pi0mass;
    double twokmasssq   = 4 * kpmass * kpmass;
    double twok0masssq  = 4 * k0mass * k0mass;
   
        double rhopipi_real = 0, rhopipi_imag = 0;
        double rhokk_real = 0, rhokk_imag = 0;

        if(s >= twopimasssq)
            rhopipi_real += (2. / 3) * ::sqrt(1 - twopimasssq / s); // Above pi+pi- threshold
        else
            rhopipi_imag += (2. / 3) * ::sqrt(-1 + twopimasssq / s);
        if(s >= twopi0masssq)
            rhopipi_real += (1. / 3) * ::sqrt(1 - twopi0masssq / s); // Above pi0pi0 threshold
        else
            rhopipi_imag += (1. / 3) * ::sqrt(-1 + twopi0masssq / s);
        if(s >= twokmasssq)
            rhokk_real += 0.5 * ::sqrt(1 - twokmasssq / s); // Above K+K- threshold
        else
            rhokk_imag += 0.5 * ::sqrt(-1 + twokmasssq / s);
        if(s >= twok0masssq)
            rhokk_real += 0.5 * ::sqrt(1 - twok0masssq / s); // Above K0K0 threshold
        else
            rhokk_imag += 0.5 * ::sqrt(-1 + twok0masssq / s);
        double A = (resmass * resmass - s) + resmass * (rhopipi_imag * g1 + rhokk_imag * g2);
        double B = resmass * (rhopipi_real * g1 + rhokk_real * g2);
        double C = 1.0 / (A * A + B * B);
        hydra::complex<double> retur(A * C, B * C);
        
    

    
    return retur;

	 }

	double fDaughter1Mass;
	double fDaughter2Mass;
	double fBachelorMass;
	double fMotherMass;
	double fRadi;

};



}  // namespace hydra


#endif /* FlatteLineShape_H_ */
