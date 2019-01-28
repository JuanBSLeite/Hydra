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
 * Utils.h
 *
 *  Created on: 18/12/2017
 *      Author: Antonio Augusto Alves Junior
 */

#ifndef Utils_H_
#define Utils_H_


#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <hydra/PhaseSpace.h>
#include <hydra/PhaseSpaceIntegrator.h>
#include <hydra/Decays.h>
#include <hydra/Complex.h>
#include <hydra/DenseHistogram.h>
#include <hydra/SparseHistogram.h>
#include <hydra/detail/Compose.h>


#include <tuple>
#include <limits>
#include <stdexcept>
#include <cassert>
#include <utility>

#ifdef _ROOT_AVAILABLE_
#include <TH3D.h>
#endif //_ROOT_AVAILABLE_

#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))

__hydra_host__ __hydra_device__
inline double dampingFactorSquare(const double &cmmom, const int &spin, const double &mRadius) {
    double square = mRadius * mRadius * cmmom * cmmom;
    double dfsq   = 1 + square; // This accounts for spin 1
    // if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.
    double dfsqres = dfsq + 8 + 2 * square + square * square;

    // Spin 3 and up not accounted for.
    // return dfsq;
    return (spin == 2) ? dfsqres : dfsq;
}


template<typename Backend, typename Model, typename Container >
size_t generate_dataset(Backend const& system, Model const& model, std::array<double, 3> const& masses, Container& decays, size_t nevents, size_t bunch_size)
{
	const double D_MASS         = masses[0];// D+ mass
	const double K_MASS         = masses[1];// K+ mass
	const double PI_MASS        = masses[2];// pi mass

	//generator
	hydra::Vector4R D(D_MASS, 0.0, 0.0, 0.0);

	// Create PhaseSpace object for B0-> K pi pi
	hydra::PhaseSpace<3> phsp{K_MASS, PI_MASS, PI_MASS};

	//allocate memory to hold the final states particles
	hydra::Decays<3, Backend > _data(bunch_size);


	do {
		phsp.SetSeed(753156);

		//generate the final state particles
		phsp.Generate(D, _data.begin(), _data.end());

		auto last = _data.Unweight(model, 1.0);

		decays.insert(decays.size()==0? decays.begin():decays.end(),
				_data.begin(), _data.begin()+last );

	} while(decays.size()<nevents );

	decays.erase(decays.begin()+nevents, decays.end());

	return decays.size();

}


template<typename Amplitude, typename Model>
double fit_fraction( Amplitude const& amp, Model const& model, std::array<double, 3> const& masses, size_t nentries)
{
	const double D_MASS         = masses[0];// D+ mass
	const double K_MASS         = masses[1];// K+ mass
	const double PI_MASS        = masses[2];// pi mass

	//--------------------
	//generator
	hydra::Vector4R D(D_MASS, 0.0, 0.0, 0.0);

	// Create PhaseSpace object for B0-> K pi pi
	hydra::PhaseSpace<3> phsp{K_MASS, PI_MASS, PI_MASS};

	//norm lambda
	auto Norm = hydra::wrap_lambda( [] __hydra_dual__ (unsigned int n, hydra::complex<double>* x){

		return hydra::norm(x[0]);
	});

	//functor
	auto functor = hydra::compose(Norm, amp);


	auto amp_int   = phsp.AverageOn(hydra::device::sys, D, functor, nentries);
	auto model_int = phsp.AverageOn(hydra::device::sys, D, model,   nentries);


	return amp_int.first/model_int.first;

}


template<typename Amplitude>
TH3D histogram_component( Amplitude const& amp, std::array<double, 3> const& masses, const char* name, size_t nentries)
{
	const double D_MASS         = masses[0];// D+ mass
	const double K_MASS         = masses[1];// K+ mass
	const double PI_MASS        = masses[2];// pi mass

	TH3D Component(name,
			";"
			"M^{2}(K^{-} #pi^{+}) [GeV^{2}/c^{4}];"
			"M^{2}(K^{-} #pi^{+}) [GeV^{2}/c^{4}];"
			"M^{2}(#pi^{+} #pi^{+}) [GeV^{2}/c^{4}]",
			100, POW2(K_MASS  + PI_MASS), POW2(D_MASS - PI_MASS),
			100, POW2(K_MASS  + PI_MASS), POW2(D_MASS - PI_MASS),
			100, POW2(PI_MASS + PI_MASS), POW2(D_MASS -  K_MASS));

	//--------------------
	//generator
	hydra::Vector4R D(D_MASS, 0.0, 0.0, 0.0);
	// Create PhaseSpace object for B0-> K pi pi
	hydra::PhaseSpace<3> phsp{K_MASS, PI_MASS, PI_MASS};

	// functor to calculate the 2-body masses
	auto dalitz_calculator = hydra::wrap_lambda(
			[] __hydra_dual__ (unsigned int n, hydra::Vector4R* p ){

		double   M2_12 = (p[0]+p[1]).mass2();
		double   M2_13 = (p[0]+p[2]).mass2();
		double   M2_23 = (p[1]+p[2]).mass2();

		return hydra::make_tuple(M2_12, M2_13, M2_23);
	});

	//norm lambda
	auto Norm = hydra::wrap_lambda(
			[] __hydra_dual__ (unsigned int n, hydra::complex<double>* x){

		return hydra::norm(x[0]);
	});

	//functor
	auto functor = hydra::compose(Norm, amp);

	hydra::Decays<3, hydra::device::sys_t > events(nentries);

	phsp.Generate(D, events.begin(), events.end());

	events.Reweight(functor);

	auto particles        = events.GetUnweightedDecays();
	auto dalitz_variables = hydra::make_range( particles.begin(), particles.end(), dalitz_calculator);
	auto dalitz_weights   = events.GetWeights();

	//model dalitz histogram
	hydra::SparseHistogram<double, 3,  hydra::device::sys_t> Hist_Component{
		{100,100,100},
		{POW2(K_MASS + PI_MASS), POW2(K_MASS + PI_MASS),  POW2(PI_MASS + PI_MASS)},
		{POW2(D_MASS - PI_MASS), POW2(D_MASS - PI_MASS), POW2(D_MASS - K_MASS)}
	};

	Hist_Component.Fill( dalitz_variables.begin(),
					dalitz_variables.end(), dalitz_weights.begin()  );

	for(auto entry : Hist_Component){

		size_t bin     = hydra::get<0>(entry);
		double content = hydra::get<1>(entry);
		unsigned int bins[3];
		Hist_Component.GetIndexes(bin, bins);
		Component.SetBinContent(bins[0]+1, bins[1]+1, bins[2]+1, content);

	}

	return Component;

}

__hydra_host__ __hydra_device__ inline
double twoBodyCMmom(double rMassSq, double d1m, double d2m) {
    // For A -> B + C, calculate momentum of B and C in rest frame of A.
    // PDG 38.16.

    double kin1 = 1 - POW2(d1m + d2m) / rMassSq;

    kin1 = kin1 >= 0 ? sqrt(kin1) : 1;

    double kin2 = 1 - POW2(d1m - d2m) / rMassSq;
    kin2        = kin2 >= 0 ? sqrt(kin2) : 1;

    return 0.5 * sqrt(rMassSq) * kin1 * kin2;
}

template <typename T,int N>
std::array<T,N> flatten(const std::array<hydra::complex<T>,N/2> &input) {
std::array<T,N> output;
for(size_t i = 0; i <= (N-2);i++){;

    output[i]   = input[i].real();
    output[i+1] = input[i].imag();

}

return output;
}


template <typename T, int N>
std::array<hydra::complex<T>,N> complex_derivative(const std::array<T,N> &x, const std::array<hydra::complex<T>,N> &y) {
	if(x.size() != y.size())
		 std::cout << "x and y must have the same diminsions!" << std::endl;

	int i, k;
	unsigned int n = N;
	std::array<hydra::complex<T>,N> u;
	std::array<hydra::complex<T>,N> y2;

	double sig, p, qn, un;
	hydra::complex<T> yp1 = 2. * (y[1] - y[0]) / (x[1] - x[0]);
	hydra::complex<T> ypn = 2. * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/
	if(yp1.real() > 0.99e30) {
		y2[0].real(0.);
		u[0].real(0.);
	} else {
		y2[0].real(-0.5);
		u[0].real(3.0 / (x[1] - x[0]) * ((y[1].real() - y[0].real()) / (x[1] - x[0]) - yp1.real()));
	}
	if(yp1.imag() > 0.99e30) {
		y2[0].imag(0.);
		u[0].imag(0.);
	} else {
		y2[0].imag(-0.5);
		u[0].imag(3.0 / (x[1] - x[0]) * ((y[1].imag() - y[0].imag()) / (x[1] - x[0]) - yp1.imag()));
	}

	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the
     * decomposed factors*/

	for(i = 1; i < (n - 1); i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p   = sig * y2[i - 1].real() + 2.0;
		y2[i].real((sig - 1.0) / p);
		u[i].real((y[i + 1].real() - y[i].real()) / (x[i + 1] - x[i])
				  - (y[i].real() - y[i - 1].real()) / (x[i] - x[i - 1]));
		u[i].real((6.0 * u[i].real() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].real()) / p);
		p = sig * y2[i - 1].imag() + 2.0;
		y2[i].imag((sig - 1.0) / p);
		u[i].imag((y[i + 1].imag() - y[i].imag()) / (x[i + 1] - x[i])
				  - (y[i].imag() - y[i - 1].imag()) / (x[i] - x[i - 1]));
		u[i].imag((6.0 * u[i].imag() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].imag()) / p);
	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn.real() > 0.99e30) {
		qn = 0.;
		un = 0.;
	} else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.real() - (y[n - 1].real() - y[n - 2].real()) / (x[n - 1] - x[n - 2]));
	}
	y2[n - 1].real((un - qn * u[n - 2].real()) / (qn * y2[n - 2].real() + 1.0));
	if(ypn.imag() > 0.99e30) {
		qn = 0.;
		un = 0.;
	} else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.imag() - (y[n - 1].imag() - y[n - 2].imag()) / (x[n - 1] - x[n - 2]));
	}
	y2[n - 1].imag((un - qn * u[n - 2].imag()) / (qn * y2[n - 2].imag() + 1.0));

	/* This is the backsubstitution loop of the tridiagonal algorithm */

	for(k = n - 2; k >= 0; k--) {
		y2[k].real(y2[k].real() * y2[k + 1].real() + u[k].real());
		y2[k].imag(y2[k].imag() * y2[k + 1].imag() + u[k].imag());
	}

	return y2;
}

#endif /* Utils_H_ */
