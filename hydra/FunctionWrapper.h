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
 * FunctionWrapper.h
 *
 *  Created on: 12/07/2016
 *      Author: Antonio Augusto Alves Junior
 */


#ifndef FUNCTIONWRAPPER_H_
#define FUNCTIONWRAPPER_H_

#if __cplusplus < 201402L // c++11
#include <hydra/detail/FunctionWrapper_cpp11.inl>
#else // c++14
#include <hydra/detail/FunctionWrapper_cpp14.inl>
#endif

#endif /* FUNCTIONWRAPPER_H_ */
