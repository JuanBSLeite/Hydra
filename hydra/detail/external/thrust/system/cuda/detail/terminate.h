/*
 *  Copyright 2008-2012 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#pragma once

#include <hydra/detail/external/thrust/detail/config.h>
#include <hydra/detail/external/thrust/system/cuda/detail/bulk.h>

HYDRA_EXTERNAL_NAMESPACE_BEGIN  namespace thrust
{
namespace system
{
namespace cuda
{
namespace detail
{


inline __hydra_device__
void terminate()
{
  thrust::system::cuda::detail::bulk_::detail::terminate();
}


__hydra_host__ __hydra_device__
inline void terminate_with_message(const char* message)
{
  thrust::system::cuda::detail::bulk_::detail::terminate_with_message(message);
}


} // end detail
} // end cuda
} // end system
} // end thrust

HYDRA_EXTERNAL_NAMESPACE_END
