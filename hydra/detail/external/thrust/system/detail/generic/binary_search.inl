/*
 *  Copyright 2008-2013 NVIDIA Corporation
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


/*! \file binary_search.inl
 *  \brief Inline file for binary_search.h
 */

#pragma once

#include <hydra/detail/external/thrust/detail/config.h>
#include <hydra/detail/external/thrust/distance.h>
#include <hydra/detail/external/thrust/binary_search.h>
#include <hydra/detail/external/thrust/iterator/zip_iterator.h>
#include <hydra/detail/external/thrust/iterator/iterator_traits.h>
#include <hydra/detail/external/thrust/binary_search.h>

#include <hydra/detail/external/thrust/for_each.h>
#include <hydra/detail/external/thrust/detail/function.h>
#include <hydra/detail/external/thrust/system/detail/generic/scalar/binary_search.h>

#include <hydra/detail/external/thrust/detail/temporary_array.h>
#include <hydra/detail/external/thrust/detail/type_traits.h>

HYDRA_EXTERNAL_NAMESPACE_BEGIN  namespace thrust
{
namespace detail
{


// XXX WAR circular #inclusion with this forward declaration
template<typename,typename> class temporary_array;


} // end detail
namespace system
{
namespace detail
{
namespace generic
{
namespace detail
{


// short names to avoid nvcc bug
struct lbf
{
  template<typename RandomAccessIterator, typename T, typename StrictWeakOrdering>
  __hydra_host__ __hydra_device__
  typename thrust::iterator_traits<RandomAccessIterator>::difference_type
    operator()(RandomAccessIterator begin, RandomAccessIterator end, const T& value, StrictWeakOrdering comp)
  {
    return thrust::system::detail::generic::scalar::lower_bound(begin, end, value, comp) - begin;
  }
};


struct ubf
{
  template<typename RandomAccessIterator, typename T, typename StrictWeakOrdering>
  __hydra_host__ __hydra_device__
  typename thrust::iterator_traits<RandomAccessIterator>::difference_type
    operator()(RandomAccessIterator begin, RandomAccessIterator end, const T& value, StrictWeakOrdering comp)
  {
    return thrust::system::detail::generic::scalar::upper_bound(begin, end, value, comp) - begin;
  }
};


struct bsf
{
  template<typename RandomAccessIterator, typename T, typename StrictWeakOrdering>
  __hydra_host__ __hydra_device__
  bool operator()(RandomAccessIterator begin, RandomAccessIterator end, const T& value, StrictWeakOrdering comp)
  {
    RandomAccessIterator iter = thrust::system::detail::generic::scalar::lower_bound(begin, end, value, comp);
    
    thrust::detail::wrapped_function<StrictWeakOrdering,bool> wrapped_comp(comp);
    
    return iter != end && !wrapped_comp(value, *iter);
  }
};


template<typename ForwardIterator, typename StrictWeakOrdering, typename BinarySearchFunction>
struct binary_search_functor
{
  ForwardIterator begin;
  ForwardIterator end;
  StrictWeakOrdering comp;
  BinarySearchFunction func;
  
  __hydra_host__ __hydra_device__
  binary_search_functor(ForwardIterator begin, ForwardIterator end, StrictWeakOrdering comp, BinarySearchFunction func)
    : begin(begin), end(end), comp(comp), func(func) {}
  
  template<typename Tuple>
  __hydra_host__ __hydra_device__
  void operator()(Tuple t)
  {
    thrust::get<1>(t) = func(begin, end, thrust::get<0>(t), comp);
  }
}; // binary_search_functor


// Vector Implementation
template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator, typename StrictWeakOrdering, typename BinarySearchFunction>
__hydra_host__ __hydra_device__
OutputIterator binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                             ForwardIterator begin, 
                             ForwardIterator end,
                             InputIterator values_begin, 
                             InputIterator values_end,
                             OutputIterator output,
                             StrictWeakOrdering comp,
                             BinarySearchFunction func)
{
  thrust::for_each(exec,
                   thrust::make_zip_iterator(thrust::make_tuple(values_begin, output)),
                   thrust::make_zip_iterator(thrust::make_tuple(values_end, output + thrust::distance(values_begin, values_end))),
                   detail::binary_search_functor<ForwardIterator, StrictWeakOrdering, BinarySearchFunction>(begin, end, comp, func));
  
  return output + thrust::distance(values_begin, values_end);
}

   

// Scalar Implementation
template<typename OutputType, typename DerivedPolicy, typename ForwardIterator, typename T, typename StrictWeakOrdering, typename BinarySearchFunction>
__hydra_host__ __hydra_device__
OutputType binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                         ForwardIterator begin,
                         ForwardIterator end,
                         const T& value, 
                         StrictWeakOrdering comp,
                         BinarySearchFunction func)
{
  // use the vectorized path to implement the scalar version
  
  // allocate device buffers for value and output
  thrust::detail::temporary_array<T,DerivedPolicy>          d_value(exec,1);
  thrust::detail::temporary_array<OutputType,DerivedPolicy> d_output(exec,1);
  
  // copy value to device
  d_value[0] = value;
  
  // perform the query
  thrust::system::detail::generic::detail::binary_search(exec, begin, end, d_value.begin(), d_value.end(), d_output.begin(), comp, func);
  
  // copy result to host and return
  return d_output[0];
}


// this functor differs from thrust::less<T>
// because it allows the types of lhs & rhs to differ
// which is required by the binary search functions
// XXX use C++14 thrust::less<> when it's ready
struct binary_search_less
{
  template<typename T1, typename T2>
  __hydra_host__ __hydra_device__
  bool operator()(const T1& lhs, const T2& rhs) const
  {
    return lhs < rhs;
  }
};

   
} // end namespace detail


//////////////////////
// Scalar Functions //
//////////////////////


template<typename DerivedPolicy, typename ForwardIterator, typename T>
__hydra_host__ __hydra_device__
ForwardIterator lower_bound(thrust::execution_policy<DerivedPolicy> &exec,
                            ForwardIterator begin,
                            ForwardIterator end,
                            const T& value)
{
  namespace p = thrust::placeholders;
  return thrust::lower_bound(exec, begin, end, value, detail::binary_search_less());
}

template<typename DerivedPolicy, typename ForwardIterator, typename T, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
ForwardIterator lower_bound(thrust::execution_policy<DerivedPolicy> &exec,
                            ForwardIterator begin,
                            ForwardIterator end,
                            const T& value, 
                            StrictWeakOrdering comp)
{
  typedef typename thrust::iterator_traits<ForwardIterator>::difference_type difference_type;
  
  return begin + detail::binary_search<difference_type>(exec, begin, end, value, comp, detail::lbf());
}


template<typename DerivedPolicy, typename ForwardIterator, typename T>
__hydra_host__ __hydra_device__
ForwardIterator upper_bound(thrust::execution_policy<DerivedPolicy> &exec,
                            ForwardIterator begin,
                            ForwardIterator end,
                            const T& value)
{
  namespace p = thrust::placeholders;
  return thrust::upper_bound(exec, begin, end, value, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename T, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
ForwardIterator upper_bound(thrust::execution_policy<DerivedPolicy> &exec,
                            ForwardIterator begin,
                            ForwardIterator end,
                            const T& value, 
                            StrictWeakOrdering comp)
{
  typedef typename thrust::iterator_traits<ForwardIterator>::difference_type difference_type;
  
  return begin + detail::binary_search<difference_type>(exec, begin, end, value, comp, detail::ubf());
}


template<typename DerivedPolicy, typename ForwardIterator, typename T>
__hydra_host__ __hydra_device__
bool binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                   ForwardIterator begin,
                   ForwardIterator end,
                   const T& value)
{
  return thrust::binary_search(exec, begin, end, value, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename T, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
bool binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                   ForwardIterator begin,
                   ForwardIterator end,
                   const T& value, 
                   StrictWeakOrdering comp)
{
  return detail::binary_search<bool>(exec, begin, end, value, comp, detail::bsf());
}


//////////////////////
// Vector Functions //
//////////////////////


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator>
__hydra_host__ __hydra_device__
OutputIterator lower_bound(thrust::execution_policy<DerivedPolicy> &exec,
                           ForwardIterator begin, 
                           ForwardIterator end,
                           InputIterator values_begin, 
                           InputIterator values_end,
                           OutputIterator output)
{
  namespace p = thrust::placeholders;
  return thrust::lower_bound(exec, begin, end, values_begin, values_end, output, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
OutputIterator lower_bound(thrust::execution_policy<DerivedPolicy> &exec,
                           ForwardIterator begin, 
                           ForwardIterator end,
                           InputIterator values_begin, 
                           InputIterator values_end,
                           OutputIterator output,
                           StrictWeakOrdering comp)
{
  return detail::binary_search(exec, begin, end, values_begin, values_end, output, comp, detail::lbf());
}


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator>
__hydra_host__ __hydra_device__
OutputIterator upper_bound(thrust::execution_policy<DerivedPolicy> &exec,
                           ForwardIterator begin, 
                           ForwardIterator end,
                           InputIterator values_begin, 
                           InputIterator values_end,
                           OutputIterator output)
{
  namespace p = thrust::placeholders;
  return thrust::upper_bound(exec, begin, end, values_begin, values_end, output, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
OutputIterator upper_bound(thrust::execution_policy<DerivedPolicy> &exec,
                           ForwardIterator begin, 
                           ForwardIterator end,
                           InputIterator values_begin, 
                           InputIterator values_end,
                           OutputIterator output,
                           StrictWeakOrdering comp)
{
  return detail::binary_search(exec, begin, end, values_begin, values_end, output, comp, detail::ubf());
}


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator>
__hydra_host__ __hydra_device__
OutputIterator binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                             ForwardIterator begin, 
                             ForwardIterator end,
                             InputIterator values_begin, 
                             InputIterator values_end,
                             OutputIterator output)
{
  namespace p = thrust::placeholders;
  return thrust::binary_search(exec, begin, end, values_begin, values_end, output, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename InputIterator, typename OutputIterator, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
OutputIterator binary_search(thrust::execution_policy<DerivedPolicy> &exec,
                             ForwardIterator begin, 
                             ForwardIterator end,
                             InputIterator values_begin, 
                             InputIterator values_end,
                             OutputIterator output,
                             StrictWeakOrdering comp)
{
  return detail::binary_search(exec, begin, end, values_begin, values_end, output, comp, detail::bsf());
}


template<typename DerivedPolicy, typename ForwardIterator, typename LessThanComparable>
__hydra_host__ __hydra_device__
thrust::pair<ForwardIterator,ForwardIterator>
equal_range(thrust::execution_policy<DerivedPolicy> &exec,
            ForwardIterator first,
            ForwardIterator last,
            const LessThanComparable &value)
{
  return thrust::equal_range(exec, first, last, value, detail::binary_search_less());
}


template<typename DerivedPolicy, typename ForwardIterator, typename T, typename StrictWeakOrdering>
__hydra_host__ __hydra_device__
thrust::pair<ForwardIterator,ForwardIterator>
equal_range(thrust::execution_policy<DerivedPolicy> &exec,
            ForwardIterator first,
            ForwardIterator last,
            const T &value,
            StrictWeakOrdering comp)
{
  ForwardIterator lb = thrust::lower_bound(exec, first, last, value, comp);
  ForwardIterator ub = thrust::upper_bound(exec, first, last, value, comp);
  return thrust::make_pair(lb, ub);
}


} // end namespace generic
} // end namespace detail
} // end namespace system
} // end HYDRA_EXTERNAL_NAMESPACE_BEGIN  namespace thrust

HYDRA_EXTERNAL_NAMESPACE_END
