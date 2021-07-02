/** @file
 *****************************************************************************

 Imeplementation of interfaces for evaluation domains.

 See evaluation_domain.hpp .

 We currently implement, and select among, three types of domains:
 - "basic radix-2": the domain has size m = 2^k and consists of the m-th roots of unity
 - "extended radix-2": the domain has size m = 2^{k+1} and consists of "the m-th roots of unity" union "a coset"
 - "step radix-2": the domain has size m = 2^k + 2^r and consists of "the 2^k-th roots of unity" union "a coset of 2^r-th roots of unity"

 *****************************************************************************
 * @author     This file is part of libfqfft, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef GET_EVALUATION_DOMAIN_TCC_
#define GET_EVALUATION_DOMAIN_TCC_

#include <libfqfft/evaluation_domain/domains/arithmetic_sequence_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>
#include <libfqfft/evaluation_domain/domains/extended_radix2_domain.hpp>
#include <libfqfft/evaluation_domain/domains/geometric_sequence_domain.hpp>
#include <libfqfft/evaluation_domain/domains/step_radix2_domain.hpp>
#include <libfqfft/evaluation_domain/evaluation_domain.hpp>
#include <libfqfft/tools/exceptions.hpp>

namespace libfqfft {

template<typename FieldT>
std::shared_ptr<evaluation_domain<FieldT> > get_evaluation_domain(const size_t min_size)
{
    std::shared_ptr<evaluation_domain<FieldT> > result;

    const size_t big = 1ul<<(libff::log2(min_size)-1);
    const size_t small = min_size - big;
    const size_t rounded_small = (1ul<<libff::log2(small));

    if ( basic_radix2_domain<FieldT>::valid_for_size(min_size) ) {
        result.reset(new basic_radix2_domain<FieldT>(min_size));        
    }
    else if ( extended_radix2_domain<FieldT>::valid_for_size(min_size) ) {
        result.reset(new extended_radix2_domain<FieldT>(min_size));
    }
    else if ( step_radix2_domain<FieldT>::valid_for_size(min_size) ) {
        result.reset(new step_radix2_domain<FieldT>(min_size));
    }
    else if ( basic_radix2_domain<FieldT>::valid_for_size(big + rounded_small) ) {
        result.reset(new basic_radix2_domain<FieldT>(big + rounded_small));
    }
    else if ( extended_radix2_domain<FieldT>::valid_for_size(big + rounded_small) ) {
        result.reset(new extended_radix2_domain<FieldT>(big + rounded_small));        
    }
    else if ( step_radix2_domain<FieldT>::valid_for_size(big + rounded_small) ) {
        result.reset(new step_radix2_domain<FieldT>(big + rounded_small));
    }
    else if ( geometric_sequence_domain<FieldT>::valid_for_size(min_size) ) {
        result.reset(new geometric_sequence_domain<FieldT>(min_size));
    }
    else if ( arithmetic_sequence_domain<FieldT>::valid_for_size(min_size) ) {
        result.reset(new arithmetic_sequence_domain<FieldT>(min_size));
    }
    else {
        throw DomainSizeException("get_evaluation_domain: no matching domain");
    }

    return result;
}

} // libfqfft

#endif // GET_EVALUATION_DOMAIN_TCC_
