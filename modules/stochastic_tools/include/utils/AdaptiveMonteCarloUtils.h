//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseUtils.h"
#include "Sampler.h"

namespace AdaptiveMonteCarloUtils
{
/* AdaptiveMonteCarloUtils contains functions that are used across the Adaptive Monte
 Carlo set of algorithms.*/

/**
 * compute the standard deviation of a data vector by only considering values from
 * a specific index.
 *
 * @param the data vector
 * @param the starting index
 */
Real computeSTD(const std::vector<Real> & data, const unsigned int & start_index);

/**
 * compute the mean of a data vector by only considering values from
 * a specific index.
 *
 * @param the data vector
 * @param the starting index
 */
Real computeMean(const std::vector<Real> & data, const unsigned int & start_index);

/**
 * return input values corresponding to the largest po percentile output values.
 *
 ****** FOR PARALLEL SUBSET SIMULATION SAMPLER ********
 *
 * @param the input vector
 * @param the output vector
 * @param the number of samples per subset
 * @param the subset index
 * @param the subset intermediate failure probability
 */
std::vector<std::vector<Real>> sortInput(const std::vector<std::vector<Real>> & inputs,
                                         const std::vector<Real> & outputs,
                                         const unsigned int samplessub,
                                         const Real subset_prob);

/**
 * return the largest po percentile output values.
 *
 ****** FOR PARALLEL SUBSET SIMULATION SAMPLER ********
 *
 * @param the output vector
 * @param the number of samples per subset
 * @param the subset index
 * @param the subset intermediate failure probability
 */
std::vector<Real> sortOutput(const std::vector<Real> & outputs,
                             const unsigned int samplessub,
                             const Real subset_prob);

/**
 * return the minimum value in a vector.
 *
 * @param the data vector
 */
Real computeMin(const std::vector<Real> & data);

/**
 * return the absolute values in a vector.
 *
 * @param the data vector
 */
std::vector<Real> computeVectorABS(const std::vector<Real> & data);

/**
 * propose a new sample using ComponentWise-MH.
 *
 * @param the current sample
 * @param the random seed
 */
Real proposeNewSampleComponentWiseMH(const Real x,
                                     const Real rnd1,
                                     const Real rnd2,
                                     const Real proposal_std = 1.0);

/**
 * propose a new sample w/ regular MH. (AdaptiveImportanceSampler)
 *
 * @param the input distributions
 * @param the random seed
 * @param the input vector from the reporter
 * @param the previously accepted samples
 */
std::vector<Real> proposeNewSampleMH(std::vector<const Distribution *> distributions,
                                     const Real rnd,
                                     const std::vector<std::vector<Real>> inputs,
                                     std::vector<std::vector<Real>> _inputs_sto);


/**
* Draw sample from multivariate normal distribution with specified mean and covariance matrix
*
* @param the number of iterations
* @param the mean vector
* @param the covariance matrix
*/
RealEigenVector sampleFromMultivariateNormal(RealEigenVector mean,
                                             RealEigenMatrix covariance_matrix,
                                             unsigned int nr_iterations = 1000);



/**
* Calculate pdf of multivariate normal distribution
*
* @param the current sample
* @param the mean vector
* @param the covariance matrix
*/
Real
MVNpdf(const RealEigenVector & x, const RealEigenVector & mean, const RealEigenMatrix & cov_mat);


} // namespace AdaptiveMonteCarloUtils
