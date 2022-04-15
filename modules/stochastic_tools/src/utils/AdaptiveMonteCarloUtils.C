//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdaptiveMonteCarloUtils.h"
#include "IndirectSort.h"
#include "libmesh/int_range.h"
#include "Normal.h"
#include "Sampler.h"

/* AdaptiveMonteCarloUtils contains functions that are used across the Adaptive Monte
 Carlo set of algorithms.*/

namespace AdaptiveMonteCarloUtils
{

Real
computeSTD(const std::vector<Real> & data, const unsigned int & start_index)
{
  if (data.size() < start_index)
    return 0.0;
  else
  {
    const Real mean = computeMean(data, start_index);
    const Real sq_diff =
        std::accumulate(data.begin() + start_index,
                        data.end(),
                        0.0,
                        [&mean](Real x, Real y) { return x + (y - mean) * (y - mean); });
    return std::sqrt(sq_diff / (data.size() - start_index));
  }
}

Real
computeMean(const std::vector<Real> & data, const unsigned int & start_index)
{
  if (data.size() < start_index)
    return 0.0;
  else
    return std::accumulate(data.begin() + start_index, data.end(), 0.0) /
           (data.size() - start_index);
}

std::vector<std::vector<Real>>
sortInput(const std::vector<std::vector<Real>> & inputs,
          const std::vector<Real> & outputs,
          const unsigned int samplessub,
          const Real subset_prob)
{
  std::vector<size_t> ind;
  Moose::indirectSort(outputs.begin(), outputs.end(), ind);

  std::vector<std::vector<Real>> out(inputs.size(), std::vector<Real>(samplessub * subset_prob));
  const size_t offset = std::ceil(samplessub * (1 - subset_prob));
  for (const auto & j : index_range(out))
    for (const auto & i : index_range(out[j]))
      out[j][i] = inputs[j][ind[i + offset]];

  return out;
}

std::vector<Real>
sortOutput(const std::vector<Real> & outputs, const unsigned int samplessub, const Real subset_prob)
{
  std::vector<size_t> ind;
  Moose::indirectSort(outputs.begin(), outputs.end(), ind);

  std::vector<Real> out(samplessub * subset_prob);
  const size_t offset = std::round(samplessub * (1 - subset_prob));
  for (const auto & i : index_range(out))
    out[i] = outputs[ind[i + offset]];

  return out;
}

Real
computeMin(const std::vector<Real> & data)
{
  return *std::min_element(data.begin(), data.end());
}

std::vector<Real>
computeVectorABS(const std::vector<Real> & data)
{
  std::vector<Real> data_abs(data.size());
  for (unsigned int i = 0; i < data.size(); ++i)
    data_abs[i] = std::abs(data[i]);
  return data_abs;
}

Real
proposeNewSampleComponentWiseMH(const Real x,
                                const Real rnd1,
                                const Real rnd2,
                                const Real proposal_std)
{
  const Real x_new = Normal::quantile(rnd1, x, 1.0);
  const Real acceptance_ratio = std::log(Normal::pdf(x_new, 0, 1)) - std::log(Normal::pdf(x, 0, 1));
  const Real new_sample = acceptance_ratio > std::log(rnd2) ? x_new : x;
  // Real val = Normal::cdf(new_sample, 0, 1);
  return new_sample;
}

std::vector<Real>
proposeNewSampleMH(std::vector<const Distribution *> distributions,
                   const Real rnd,
                   const std::vector<std::vector<Real>> inputs,
                   std::vector<std::vector<Real>> _inputs_sto)
{
  Real acceptance_ratio = 0.0;
  unsigned int d = distributions.size();
  std::vector<Real> _prev_value(d);
  std::vector<Real> new_proposal(d);

  for (dof_id_type i = 0; i < d; ++i)
  {
    _prev_value[i] = Normal::quantile(distributions[i]->cdf(inputs[i][0]), 0, 1);
    acceptance_ratio += std::log(Normal::pdf(_prev_value[i], 0, 1)) -
                        std::log(Normal::pdf(_inputs_sto[i].back(), 0, 1));
  }

  if (acceptance_ratio > std::log(rnd))
  {
    for (dof_id_type i = 0; i < d; ++i)
      new_proposal[i] = _prev_value[i];
  }
  else
  {
    for (dof_id_type i = 0; i < d; ++i)
      new_proposal[i] = _inputs_sto[i].back();
  }
  return new_proposal;
}


RealEigenVector sampleFromMultivariateNormal(RealEigenVector mean, RealEigenMatrix covariance_matrix, unsigned int nr_iterations) //, , Eigen::VectorXd mean, Eigen::MatrixXd covariance_matrix
{
  unsigned int d = mean.rows();  // sample dimension

  // Generate x from the N(0, I) Distribution
  RealEigenVector x(d);
  RealEigenVector sum(d);
  sum.setZero();

  for (unsigned int i = 0; i < nr_iterations; i++)
  {
    x.setRandom();
    x = 0.5 * (x + RealEigenVector::Ones(d));
    sum += x;
  }

  sum = sum - (static_cast<double>(nr_iterations) / 2) * RealEigenVector::Ones(d);
  x = sum / (std::sqrt(static_cast<double>(nr_iterations) / 12));

  // Find the eigen vectors of the covariance matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(covariance_matrix);
  RealEigenMatrix eigenvectors = eigen_solver.eigenvectors().real();

  // Find eigenvalues of the covariance matrix
  RealEigenMatrix eigenvalues = eigen_solver.eigenvalues().real().asDiagonal();

  // Find the transformation matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(eigenvalues);
  RealEigenMatrix sqrt_eigenvalues = es.operatorSqrt();
  RealEigenMatrix Q = eigenvectors * sqrt_eigenvalues;

  return Q * x + mean;
}


Real
MVNpdf(const RealEigenVector & x, const RealEigenVector & mean, const RealEigenMatrix & cov_mat)
{
  // Compute PDF here
  const Real norm = 1.0 / std::sqrt(MathUtils::pow(2.0 * M_PI, x.size()) * cov_mat.determinant());
  const Real quadform = (x - mean).transpose() * cov_mat.inverse() * (x - mean);
  return norm * std::exp(-0.5 * quadform);
}


} // namespace AdaptiveMonteCarloUtils
