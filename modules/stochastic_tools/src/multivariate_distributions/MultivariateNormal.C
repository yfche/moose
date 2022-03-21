//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultivariateNormal.h"
#include "math.h"
#include "libmesh/utility.h"

registerMooseObject("StochasticToolsApp", MultivariateNormal);

InputParameters
MultivariateNormal::validParams()
{
  InputParameters params = MultivariateDistribution::validParams();
  // InputParameters params = MooseObject::validParams();
  params.addClassDescription("MultivariateNormal distribution");
  params.addRequiredParam<Eigen::VectorXd>("mean", "Mean (or expectation) of the distribution.");
  params.addRequiredRangeCheckedParam<Eigen::MatrixXd>(
      "covariance_matrix", "covariance_matrix>0", "Covariance matrix of the distribution ");
  // params.registerBase("MultivariateNormalDistribution");
  // params.registerSystemAttributeName("MultivariateNormalDistribution");
  return params;
}

MultivariateNormal::MultivariateNormal(const InputParameters & parameters): MultivariateDistribution(parameters),
    _mean(getParam<Eigen::VectorXd>("mean")),
    _covariance_matrix(getParam<Eigen::MatrixXd>("covariance_matrix"))
{
}

Real
MultivariateNormal::pdf(const Eigen::VectorXd & x, const Eigen::VectorXd & mean, const Eigen::MatrixXd & covariance_matrix)
{
  Real n = x.rows();
  Real sqrt2pi = std::sqrt(2 * M_PI);
  Real quadform = (x - mean).transpose() * covariance_matrix.inverse() * (x - mean);
  Real norm = std::pow(sqrt2pi, -n) * std::pow(covariance_matrix.determinant(), -0.5);
  return norm * exp(-0.5 * quadform);
}

Real
MultivariateNormal::pdf(const Eigen::VectorXd & x) const
{
  return pdf(x, _mean, _covariance_matrix);
}

//
// Real
// MultivariateNormal::cdf(const Real & x, const Real & mean, const Real & std_dev)
// {
//   return 0.5 * (1.0 + std::erf((x - mean) / (std_dev * std::sqrt(2.0))));
// }
//
// Real
// MultivariateNormal::quantile(const Real & p, const Real & mean, const Real & std_dev)
// {
//   const Real x = (p < 0.5 ? p : 1.0 - p);
//   const Real y = std::sqrt(-2.0 * std::log(x));
//   const Real y2 = y * y;
//   const Real y3 = y2 * y;
//   const Real y4 = y3 * y;
//   const Real sgn = (p - 0.5 < 0.0 ? -1.0 : 1.0);
//   const Real Zp = sgn * (y + (_a[0] + _a[1] * y + _a[2] * y2 + _a[3] * y3 + _a[4] * y4) /
//                                  (_b[0] + _b[1] * y + _b[2] * y2 + _b[3] * y3 + _b[4] * y4));
//   return Zp * std_dev + mean;
// }
//

//
// Real
// MultivariateNormal::cdf(const Real & x) const
// {
//   return 0;
//   // return cdf(x, _mean, _standard_deviation);
// }
//
// Real
// MultivariateNormal::quantile(const Real & p) const
// {
//   return 0;
//   // return quantile(p, _mean, _standard_deviation);
// }
