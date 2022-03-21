//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObject.h"
#include "PerfGraphInterface.h"

/**
 * All MultivarateDistributions should inherit from this class
 */
class MultivariateDistribution : public MooseObject, public PerfGraphInterface
{
public:
  static InputParameters validParams();

  MultivariateDistribution(const InputParameters & parameters);
  /**
   * Compute the probability with given probability distribution function (PDF) at x
   */
  virtual Real pdf(const Eigen::VectorXd & x) const = 0;

  // /**
  //  * Compute the cumulative probability with given cumulative probability distribution (CDF) at x
  //  */
  // virtual Real cdf(const Real & x) const = 0;
  //
  // /**
  //  * Compute the inverse CDF (quantile function) value for given variable value y
  //  */
  // virtual Real quantile(const Real & y) const = 0;
  //
  /**
   * Compute the median of the distribution
   */
  virtual Real median() const;
};
