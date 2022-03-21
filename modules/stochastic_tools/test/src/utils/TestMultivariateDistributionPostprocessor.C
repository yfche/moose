//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TestMultivariateDistributionPostprocessor.h"
#include "MultivariateDistribution.h"

registerMooseObject("StochasticToolsTestApp", TestMultivariateDistributionPostprocessor);

InputParameters
TestMultivariateDistributionPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<MultivariateDistributionName>(
      "distribution", "The distribution which supplies the postprocessor value.");
  params.addRequiredParam<Eigen::VectorXd>(
      "value", "A value to pass to the cdf, pdf, or quantile function of the given distribution.");

  MooseEnum method("pdf");  // "cdf pdf quantile"
  params.addParam<MooseEnum>("method", method, "The distribution method to call.");
  return params;
}

TestMultivariateDistributionPostprocessor::TestMultivariateDistributionPostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _distribution(getMultivariateDistribution("distribution")),
    _value(getParam<Eigen::VectorXd>("value")),
    _distribution_method(getParam<MooseEnum>("method"))
{
}

PostprocessorValue
TestMultivariateDistributionPostprocessor::getValue()
{
  if (_distribution_method == "pdf")
    return _distribution.pdf(_value);
  // else if (_distribution_method == "cdf")
  //   return _distribution.cdf(_value);
  // else if (_distribution_method == "quantile")
  //   return _distribution.quantile(_value);
  mooseError("This should be possible to reach, what did you do!");
}
