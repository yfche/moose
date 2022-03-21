//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultivariateDistribution.h"
#include "MooseRandom.h"

InputParameters
MultivariateDistribution::validParams()
{
  InputParameters params = MooseObject::validParams();
  params.addRequiredParam<std::string>("type", "class/type name identifying the multivariate distribution");
  params.registerBase("MultivariateDistribution");
  params.registerSystemAttributeName("MultivariateDistribution");
  return params;
}

MultivariateDistribution::MultivariateDistribution(const InputParameters & parameters)
  : MooseObject(parameters), PerfGraphInterface(this)
{
}

Real
MultivariateDistribution::median() const
{
  mooseError("The multivariate distribution '", name(), "' must override the median method.");
  return 0;
}
