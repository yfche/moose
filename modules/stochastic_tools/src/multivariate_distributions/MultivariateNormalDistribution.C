//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultivariateNormalDistribution.h"

registerMooseObjectReplaced("StochasticToolsApp", MultivariateNormalDistribution, "03/14/2022 00:00", MultivariateNormal);

InputParameters
MultivariateNormalDistribution::validParams()
{
  return MultivariateNormal::validParams();
}

MultivariateNormalDistribution::MultivariateNormalDistribution(const InputParameters & parameters) : MultivariateNormal(parameters) {}
