//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultivariateDistributionInterface.h"
#include "MultivariateDistribution.h"
#include "SubProblem.h"
#include "MooseTypes.h"
#include "FEProblemBase.h"

InputParameters
MultivariateDistributionInterface::validParams()
{
  InputParameters params = emptyInputParameters();
  return params;
}

MultivariateDistributionInterface::MultivariateDistributionInterface(const MooseObject * moose_object)
  : _dni_params(moose_object->parameters()),
    _dni_feproblem(*_dni_params.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    _dni_moose_object_ptr(moose_object)
{
}

const MultivariateDistribution &
MultivariateDistributionInterface::getMultivariateDistribution(const std::string & name) const
{
  MultivariateDistributionName dist_name = _dni_params.get<MultivariateDistributionName>(name);
  return _dni_feproblem.getMultivariateDistribution(dist_name);
}

const MultivariateDistribution &
MultivariateDistributionInterface::getMultivariateDistributionByName(const MultivariateDistributionName & name) const
{
  return _dni_feproblem.getMultivariateDistribution(name);
}
