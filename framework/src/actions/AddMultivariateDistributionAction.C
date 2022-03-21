//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddMultivariateDistributionAction.h"
#include "FEProblem.h"

registerMooseAction("MooseApp", AddMultivariateDistributionAction, "add_multivariate_distribution");

InputParameters
AddMultivariateDistributionAction::validParams()
{
  InputParameters params = MooseObjectAction::validParams();
  params.addClassDescription("Add a MultivariateDistribution object to the simulation.");
  return params;
}

AddMultivariateDistributionAction::AddMultivariateDistributionAction(InputParameters params) : MooseObjectAction(params) {}

void
AddMultivariateDistributionAction::act()
{
  _problem->addMultivariateDistribution(_type, _name, _moose_object_pars);
}
