//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InputParameters.h"
#include "ParallelUniqueId.h"
#include "MooseObject.h"

// Forward declarations
class MultivariateDistribution;
class FEProblemBase;

/**
 * Interface for objects that need to use distributions
 *
 * Inherit from this class at a very low level to make the getMultivariateDistribution method available.
 */
class MultivariateDistributionInterface
{
public:
  static InputParameters validParams();

  MultivariateDistributionInterface(const MooseObject * moose_object);

  ///@{
  /**
   * Get a multivariate distribution with a given name
   * @param name The name of the parameter key of the multivariate distribution to retrieve
   * @return The multivariate distribution with name associated with the parameter 'name'
   */
  const MultivariateDistribution & getMultivariateDistribution(const std::string & name) const;

  template <typename T>
  const T & getMultivariateDistribution(const std::string & name) const;
  ///@}

  ///@{
  /**
   * Get a multivariate distribution with a given name
   * @param name The name of the multivariate distribution to retrieve
   * @return The multivariate distribution with name 'name'
   */
  const MultivariateDistribution & getMultivariateDistributionByName(const MultivariateDistributionName & name) const;

  template <typename T>
  const T & getMultivariateDistributionByName(const std::string & name) const;
  ///@}

private:
  /// Parameters of the object with this interface
  const InputParameters & _dni_params;

  /// Reference to FEProblemBase instance
  FEProblemBase & _dni_feproblem;

  /// Pointer to the MooseObject
  const MooseObject * const _dni_moose_object_ptr;
};

template <typename T>
const T &
MultivariateDistributionInterface::getMultivariateDistribution(const std::string & name) const
{
  try
  {
    const T & dist = dynamic_cast<const T &>(getMultivariateDistribution(name));
    return dist;
  }
  catch (std::bad_cast & exception)
  {
    MultivariateDistributionName dist_name = _dni_params.get<MultivariateDistributionName>(name);
    mooseError("The '",
               _dni_moose_object_ptr->name(),
               "' object failed to retrieve '",
               dist_name,
               "' multivariate distribution with the desired type.");
  }
}

template <typename T>
const T &
MultivariateDistributionInterface::getMultivariateDistributionByName(const std::string & name) const
{
  try
  {
    const T & dist = dynamic_cast<const T &>(getMultivariateDistribution(name));
    return dist;
  }
  catch (std::bad_cast & exception)
  {
    mooseError("The '",
               _dni_moose_object_ptr->name(),
               "' object failed to retrieve '",
               name,
               "' multivariate distribution with the desired type.");
  }
}
