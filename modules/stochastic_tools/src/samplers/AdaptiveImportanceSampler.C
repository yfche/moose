//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdaptiveImportanceSampler.h"
#include "AdaptiveMonteCarloUtils.h"
#include "Distribution.h"
#include "Normal.h"
#include "Uniform.h"

registerMooseObjectAliased("StochasticToolsApp", AdaptiveImportanceSampler, "AdaptiveImportance");

InputParameters
AdaptiveImportanceSampler::validParams()
{
  InputParameters params = Sampler::validParams();
  params.addClassDescription("Adaptive Importance Sampler.");
  params.addRequiredParam<std::vector<DistributionName>>(
      "distributions",
      "The distribution names to be sampled, the number of distributions provided defines the "
      "number of columns per matrix.");
  params.addRequiredParam<ReporterName>("inputs_reporter", "Reporter with input parameters.");
  params.addRequiredParam<std::vector<Real>>("proposal_std",
                                             "Standard deviations of the proposal distributions");
  params.addRequiredParam<Real>("output_limit", "Limiting values of the VPPs");
  params.addRequiredParam<std::vector<Real>>(
      "initial_values", "Initial input values to get the importance sampler started");
  params.addRequiredParam<int>("num_samples_train",
                               "Number of samples to learn the importance distribution");
  params.addRequiredParam<Real>(
      "std_factor", "Factor to be multiplied to the standard deviation of the importance samples");
  params.addParam<bool>("use_absolute_value", false, "Use absolute value of the sub app output");
  params.addParam<unsigned int>(
      "num_random_seeds",
      100000,
      "Initialize a certain number of random seeds. Change from the default only if you have to.");
  MooseEnum method("MH AdaptiveMH", "MH");
  params.addParam<MooseEnum>(
      "method", method, "The method to generate new samples in Markov chain.");
  // // Yifeng
  // params.addParam<unsigned int>("num_samples_adaption", 10,
  //                               "Number of samples to generate before performing adaption.");

  return params;
}

AdaptiveImportanceSampler::AdaptiveImportanceSampler(const InputParameters & parameters)
  : Sampler(parameters),
    ReporterInterface(this),
    _proposal_std(getParam<std::vector<Real>>("proposal_std")),
    _initial_values(getParam<std::vector<Real>>("initial_values")),
    _output_limit(getParam<Real>("output_limit")),
    _num_samples_train(getParam<int>("num_samples_train")),
    _std_factor(getParam<Real>("std_factor")),
    _use_absolute_value(getParam<bool>("use_absolute_value")),
    _num_random_seeds(getParam<unsigned int>("num_random_seeds")),
    _sampling_method(getParam<MooseEnum>("method")),
    // _num_samples_adaption(getParam<unsigned int>("num_samples_adaption")), // Yifeng
    _step(getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")->timeStep()),
    _inputs(getReporterValue<std::vector<std::vector<Real>>>("inputs_reporter"))
{
  // Filling the `distributions` vector with the user-provided distributions.
  for (const DistributionName & name : getParam<std::vector<DistributionName>>("distributions"))
    _distributions.push_back(&getDistributionByName(name));

  /* Adaptive Importance Sampling (AdaptiveImportanceSampler) relies on a Markov Chain Monte Carlo
     (MCMC) algorithm. As such, in MOOSE, any use of MCMC algorithms requires that the `num_steps`
     parameter in the main App's executioner would control the total number of samples. Therefore,
     the `num_rows` parameter typically used by exisiting non-MCMC samplers to set the total number
     of samples has no use here and is fixed to 1.*/
  setNumberOfRows(1);

  // Setting the number of columns in the sampler matrix (equal to the number of distributions).
  setNumberOfCols(_distributions.size());

  /* `inputs_sto` is a member variable that aids in forming the importance distribution.
     One dimension of this variable is equal to the number of distributions. The other dimension
     of the variable, at the last step, is equal to the number of samples the user desires.*/
  _inputs_sto.resize(_distributions.size());

  // Mapping all the input distributions to a standard normal space
  for (unsigned int i = 0; i < _distributions.size(); ++i)
    _inputs_sto[i].push_back(Normal::quantile(_distributions[i]->cdf(_initial_values[i]), 0, 1));

  /* `prev_value` is a member variable for tracking the previously accepted samples in the
     MCMC algorithm and proposing the next sample.*/
  _prev_value.resize(_distributions.size());

  // `check_step` is a member variable for ensuring that the MCMC algorithm proceeds in a sequential
  // fashion.
  _check_step = 0;

  // Storage for means of input values for proposing the next sample
  _mean_sto.resize(_distributions.size());

  // Storage for standard deviations of input values for proposing the next sample
  _std_sto.resize(_distributions.size());

  setNumberOfRandomSeeds(_num_random_seeds);

  // AdaptiveOptimalScaling
  if (_sampling_method == "AdaptiveMH")
  {
    _opt_rate = 0.44;
    _COV = 1 * RealEigenMatrix::Identity(_distributions.size(), _distributions.size());
    Moose::out << "_COV initial: " << _COV << std::endl;
    _C0 = 1;
    _Log_lambda.resize(1);
    _Gamma.resize(1);
    // _MU.resize(_distributions.size());
    _MU = RealEigenVector::Zero(_distributions.size());
    Moose::out << "_MU initial: " << _MU << std::endl;
    _Log_lambda[0] = std::log(std::pow(2.38, 2) / _distributions.size());
    _Gamma[0] = _C0;
    // std::fill(_MU.begin(), _MU.end(), 0.0);
  }

}

const std::vector<Real> &
AdaptiveImportanceSampler::getInitialValues() const
{
  return _initial_values;
}

const int &
AdaptiveImportanceSampler::getNumSamplesTrain() const
{
  return _num_samples_train;
}

const bool &
AdaptiveImportanceSampler::getUseAbsoluteValue() const
{
  return _use_absolute_value;
}

const Real &
AdaptiveImportanceSampler::getOutputLimit() const
{
  return _output_limit;
}

Real
AdaptiveImportanceSampler::computeSample(dof_id_type /*row_index*/, dof_id_type col_index)
{
  const bool sample = _step > 1 && col_index == 0 && _check_step != _step;
  if (_step <= _num_samples_train)
  {
    /* This is the importance distribution training step. Markov Chains are set up
       to sample from the importance region or the failure region using the Metropolis
       algorithm. Given that the previous sample resulted in a model failure, the next
       sample is proposed such that it is very likely to result in a model failure as well.
       The `initial_values` and `proposal_std` parameters provided by the user affects the
       formation of the importance distribution. */
    if (_step == 1)
    {
      for (dof_id_type i = 0; i < _distributions.size(); ++i)
        _MU(i) = _inputs_sto[i].back();
    }


    if (sample)
    {
      if (_sampling_method == "MH")
      {
        std::vector<Real> new_sample = AdaptiveMonteCarloUtils::proposeNewSampleMH(
            _distributions, getRand(_step), _inputs, _inputs_sto);
        for (dof_id_type i = 0; i < _distributions.size(); ++i)
          _inputs_sto[i].push_back(new_sample[i]);
      }

      else if (_sampling_method == "AdaptiveMH")
      {
        // Moose::out << "inputs0: " << _inputs[0].size() << std::endl;
        // Moose::out << "inputs1: " << _inputs[1].size() << std::endl;

        Real Lambda = std::exp(_Log_lambda.back());
        Moose::out << "Lambda: " << Lambda << std::endl;
        RealEigenVector current_sample(_distributions.size());
        RealEigenVector propose_sample(_distributions.size());
        for (dof_id_type i = 0; i < _distributions.size(); ++i)
        {
          current_sample(i) = _inputs_sto[i].back();
          propose_sample(i) = _distributions[i]->cdf(_inputs[i][0]);
          // current_sample(i) = Normal::quantile(_distributions[i]->cdf(_inputs[i].back()), 0, 1);
        }

        Moose::out << "current_sample: " << current_sample << std::endl;
        // RealEigenVector generate_sample = AdaptiveMonteCarloUtils::sampleFromMultivariateNormal(current_sample, Lambda * _COV);
        RealEigenVector generate_sample = AdaptiveMonteCarloUtils::sampleFromMultivariateNormal(propose_sample, Lambda * _COV);

        RealEigenVector mean(_distributions.size());
        mean.setZero();
        RealEigenMatrix covariance = RealEigenMatrix::Identity(_distributions.size(), _distributions.size());

        Real prob_eval_new = AdaptiveMonteCarloUtils::MVNpdf(generate_sample, mean, covariance);
        Real prob_eval_old = AdaptiveMonteCarloUtils::MVNpdf(current_sample, mean, covariance);
        Real acceptance_ratio = std::log(prob_eval_new) - std::log(prob_eval_old);


        RealEigenVector new_sample(_distributions.size());
        if (acceptance_ratio > std::log(getRand(_step)))  // previously ">"
          new_sample = generate_sample;
        else
          new_sample = current_sample;

        for (dof_id_type i = 0; i < _distributions.size(); ++i)
          _inputs_sto[i].push_back(new_sample[i]);

        // Adaption
        // _Gamma.push_back(_C0 * std::pow(_step, -1/(1+ std::exp(_Log_lambda[0]))));
        _Gamma.push_back(_C0 * std::pow(_step, -1));
        Moose::out << "_Gamma: " << _Gamma.back() << std::endl;
        Real ALPHA = std::min<Real>(1, exp(-acceptance_ratio));
        Moose::out << "ALPHA: " << ALPHA << std::endl;
        _Log_lambda.push_back(_Log_lambda.back() + _Gamma.back() * (ALPHA - _opt_rate));
        Moose::out << "_Log_lambda: " << _Log_lambda.back() << std::endl;
        _MU = _MU + _Gamma.back() * (new_sample - _MU);

        RealEigenMatrix SS = (new_sample - _MU) * (new_sample - _MU).transpose();
        // Moose::out << "(new_sample - _MU).transpose: " << (new_sample - _MU).transpose() << std::endl;
        // Moose::out << "new_sample - _MU: " << new_sample - _MU << std::endl;
        Moose::out << "new_sample: " << new_sample << std::endl;
        Moose::out << "_MU: " << _MU << std::endl;
        Moose::out << "SS: " << SS << std::endl;
        _COV = _COV + _Gamma.back() * (SS - _COV);

        Moose::out << "_COV: " << _COV << std::endl;
      }

    }
    _prev_value[col_index] =
        Normal::quantile(getRand(_step), _inputs_sto[col_index].back(), _proposal_std[col_index]);
  }
  else if (sample)
  {
    /* This is the importance sampling step using the importance distribution created
       in the previous step. Once the importance distribution is known, sampling from
       it is similar to a regular Monte Carlo sampling. */
    for (dof_id_type i = 0; i < _distributions.size(); ++i)
    {
      if (_step == _num_samples_train + 1)
      {
        _mean_sto[i] = AdaptiveMonteCarloUtils::computeMean(_inputs_sto[i], 1);
        _std_sto[i] = AdaptiveMonteCarloUtils::computeSTD(_inputs_sto[i], 1);
      }
      _prev_value[i] = (Normal::quantile(getRand(_step), _mean_sto[i], _std_factor * _std_sto[i]));
    }
  }

  _check_step = _step;
  return _distributions[col_index]->quantile(Normal::cdf(_prev_value[col_index], 0, 1));
}
