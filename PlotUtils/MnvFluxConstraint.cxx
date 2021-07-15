/*
 * MnvFluxConstraint.cxx
 *
 *  Implementation file for MnvFluxConstraint classes.
 *
 *  Created on: Jul 14, 2014
 *      Author: J. Wolcott <jwolcott@fnal.gov>
 */

#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "MnvFluxConstraint.h"

using namespace MAT;

// static members must be initialized outside the declaration
bool MnvHistoConstrainer::s_disableTH1AddDir = true;

// ---------------------------------
void MnvHistoConstrainer::DisableAllConstraints()
{
	m_constraintsInUse.clear();
}

// ---------------------------------
void MnvHistoConstrainer::DisableConstraint(const std::string& constraintName) EXCEPTION_SPEC (ConstraintAccessError)
{
	m_constraintsInUse.erase(constraintName);
}

// ---------------------------------
void MnvHistoConstrainer::EnableAllConstraints()
{
	// this method is somewhat slow since it has to look through
	// all the entries in the master map.  but it shouldn't
	// be done often, so it's ok.
	for (std::map<UnivID, std::map<std::string, double> >::const_iterator it_univ = m_constraintWgts.begin();
		it_univ != m_constraintWgts.end();
		++it_univ)
	{
		for (std::map<std::string, double>::const_iterator it_constraint = it_univ->second.begin();
			it_constraint != it_univ->second.end();
			++it_constraint)
		{
			m_constraintsInUse.insert(it_constraint->first);
		} // for (it_constraint)
	} // for (it_univ)
}

// ---------------------------------
void MnvHistoConstrainer::EnableConstraint(const std::string& constraintName) EXCEPTION_SPEC (ConstraintAccessError)
{
	m_constraintsInUse.insert(constraintName);
}

// ---------------------------------
/// Is the error band 'errBandName' constrained by any of the constraints in the list 'constraintsToUse'?
bool MnvHistoConstrainer::ErrBandIsConstrained(const std::string& errBandName, const std::set<std::string>& constraintsToUse) const
{
	for (std::set<std::string>::const_iterator it_constraint = constraintsToUse.begin();
		it_constraint != constraintsToUse.end();
		++it_constraint)
	{
		std::pair<std::multimap<std::string, std::string>::const_iterator, std::multimap<std::string, std::string>::const_iterator> range = m_constraintErrBands.equal_range(*it_constraint);
		for (std::multimap<std::string, std::string>::const_iterator it_constrErr = range.first;
			it_constrErr != range.second;
			++it_constrErr)
		{
			if (it_constrErr->second == errBandName)
				return true;
		}
	}
	return false;
}

// ---------------------------------
MnvHistoConstrainer::SpectatorVarianceStrategy MnvHistoConstrainer::GetSpectatorCorrectionStrategy(const std::string & errName) const
	EXCEPTION_SPEC(MissingSpectatorStrategyError)
{
	if (m_spectatorErrBandStrategies.find(errName) == m_spectatorErrBandStrategies.end())
		throw MissingSpectatorStrategyError( ("Requested spectator error band '" + errName + "' has no strategy specified").c_str() );
	
	return m_spectatorErrBandStrategies.at(errName);
}

// ---------------------------------
const std::map<std::string, MnvHistoConstrainer::SpectatorVarianceStrategy> & MnvHistoConstrainer::GetSpectatorCorrectionStrategies() const
{
	return m_spectatorErrBandStrategies;
}

// ---------------------------------
void MnvHistoConstrainer::InsertUnivConstraint(const std::string& univGroupName, unsigned int univNum, const std::string & constraintName, double wgt)
	EXCEPTION_SPEC (ConstraintLoadError)
{
	UnivID univID(univGroupName, univNum);
	if (m_constraintWgts.find(univID) == m_constraintWgts.end())
		m_constraintWgts[univID] = std::map<std::string, double>();

	m_constraintErrBands.insert( std::make_pair(constraintName, univGroupName) );

	if (m_constraintWgts[univID].find(constraintName) != m_constraintWgts[univID].end())
	{
		std::stringstream ss;
		ss << "Universe '" << univGroupName << "' number " << univNum;
		ss << " already contains weight for constraint '" << constraintName << "': " << m_constraintWgts[univID][constraintName];
		throw ConstraintLoadError( ss.str().c_str() );
	}

	m_constraintWgts[univID][constraintName] = wgt;

	// add all loaded constraints to the working list by default
	m_constraintsInUse.insert(constraintName);
}

// ---------------------------------
/// Load a constraint from a file.
///
/// File format should be as follows (columns separated by whitespace):
///
/// univ_name univ_num weight
/// univ_name univ_num weight
///
/// etc.
void MnvHistoConstrainer::LoadConstraint(const std::string& constraintName, const std::string& constraintFilePath) EXCEPTION_SPEC (ConstraintLoadError)
{
	std::ifstream constraintFile(constraintFilePath.c_str());
	if (!constraintFile.is_open())
		throw ConstraintLoadError( ("Can't load constraint file: " + constraintFilePath).c_str() );

	std::string line;
	while (getline(constraintFile, line))
	{
		if (line.size() == 0)
			continue;

		std::stringstream ss(line);
		std::string univGroup;
		unsigned int univID;
		double wgt;
		ss >> univGroup;

		// make sure if file is garbled, user is informed
		if (ss.bad() || ss.fail())
			throw ConstraintLoadError( ("Can't parse text in file: " + constraintFilePath).c_str() );

		// deal with comments, which begin with '#'
		if (univGroup[0] == '#')
			continue;

		ss >> univID;
		ss >> wgt;

		// again make sure if file is garbled, user is informed
		// (note that we don't use ss.good() because it returns false
		//  if the end of the stream is reached, which SHOULD happen here)
		if (ss.bad() || ss.fail() || !ss.eof())
			throw ConstraintLoadError( ("Can't parse text in file: " + constraintFilePath).c_str() );

		this->InsertUnivConstraint(univGroup, univID, constraintName, wgt);
	}
} // MnvHistoConstrainer::LoadConstraint()

// ---------------------------------
void MnvHistoConstrainer::SetSpectatorCorrectionStrategy(const std::string & errName, MnvHistoConstrainer::SpectatorVarianceStrategy strategy)
	EXCEPTION_SPEC(SpectatorConstraintCollisionError)
{
	if ( this->ErrBandIsConstrained(errName, m_constraintsInUse) )
		throw SpectatorConstraintCollisionError( ("Cannot give error band '" + errName + "' a spectator correction strategy because it is already in use as a constraint").c_str() );
		
	m_spectatorErrBandStrategies[errName] = strategy;
}

// ---------------------------------
void MnvHistoConstrainer::SetSpectatorCorrectionStrategies(const std::map<std::string, MnvHistoConstrainer::SpectatorVarianceStrategy> & strategies)
	EXCEPTION_SPEC(SpectatorConstraintCollisionError)
{
	for (std::map<std::string, SpectatorVarianceStrategy>::const_iterator it_strat = strategies.begin();
		it_strat != strategies.end();
		++it_strat)
	{
		this->SetSpectatorCorrectionStrategy(it_strat->first, it_strat->second);
	}
}

