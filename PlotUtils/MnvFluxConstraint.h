/*
 * MnvFluxConstraint.h
 *
 *  Applies flux constraint(s) to a MnvH?D histogram
 *  (including its universes).
 *
 *  Created on: Jul 14, 2014
 *      Author: J. Wolcott <jwolcott@fnal.gov>
 */

#ifndef MNVFLUXCONSTRAINT_H_
#define MNVFLUXCONSTRAINT_H_

//We can't have "dynamic exception specifications" in c++17.
//They really didn't do much before that as far as I can tell,
//but I'd rather not just remove them altogether.  So, I'll define
//a macro that turns into a dynamic exception specification if
//we're pre-c++-17 and noexcept(false) otherwise.
//
//Found the right variables to check for c++17 at https://stackoverflow.com/questions/47284705/c1z-dynamic-exception-specification-error
#if __cplusplus >= 201703L
  #define EXCEPTION_SPEC(TYPE) noexcept(false)
#else
  #define EXCEPTION_SPEC(TYPE) throw(TYPE)
#endif //__cplusplus >= 201703L

#include <map>
#include <string>
#include <set>

#include "MnvH1D.h"
#include "MnvH2D.h"
#include "MnvVertErrorBand.h"
#include "MnvVertErrorBand2D.h"
#include "Exceptions.h"

namespace PlotUtils
{
	// some forward declarations
//	class MnvH1D;
//	class MnvH2D;
//	class MnvVertErrorBand;
//	class MnvVertErrorBand2D;
	class MnvVertErrorBand3D;


	/*!
	 * @brief A class which can apply a flux constraint to an MnvH1D, MnvH2D
	 * @author J. Wolcott <jwolcott@fnal.gov>
	 * For more information, see http://minerva-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=10281
	 *
	 * \warning This class sets TH1::AddDirectory(false) because it does a lot of copying of histograms
	 * and is painfully slow otherwise.  If you depend on retrieving histograms from the current directory
	 * from within the macro in which an MnvHistoConstrainer is being used, you may use
	 * the static MnvHistoConstrainer::DisableTH1AddDirectory() before instantiating an MnvHistoConstrainer.
	 */
	class MnvHistoConstrainer
	{
		public:
			/// custom container for universe IDs which can be used as a std::map key
			struct UnivID
			{
				UnivID(const std::string & univGroup, unsigned int univID) :
					groupName(univGroup), id(univID)
				{};

				bool operator<(const UnivID& other) const { return groupName < other.groupName || (groupName == other.groupName && id < other.id); };

				std::string groupName;
				unsigned int id;
			};
			
			MnvHistoConstrainer() { if (s_disableTH1AddDir) TH1::AddDirectory(false); };
			virtual ~MnvHistoConstrainer() {};

			/// enumeration specifying which strategy to use when adjusting 'spectator' error bands
			enum SpectatorVarianceStrategy { PRESERVE_ABSOLUTE_ERR, PRESERVE_FRACTIONAL_ERR };

			///@{
			/// @name public interface

			/// workhorse method: given an MnvH?D, returns a NEW MnvH?D filled
			/// with the result of applying the constraint
			template<class MnvHistoType, class MnvVertErrBandType>
			MnvHistoType * ConstrainHisto(const MnvHistoType * inHisto) const;

			/// workhorse method with option to specify which constraint(s) to apply
			template<class MnvHistoType, class MnvVertErrBandType>
			MnvHistoType * ConstrainHisto(const MnvHistoType * inHisto, const std::set<std::string> & constraintsToUse) const;

			///@{
			/// @name Constraint management
			void LoadConstraint(const std::string & constraintName, const std::string & constraintFilePath) EXCEPTION_SPEC (ConstraintLoadError);

			void DisableAllConstraints();
			void DisableConstraint(const std::string & constraintName) EXCEPTION_SPEC (ConstraintAccessError);
			void EnableAllConstraints();
			void EnableConstraint(const std::string & constraintName) EXCEPTION_SPEC (ConstraintAccessError);
			///@}
			///@}
			
			///@{
			/// @name Spectator error strategy management
			
			/// get a single strategy for a specific error band
			SpectatorVarianceStrategy GetSpectatorCorrectionStrategy(const std::string & errName) const EXCEPTION_SPEC(MissingSpectatorStrategyError);
			
			/// retrieve all the strategies this object knows about
			const std::map<std::string, SpectatorVarianceStrategy> & GetSpectatorCorrectionStrategies() const;
			
			/// set the strategy for one spectator error band.
			/// note that if this error name is in the constraint list, an exception will be thrown
			/// (can't be a spectator and a constraint simultaneously).
			void SetSpectatorCorrectionStrategy(const std::string & errName, SpectatorVarianceStrategy strategy) EXCEPTION_SPEC(SpectatorConstraintCollisionError);
			
			/// set the strategies for multiple spectator error bands.
			void SetSpectatorCorrectionStrategies(const std::map<std::string, SpectatorVarianceStrategy> & strategies) EXCEPTION_SPEC(SpectatorConstraintCollisionError);
			
			///@}

			/// Use this method to specify whether MnvHistoConstrainer should use TH1::AddDirectory(false)
			/// to disable the automatic adding of histograms to the current (ROOT) gDirectory.
			/// Don't use it unless you need to -- performance will suffer substantially.
			static void DisableTH1AddDirectory(bool doDisable) { s_disableTH1AddDir = doDisable; };

		private:
			///@{
			/// @name algorithm internals
			void InsertUnivConstraint(const std::string & univGroupName, unsigned int univID, const std::string & constraintName, double wgt) EXCEPTION_SPEC (ConstraintLoadError);
			
			template <class MnvHistoType, class MnvVertErrBandType>
			void CalcErrBandMean(const MnvHistoType * inHisto, MnvHistoType * outHisto, const std::set<std::string> & constraintsToUse, bool weighted) const;

			template <class MnvHistoType, class MnvVertErrBandType>
			void CorrectCV(const MnvHistoType * inHisto, MnvHistoType* outHisto, const MnvHistoType * correctedAvgs, const std::set<std::string> & constraintsToUse) const;

			template <class MnvHistoType, class MnvVertErrBandType>
			void CorrectFluxUniv(MnvHistoType * histoToCorrect, const std::set<std::string>& constraintsToUse) const;

			template<class MnvHistoType>
			void CorrectSpectatorUniv(const MnvHistoType * inHisto, MnvHistoType* outHisto, const std::set<std::string>& constraintsToUse) const;

			bool ErrBandIsConstrained(const std::string& errBandName, const std::set<std::string>& constraintsToUse) const;
			///@}

			///@{
			/// @name data members

			/// the names of the error bands adjusted by each of the constraints (filled automatically while loading constraints)
			std::multimap<std::string, std::string> m_constraintErrBands;

			/// which constraints will we apply by default to a histogram?
			std::set<std::string> m_constraintsInUse;

			/// the weights for each of the constraints, organized by universe
			std::map<UnivID, std::map<std::string, double> > m_constraintWgts;
			
			/// what aspect of  should be preserved when adjusting 'spectator' error bands: absolute or fractional variance?
			std::map<std::string, SpectatorVarianceStrategy> m_spectatorErrBandStrategies;

			/// should I disable default "add TH1 to current directory" behavior of ROOT?  (default: true)
			static bool s_disableTH1AddDir;
			///@}
	};  // class MnvHistoConstrainer

	// need to explicitly instantiate the template for the cases we'll use in CINT/PyROOT
	// (that way the dictionary generation is complete)
	template MnvH1D * MnvHistoConstrainer::ConstrainHisto<MnvH1D, MnvVertErrorBand>(const MnvH1D* inHisto) const;
	template MnvH2D * MnvHistoConstrainer::ConstrainHisto<MnvH2D, MnvVertErrorBand2D>(const MnvH2D* inHisto) const;

	// /////////////////////////////////////////////////////////////////////////////////////////
	//
	//     templated methods
	//     (must live in the header because they need to be explicitly instantiated, per above)
	//
	// /////////////////////////////////////////////////////////////////////////////////////////

	/// Compute the mean of the histograms in an error band using the constraint values.
	template <class MnvHistoType, class MnvVertErrBandType>
	void MnvHistoConstrainer::CalcErrBandMean(const MnvHistoType * inHisto, MnvHistoType * outHisto, const std::set<std::string> & constraintsToUse, bool weighted) const
	{
		// start fresh.
		outHisto->ClearAllErrorBands();
		outHisto->Reset();
	
		double totalUnivWgtSum = 0;

		std::vector<std::string> vertErrorBandNames(inHisto->GetVertErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = vertErrorBandNames.begin();
			it_eb != vertErrorBandNames.end();
			++it_eb)
		{
			// we only want to consider the flux-type errors constrained by this constraint in this new histogram
			if (!this->ErrBandIsConstrained(*it_eb, constraintsToUse))
				continue;
				
			const MnvVertErrBandType * errBand = inHisto->GetVertErrorBand(*it_eb);
			for (unsigned int univ_idx = 0; univ_idx < errBand->GetNHists(); ++univ_idx)
			{
				UnivID univID(*it_eb, univ_idx);
				if (m_constraintWgts.find(univID) == m_constraintWgts.end())
					continue;

				// this is just a working copy of this universe within the error band
				MnvHistoType errBandUniv(*errBand->GetHist(univ_idx));

				double univConstraint = 1;
				if (weighted)
				{
					const std::map<std::string, double> & constraints((m_constraintWgts.find(univID))->second);
					for (std::set<std::string>::const_iterator it_constraint = constraintsToUse.begin();
						it_constraint != constraintsToUse.end();
						++it_constraint)
					{
						if (constraints.find(*it_constraint) == constraints.end())
							continue;

						univConstraint *= constraints.find(*it_constraint)->second;
					}
				}

				if (univConstraint != 1)
					errBandUniv.Scale(univConstraint);
				outHisto->Add(&errBandUniv);
				totalUnivWgtSum += univConstraint;
			} // for (univ_idx)
		} // for (it_eb)

		if (totalUnivWgtSum == 0)
			throw NoFluxUnivError("No universes were selected by this constraint in calculating the mean!");

		outHisto->Scale(1./totalUnivWgtSum);
	} // MnvHistoConstrainer::CalcErrBandMean

	/// This variant just uses the default "constraints to use" in m_constraintsInUse
	/// and passes it along to the real version
	template<class MnvHistoType, class MnvVertErrBandType>
	MnvHistoType* MnvHistoConstrainer::ConstrainHisto(const MnvHistoType* inHisto) const
	{
		return this->ConstrainHisto<MnvHistoType, MnvVertErrBandType>(inHisto, m_constraintsInUse);
	}

	/// Apply the requested flux constraints to a histogram.
	template<class MnvHistoType, class MnvVertErrBandType>
	MnvHistoType* MnvHistoConstrainer::ConstrainHisto(const MnvHistoType* inHisto, const std::set<std::string> & constraintsToUse) const
	{
		// you might think it's inefficient to deep-copy all the error bands
		// just to clear them a moment later, and you'd be right.
		// but this enables me to copy the input histogram's bin structure 
		// and not have to include any more template arguments to delete the other error bands
		// (some of which are different types like MnvLatErrorBand/MnvLatErrorBand{2,3}D) one-by-one.
		// we'll need these corrections for both the CV and flux bands, so store it.
		MnvHistoType * weightedMeans = new MnvHistoType( inHisto->GetCVHistoWithStatError() );  // just a working container
		weightedMeans->ClearAllErrorBands();
		weightedMeans->Reset();
		
		MnvHistoType * newHisto = new MnvHistoType(*inHisto);
		
		this->CalcErrBandMean<MnvHistoType, MnvVertErrBandType>(inHisto, weightedMeans, constraintsToUse, true);  // a weighted mean
		
		// calculate the corrections...
		this->CorrectFluxUniv<MnvHistoType, MnvVertErrBandType>(newHisto, constraintsToUse);
		this->CorrectCV<MnvHistoType, MnvVertErrBandType>(inHisto, newHisto, weightedMeans, constraintsToUse);
		this->CorrectSpectatorUniv<MnvHistoType>(inHisto, newHisto, constraintsToUse);

		// clean up
		delete weightedMeans;

		// finally, send back the answer
		return newHisto;
	}

	/// @brief Determine the correction factors in each bin for the CV
	/// (and non-flux-type error bands) based on the constraints in use,
	/// then apply them to the given histogram.
	///
	/// The corrections are calculated based on the method outlined
	/// in DocDB 10076, p. 6.
	template<class MnvHistoType, class MnvVertErrBandType>
	void MnvHistoConstrainer::CorrectCV(const MnvHistoType * inHisto, MnvHistoType* outHisto, const MnvHistoType * correctedAvgs, const std::set<std::string> & constraintsToUse) const
	{
		// first clone the histogram.  we'll be returning a new one.
		MnvHistoType * newHisto = new MnvHistoType(*inHisto);
		newHisto->SetName( (std::string(inHisto->GetName()) + "_fluxconstrained").c_str() );

		// this histogram will be used to store the unconstrained averages,
		// which we'll compare to the constrained ones to get the corrections.
		// (it doesn't need a name change since it is internal only.)
		MnvHistoType * uncorrectedAvgs = new MnvHistoType(*inHisto);
		uncorrectedAvgs->ClearAllErrorBands();
		uncorrectedAvgs->Reset();
		this->CalcErrBandMean<MnvHistoType, MnvVertErrBandType>(inHisto, uncorrectedAvgs, constraintsToUse, false); // the unweighted mean.  here we throwing away the return value since it is just 1

		// copy the corrected averages histogram since we'll be working with it further
		MnvHistoType * correctionHisto = new MnvHistoType(*correctedAvgs);
		for (int bin_idx = 0; bin_idx < correctionHisto->fN; bin_idx++ )
			correctionHisto->SetBinError(bin_idx, 0);

		// the correction factors are the ratio of these two.
		correctionHisto->Divide(correctedAvgs, uncorrectedAvgs);

		// now, ensure that all the error bands contain unity so that the multiplication
		// does the correct thing
		correctionHisto->AddMissingErrorBandsAndFillWithCV(*inHisto);
		std::vector<std::string> vertErrorBandNames(inHisto->GetVertErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = vertErrorBandNames.begin();
			it_eb != vertErrorBandNames.end();
			++it_eb)
		{
			MnvVertErrBandType * errBand = correctionHisto->GetVertErrorBand(*it_eb);
			for (unsigned int univ_idx = 0; univ_idx < errBand->GetNHists(); ++univ_idx)
			{
				errBand->GetHist(univ_idx)->Reset();
				// this is delightfully ugly.
				// to avoid having to include ANOTHER template parameter
				// (i.e., the type of ROOT histogram in the error bands),
				// we take advantage of the fact that ROOT's histograms
				// are all linearized arrays at heart, and that the arrays
				// are publicly accessible.
				for (int bin_idx = 0; bin_idx < errBand->GetHist(univ_idx)->fN; bin_idx++ )
					errBand->GetHist(univ_idx)->fArray[bin_idx] = 1.0;
			}
		}
		std::vector<std::string> latErrorBandNames(inHisto->GetLatErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = latErrorBandNames.begin();
			it_eb != latErrorBandNames.end();
			++it_eb)
		{
			for (unsigned int univ_idx = 0; univ_idx < correctionHisto->GetLatErrorBand(*it_eb)->GetNHists(); ++univ_idx)
			{
				correctionHisto->GetLatErrorBand(*it_eb)->GetHist(univ_idx)->Reset();
				// same trick as above
				for (int bin_idx = 0; bin_idx < correctionHisto->GetLatErrorBand(*it_eb)->GetHist(univ_idx)->fN; bin_idx++ )
					correctionHisto->GetLatErrorBand(*it_eb)->GetHist(univ_idx)->fArray[bin_idx] = 1.0;
			}
		}
		
		// apply the constraint (to the CV only since the error bands all contain unity).
		outHisto->Multiply(outHisto, correctionHisto);
		
		// clean up
		delete uncorrectedAvgs;
		delete correctionHisto;

	} // MnvHistoConstrainer::CorrectCV()

	/// @brief Add the consolidated weights from the applicable constraints to the appropriate error bands' weight lists.
	///
	/// Since the error bands themselves know how to correctly calculate weighted covariance matrices,
	/// we simply add the weights to their lists.
	template<class MnvHistoType, class MnvVertErrBandType>
	void MnvHistoConstrainer::CorrectFluxUniv(MnvHistoType * histoToCorrect, const std::set<std::string>& constraintsToUse) const
	{
		std::vector<std::string> vertErrorBandNames(histoToCorrect->GetVertErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = vertErrorBandNames.begin();
			it_eb != vertErrorBandNames.end();
			++it_eb)
		{
			if (! this->ErrBandIsConstrained(*it_eb, constraintsToUse))
				continue;
				
//			std::cout << " Constraining error band: " << *it_eb << std::endl;

			MnvVertErrBandType * errBand = histoToCorrect->GetVertErrorBand(*it_eb);
			for (unsigned int univ_idx = 0; univ_idx < errBand->GetNHists(); ++univ_idx)
			{
				UnivID univID(*it_eb, univ_idx);
				if (m_constraintWgts.find(univID) == m_constraintWgts.end())
					continue;

				double univConstraint = 1;
				const std::map<std::string, double> & constraints((m_constraintWgts.find(univID))->second);
				for (std::set<std::string>::const_iterator it_constraint = constraintsToUse.begin();
					it_constraint != constraintsToUse.end();
					++it_constraint)
				{
					if (constraints.find(*it_constraint) == constraints.end())
						continue;

					univConstraint *= constraints.find(*it_constraint)->second;
				}

//				std::cout<< "  applying weight: " << univConstraint << " to universe: " << univ_idx << std::endl;
				errBand->SetUnivWgt(univ_idx, univConstraint);
			} // for (univ_idx)
		} // for (it_eb)

	} // MnvHistoConstrainer::CorrectFluxUniv()
	
	template<class MnvHistoType>
	void MnvHistoConstrainer::CorrectSpectatorUniv(const MnvHistoType * inHisto, MnvHistoType* outHisto, const std::set<std::string>& constraintsToUse) const
	{
		// first, determine the effect that the correction had on the CV.
		// we need to store this two ways, because depending on the strategy
		// of the error band, it needs to be applied differently:
		//  - for error bands using the PRESERVE_ABSOLUTE_ERR strategy, we'll be ADDING the difference;
		//  - for error bands using the PRESERVE_FRACTIONAL_ERR strategy, we'll be MULTIPLYING by the quotient.
		MnvHistoType * differenceHisto = new MnvHistoType(*outHisto);
		differenceHisto->Add(inHisto, -1);
		differenceHisto->ClearAllErrorBands();

		MnvHistoType * quotientHisto = new MnvHistoType(*outHisto);
		quotientHisto->Divide(outHisto, inHisto);
		quotientHisto->ClearAllErrorBands();
		
		// we only want to use the difference on the error bands,
		// since the CV is already corrected.  propagate the difference
		// to the errors, then set the CV of the 'difference' histogram to have identity
		// (0 for the differenceHisto and 1 for the quotientHisto).
		differenceHisto->AddMissingErrorBandsAndFillWithCV(*inHisto);
		quotientHisto->AddMissingErrorBandsAndFillWithCV(*inHisto);
		for (int bin_idx = 0; bin_idx < differenceHisto->fN; bin_idx++ )
		{
			differenceHisto->SetBinContent(bin_idx, 0);
			differenceHisto->SetBinError(bin_idx, 0);
			
			quotientHisto->SetBinContent(bin_idx, 1);
			quotientHisto->SetBinError(bin_idx, 0);
		}

		// remember, each error band has its own copy of the CV too........
		// also, make sure the 'difference' and 'quotient' histograms
		// only contain corrections for the error bands they're actually meant to correct.
		std::vector<std::string> latErrorBandNames(inHisto->GetLatErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = latErrorBandNames.begin();
			it_eb != latErrorBandNames.end();
			++it_eb)
		{
			if ( this->ErrBandIsConstrained(*it_eb, constraintsToUse) )
			{
				delete quotientHisto->PopLatErrorBand(*it_eb);
				delete differenceHisto->PopLatErrorBand(*it_eb);
				continue;
			}

			if ( m_spectatorErrBandStrategies.find(*it_eb) == m_spectatorErrBandStrategies.end())
				throw MissingSpectatorStrategyError( ("Error strategy for spectator lateral error '" + (*it_eb) + "' was not specified").c_str() );


			MnvHistoType * histoToAdjust;
			double identity;
			if ( m_spectatorErrBandStrategies.at(*it_eb) == PRESERVE_ABSOLUTE_ERR)
			{
				delete quotientHisto->PopLatErrorBand(*it_eb);
				histoToAdjust = differenceHisto;
				identity = 0;
			}
			else
			{
				delete differenceHisto->PopLatErrorBand(*it_eb);
				histoToAdjust = quotientHisto;
				identity = 1;
			}

			for (int bin_idx = 0; bin_idx < outHisto->GetLatErrorBand(*it_eb)->fN; bin_idx++)
				histoToAdjust->GetLatErrorBand(*it_eb)->fArray[bin_idx] = identity;
		}
		std::vector<std::string> vertErrorBandNames(inHisto->GetVertErrorBandNames());
		for (std::vector<std::string>::const_iterator it_eb = vertErrorBandNames.begin();
			it_eb != vertErrorBandNames.end();
			++it_eb)
		{
			if ( this->ErrBandIsConstrained(*it_eb, constraintsToUse) )
			{
				delete quotientHisto->PopVertErrorBand(*it_eb);
				delete differenceHisto->PopVertErrorBand(*it_eb);
				continue;
			}

			if ( m_spectatorErrBandStrategies.find(*it_eb) == m_spectatorErrBandStrategies.end())
				throw MissingSpectatorStrategyError( ("Error strategy for spectator vertical error '" + (*it_eb) + "' was not specified").c_str() );

			MnvHistoType * histoToAdjust;
			double identity;
			if ( m_spectatorErrBandStrategies.at(*it_eb) == PRESERVE_ABSOLUTE_ERR)
			{
				delete quotientHisto->PopVertErrorBand(*it_eb);
				histoToAdjust = differenceHisto;
				identity = 0;
			}
			else
			{
				delete differenceHisto->PopVertErrorBand(*it_eb);
				histoToAdjust = quotientHisto;
				identity = 1;
			}

			for (int bin_idx = 0; bin_idx < outHisto->GetVertErrorBand(*it_eb)->fN; bin_idx++)
				histoToAdjust->GetVertErrorBand(*it_eb)->fArray[bin_idx] = identity;
		}

		// once more, to fill in any error bands that were removed above.
		// this way they won't do anything in the correction step below
		// (since at this point the CV contains identity, and thus identity
		//  will be propagated to them).
		differenceHisto->AddMissingErrorBandsAndFillWithCV(*inHisto);
		quotientHisto->AddMissingErrorBandsAndFillWithCV(*inHisto);

//		differenceHisto->Print("all");
		
		// finally, correct the error bands by adding the difference.
		outHisto->Multiply(outHisto, quotientHisto);
		outHisto->Add(differenceHisto);

		delete differenceHisto;
		delete quotientHisto;
	} // MnvHistoConstrainer::CorrectSpectatorUniv()

}   // namespace PlotUtils


#endif /* MNVFLUXCONSTRAINT_H_ */
