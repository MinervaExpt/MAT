#ifndef FluxReweighterWithWiggleFit_h
#define FluxReweighterWithWiggleFit_h

namespace PlotUtils {

    class FluxReweighter;
    class MnvH1D;
    
}

namespace PlotUtils {

    class FluxReweighterWithWiggleFit : public PlotUtils::FluxReweighter {
      public:
        FluxReweighterWithWiggleFit(int nuPDG,
                                    bool applyNuEConstraint,
                                    enum EPlaylist playlist,
                                    enum EFluxVersion fluxVersion,
                                    enum EG4NumiVersion g4NumiVersion);
        
        FluxReweighterWithWiggleFit(int nuPDG,
                                    bool applyNuEConstraint,
                                    std::string playlist_str,
                                    enum EFluxVersion fluxVersion,
                                    enum EG4NumiVersion g4NumiVersion);

         //! Get the central value flux weight.  Enu is true neutrino energy (GeV).
        virtual double GetFluxCVWeight( double Enu, int nuPDG );

        //! Get the flux error weights for all flux universes.  Enu is true neutrino energy (GeV).
        virtual std::vector<double> GetFluxErrorWeights( double Enu, int nuPDG );

        //! Get the flux error weight for a single flux universe.  Enu is true neutrino energy (GeV).
        virtual double GetFluxErrorWeight( double Enu, int nuPDG, unsigned int universe );

    
        MnvH1D* GetFluxWiggleMnvH1D(); //amits fit	


        virtual ~FluxReweighterWithWiggleFit();
        
      private:
        MnvH1D* MultiplyHists(MnvH1D* h1, MnvH1D* h2); // arguments actually are not modified

        MnvH1D* m_wiggleFitReweightHist;

    };

}
#endif
