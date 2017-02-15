#ifndef FLASHgg_DiLeptonDiJetBase_h
#define FLASHgg_DiLeptonDiJetBase_h

#include "dafne/Systematics/interface/BaseSystMethod_DiLeptonDiJet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "flashgg/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <memory>

namespace flashgg {


    template <class param_var>
    class ElectronFromDiLeptonDiJetBase : public BaseSystMethod<DiLeptonDiJetCandidate, param_var>
    {

    public:
        ElectronFromDiLeptonDiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv );

        void applyCorrection( DiLeptonDiJetCandidate &y, param_var syst_shift ) override;
        float makeWeight( const DiLeptonDiJetCandidate &y, param_var syst_shift ) override;
        std::string shiftLabel( param_var ) const override;
        void eventInitialize( const edm::Event &iEvent, const edm::EventSetup & iSetup ) override;

        void setRandomEngine( CLHEP::HepRandomEngine &eng ) override
        {
            //            std::cout << " ElectronFromDiLeptonDiJetBase::setRandomEngine " << std::endl;
            BaseSystMethod<DiLeptonDiJetCandidate, param_var>::setRandomEngine( eng );
            electron_corr_->setRandomEngine( eng );
        }
        
    protected:
        bool debug_;

    private:
        std::unique_ptr<BaseSystMethod<flashgg::Electron, param_var> > electron_corr_;
        std::unique_ptr<BaseSystMethod<flashgg::Electron, param_var> > electron_corr2_;
    };

    template<class param_var>
    ElectronFromDiLeptonDiJetBase<param_var>::ElectronFromDiLeptonDiJetBase( const edm::ParameterSet &conf, edm::ConsumesCollector && iC, const GlobalVariablesComputer * gv ) :
        BaseSystMethod<DiLeptonDiJetCandidate, param_var>::BaseSystMethod( conf, std::forward<edm::ConsumesCollector>(iC) ),
        debug_( conf.getUntrackedParameter<bool>( "Debug", false ) )
    {
        std::string electronMethodName = conf.getParameter<std::string>( "ElectronMethodName" );
        electron_corr_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf, std::forward<edm::ConsumesCollector>(iC), gv ) );
        if(conf.exists("BinList2"))  //if defined, BinList2 gives bins for sublead, lead uses BinList
            {
                edm::ParameterSet conf2;// =  conf.clone();
                
                conf2.copyFrom(conf,"ElectronMethodName");
                conf2.copyFrom(conf,"MethodName");
                conf2.copyFrom(conf,"Label");
                conf2.copyFrom(conf,"NSigmas");
                conf2.copyFrom(conf,"OverallRange");
                conf2.copyFrom(conf,"Debug");
                conf2.copyFrom(conf,"ApplyCentralValue");
                const auto &pset = conf.getParameterSet( "BinList2" );
                conf2.addParameter<edm::ParameterSet>("BinList", pset);
                std::string binListName = "BinList";
                conf2.insertParameterSet(true,binListName, *(conf.retrieveUnknownParameterSet("BinList2")));
                electron_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf2, std::forward<edm::ConsumesCollector>(iC),  gv ) );
                
            }
        else { //if BinList2 is not defined, use BinList for both lead and sublead photons
            electron_corr2_.reset( FlashggSystematicMethodsFactory<flashgg::Electron, param_var>::get()->create( electronMethodName, conf, std::forward<edm::ConsumesCollector>(iC),  gv ) );
        }
        this->setMakesWeight( electron_corr_->makesWeight() );
    }

    template<class param_var>
    std::string ElectronFromDiLeptonDiJetBase<param_var>::shiftLabel( param_var syst_value ) const
    {
        return electron_corr_->shiftLabel( syst_value );
    }

    template<typename param_var>
    float ElectronFromDiLeptonDiJetBase<param_var>::makeWeight( const DiLeptonDiJetCandidate &y, param_var syst_shift )
    {
        if( debug_ ) {
            std::cout << "START OF ElectronFromDiLeptonDiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
                      // << y.mass() << " " << y.pt() << " " << y.leadingPhoton()->energy() << " " << y.subLeadingPhoton()->energy() << " "
                      // << y.leadingPhoton()->eta() << " " << y.subLeadingPhoton()->eta() 
                      << std::endl;
        }

        float weight1 = 1.;
        float weight2 = 1.;

        if (y.isEEJJ() || y.isEETT()) { 
            weight1 = electron_corr_->makeWeight( *(y.leadingEle()), syst_shift );
            weight2 = electron_corr2_->makeWeight( *(y.subLeadingEle()), syst_shift );
        }

        if (o.isEMJJ()) weight1 = electron_corr_->makeWeight( *(y.electron()), syst_shift );

        float dieleweight = weight1*weight2;
        if( debug_ ) {
            std::cout << "END OF ElectronFromDiLeptonDiJet::makeWeight M PT E1 E2 ETA1 ETA2 "
                      << " weight1=" << weight1 << " weight2=" << weight2 << " dieleweight=" << dieleweight << std::endl;
        }
        return dieleweight;
    }

    template<class param_var>
    void ElectronFromDiLeptonDiJetBase<param_var>::applyCorrection( DiLeptonDiJetCandidate &y, param_var syst_shift )
    {
        if( debug_ ) {
            std::cout << "START OF ElectronFromDiLeptonDiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
                      // << y.mass() << " " << y.pt() << " " << y.leadingPhoton()->energy() << " " << y.subLeadingPhoton()->energy() << " "
                      // << y.leadingPhoton()->eta() << " " << y.subLeadingPhoton()->eta() 
                      << std::endl;
        }
        // y.makePhotonsPersistent();
        
        if (y.isEEJJ() || y.isEETT()) {      
            electron_corr_->applyCorrection( *(y.leadingEle()), syst_shift );
            electron_corr2_->applyCorrection( *(y.subLeadingEle()), syst_shift );
        }

        if (o.isEMJJ()) electron_corr_->applyCorrection( *(y.electron()), syst_shift );

        y.computeP4AndOrder();
        if( debug_ ) {
            std::cout << "END OF ElectronFromDiLeptonDiJet::applyCorrection M PT E1 E2 ETA1 ETA2 "
                      // << y.mass() << " " << y.pt() << " " << y.leadingPhoton()->energy() << " " << y.subLeadingPhoton()->energy() << " "
                      // << y.leadingPhoton()->eta() << " " << y.subLeadingPhoton()->eta() 
                      << std::endl;
        }
    }

    template<class param_var>
    void ElectronFromDiLeptonDiJetBase<param_var>::eventInitialize( const edm::Event &ev, const edm::EventSetup & es )
    {
        if( debug_ ) {
            std::cout << "calling event initialize for both photons " << std::endl;
        }
        electron_corr_->eventInitialize( ev, es );
        electron_corr2_->eventInitialize( ev, es );
    }
}

#endif

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4