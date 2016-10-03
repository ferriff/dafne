// DiLeptonDiJetTagTestAnalyzer.cc by S. Zenz
//
// * Tests getting tags out of the event, sorting them, casting them
// * Dumps debugging output to the screen
// * Useful for quick tests of code changes, and should be kept up-to-date as tags are added/changed
// * Should NOT be included in productions
//
// Adapted from globelikeTreeMakerWithTagSorter code by L. D. Corpe, which was
// Adapted from the flashggCommissioning tree maker code  by C. Favaro et al.

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "dafne/DataFormats/interface/DiLeptonDiJetTagBase.h"
#include "dafne/DataFormats/interface/DiEleDiJetTag.h"
#include "dafne/DataFormats/interface/DiMuDiJetTag.h"
#include "dafne/DataFormats/interface/DiEleDiTrackTag.h"
#include "dafne/DataFormats/interface/DiMuDiTrackTag.h"


using namespace std;
using namespace edm;


namespace flashgg {

    class DiLeptonDiJetTagTestAnalyzer : public edm::EDAnalyzer
    {
    public:
        explicit DiLeptonDiJetTagTestAnalyzer( const edm::ParameterSet & );
        ~DiLeptonDiJetTagTestAnalyzer();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    private:

        edm::Service<TFileService> fs_;

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::View<flashgg::DiLeptonDiJetTagBase> > TagSorterToken_;
        bool expectMultiples_;
    };

    DiLeptonDiJetTagTestAnalyzer::DiLeptonDiJetTagTestAnalyzer( const edm::ParameterSet &iConfig ):
        TagSorterToken_( consumes<edm::View<flashgg::DiLeptonDiJetTagBase> >( iConfig.getParameter<InputTag> ( "TagSorter" ) ) ),
        expectMultiples_( iConfig.getUntrackedParameter<bool>( "ExpectMultiples", false) )
    {
    }

    DiLeptonDiJetTagTestAnalyzer::~DiLeptonDiJetTagTestAnalyzer()
    {
    }

    void
    DiLeptonDiJetTagTestAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {
        Handle<edm::View<flashgg::DiLeptonDiJetTagBase> > TagSorter;
        iEvent.getByToken( TagSorterToken_, TagSorter );

        if (!expectMultiples_) {
            assert (TagSorter.product()->size() <= 1);
            if ( TagSorter.product()->size() == 0) std::cout << "[NO TAG]" << std::endl;
        } else {
            cout << "Multiple tags allowed and we have a total of " << TagSorter.product()->size() << endl;
        }

        for ( auto tag = TagSorter.product()->begin() ; tag != TagSorter.product()->end() ; tag++ ) {
            const flashgg::DiLeptonDiJetTagBase *chosenTag = &*( tag );

            const DiEleDiJetTag *dieledijet = dynamic_cast<const DiEleDiJetTag *>( chosenTag );
            if (dieledijet != NULL) {
                cout << "[DIELEDIJET] Category " << dieledijet->categoryNumber()
                << " with njets =" << dieledijet->jets().size()  
                //" jet pt=" << dieledijet->jet()->pt() << " dipho mass=" << zplusjet->diPhoton()->mass() 
                << " and nelectrons=" << dieledijet->electrons().size()
                << endl;
            }

            const DiMuDiJetTag *dimudijet = dynamic_cast<const DiMuDiJetTag *>( chosenTag );
            if (dimudijet != NULL) {
                cout << "[DIMUDIJET] Category " << dimudijet->categoryNumber()
                << " with njets =" << dimudijet->jets().size()  
                //" jet pt=" << dimudijet->jet()->pt() << " dipho mass=" << zplusjet->diPhoton()->mass() 
                << " and nmuons=" << dimudijet->muons().size()
                << endl;
            }

            const DiEleDiTrackTag *dieleditrack = dynamic_cast<const DiEleDiTrackTag *>( chosenTag );
            if (dieleditrack != NULL) {
                cout << "[DIELEDITRACK] Category " << dieleditrack->categoryNumber()
                << " with ntracks =" << dieleditrack->tracks().size()
                << " and nelectrons=" << dieleditrack->electrons().size()
                << endl;
            }

            const DiMuDiTrackTag *dimuditrack = dynamic_cast<const DiMuDiTrackTag *>( chosenTag );
            if (dimuditrack != NULL) {
                cout << "[DIMUDITRACK] Category " << dimuditrack->categoryNumber()
                << " with ntracks =" << dimuditrack->tracks().size() 
                << " and nmuons=" << dimuditrack->muons().size()
                << endl;
            }

            // << " with lead jet pt eta "
            //               << vbftag->leadingJet().pt() << " " << vbftag->leadingJet().eta()
            //               << " and sublead jet eta " << vbftag->subLeadingJet().pt() << " " << vbftag->subLeadingJet().eta() << " mass=" << vbftag->diPhoton()->mass()


            if( dieledijet == NULL && dimudijet == NULL && dieleditrack == NULL && dimuditrack == NULL ) {
                std::cout << "[FAILED TO CONVERT TAG] with SumPt " << chosenTag->sumPt() << std::endl;
            }

        } 
    } 

    void
    DiLeptonDiJetTagTestAnalyzer::beginJob()
    {
    }

    void
    DiLeptonDiJetTagTestAnalyzer::endJob()
    {
    }

    void
    DiLeptonDiJetTagTestAnalyzer::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

} 

typedef flashgg::DiLeptonDiJetTagTestAnalyzer FlashggDiLeptonDiJetTagTestAnalyzer;
DEFINE_FWK_MODULE( FlashggDiLeptonDiJetTagTestAnalyzer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4