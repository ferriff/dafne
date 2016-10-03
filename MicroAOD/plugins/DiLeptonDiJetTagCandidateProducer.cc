#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"

#include "dafne/DataFormats/interface/DiLeptonDiJetTagCandidate.h"

#include <vector>
#include <string>

using namespace edm;
using namespace std;

namespace flashgg {

    class DiLeptonDiJetTagCandidateProducer : public EDProducer
    {

    public:
        DiLeptonDiJetTagCandidateProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        edm::EDGetTokenT<edm::OwnVector<flashgg::DiLeptonDiJetTagBase> > tagSorterToken_;

    };

    DiLeptonDiJetTagCandidateProducer::DiLeptonDiJetTagCandidateProducer( const ParameterSet &iConfig ) :
        tagSorterToken_( consumes<edm::OwnVector<flashgg::DiLeptonDiJetTagBase> >( iConfig.getUntrackedParameter<InputTag> ( "TagSorter")))
    {
        produces<std::vector<flashgg::DiLeptonDiJetTagCandidate> >();
    }

    void DiLeptonDiJetTagCandidateProducer::produce( Event &evt, const EventSetup & )
    {
        auto_ptr<std::vector<flashgg::DiLeptonDiJetTagCandidate> > tagsColl( new std::vector<flashgg::DiLeptonDiJetTagCandidate> );
        
        Handle<edm::OwnVector<flashgg::DiLeptonDiJetTagBase> > tagSorter;
        evt.getByToken( tagSorterToken_, tagSorter );

        const flashgg::DiLeptonDiJetTagBase *chosenTag = &*( tagSorter.product()->begin() );
        DiLeptonDiJetTagCandidate tags( chosenTag ); 
        tagsColl->push_back( tags );
 
        evt.put( tagsColl );                               
    }
}

typedef flashgg::DiLeptonDiJetTagCandidateProducer FlashggDiLeptonDiJetTagCandidateProducer;
DEFINE_FWK_MODULE( FlashggDiLeptonDiJetTagCandidateProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4