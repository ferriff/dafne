#ifndef FLASHgg_DiLeptonDiJetTagBase_h
#define FLASHgg_DiLeptonDiJetTagBase_h

#include "dafne/DataFormats/interface/DiLeptonDiJetCandidate.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"

namespace flashgg {

    class DiLeptonDiJetTagBase : public WeightedObject
    {
    public:
        enum tag_t { kUndefined = 0, kDiEleDiJet, kDiMuDiJet, kDiEleDiTrack, kDiMuDiTrack };

        DiLeptonDiJetTagBase();
        virtual ~DiLeptonDiJetTagBase(); 
        DiLeptonDiJetTagBase( edm::Ptr<DiLeptonDiJetCandidate> );


        const edm::Ptr<DiLeptonDiJetCandidate> diLeptonDiJet() const { return diLeptonDiJet_; }

        int diLeptonDiJetIndex() const {return diLeptonDiJetIndex_;}
        void setDiLeptonDiJetIndex( int i ) { diLeptonDiJetIndex_ = i; }

        float sumPt() const { return this->diLeptonDiJet()->sumPt() ;}

        bool operator <( const DiLeptonDiJetTagBase &b ) const;

        virtual DiLeptonDiJetTagBase *clone() const { return ( new DiLeptonDiJetTagBase( *this ) ); }

        void setCategoryNumber( int value ) { category_number_ = value; }
        int categoryNumber() const { return category_number_; }
        operator int() const { return categoryNumber(); }

        void setTagTruth( const edm::Ptr<TagTruthBase> value ) { truth_ = value; }
        const edm::Ptr<TagTruthBase> tagTruth() const { return truth_; }

        void setSystLabel( const std::string label ) { systLabel_ = label; }
        std::string systLabel() const { return systLabel_; }
        bool hasSyst( const string &label ) const { return ( systLabel_ == label );}

        // void setIsGold ( int runNumber );
        // void setIsGoldMC( bool isGold ) { isGold_ = isGold; }
        // bool isGold() const { return isGold_; }

        virtual DiLeptonDiJetTagBase::tag_t tagEnum() const { return DiLeptonDiJetTagBase::kUndefined; }


        unsigned nOtherTags() const { 
            assert(otherTagTypes_.size() == otherTagCategories_.size());
            assert(otherTagTypes_.size() == otherTagIndices_.size());
            return otherTagTypes_.size(); 
        }

        void addOtherTag( const DiLeptonDiJetTagBase& other ) { 
            otherTagTypes_.push_back(other.tagEnum());
            otherTagCategories_.push_back(other.categoryNumber());
            otherTagIndices_.push_back(other.diLeptonDiJetIndex());
        }

        void addOtherTags( std::vector<std::tuple<DiLeptonDiJetTagBase::tag_t,int,int> > others ) { 
            for (unsigned i = 0 ; i < others.size() ; i++) {
                otherTagTypes_.push_back(std::get<0>(others[i]));
                otherTagCategories_.push_back(std::get<1>(others[i]));
                otherTagIndices_.push_back(std::get<2>(others[i]));
            }
        }

        DiLeptonDiJetTagBase::tag_t otherTagType( unsigned i ) const { return otherTagTypes_[i]; }

        int otherTagCategory( unsigned i ) const { return otherTagCategories_[i]; }

        int otherTagDiLeptonDiJetIndex ( unsigned i ) const { return otherTagIndices_[i]; }


    private:
        edm::Ptr<DiLeptonDiJetCandidate> diLeptonDiJet_;
        int diLeptonDiJetIndex_;
        int category_number_;
        edm::Ptr<TagTruthBase> truth_;
        string systLabel_;
        // bool isGold_;
        std::vector<std::tuple<DiLeptonDiJetTagBase::tag_t,int,int> > otherTags_; // (type,category,diphoton index)   
        std::vector<DiLeptonDiJetTagBase::tag_t> otherTagTypes_;
        std::vector<int> otherTagCategories_;
        std::vector<int> otherTagIndices_;
    };

}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

