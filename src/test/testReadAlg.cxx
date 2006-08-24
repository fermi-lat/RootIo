#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/Digi/CalDigi.h"
#include "idents/CalXtalId.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalCluster.h"   
#include "Event/Recon/CalRecon/CalXtalRecData.h"   

#include "LdfEvent/Gem.h"
#include "LdfEvent/LsfMetaEvent.h"
#include "LdfEvent/LsfCcsds.h"

#include "OnboardFilterTds/FilterStatus.h"

#include <map>

/** @class testReadAlg
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/testReadAlg.cxx,v 1.16 2006/06/23 07:17:34 heather Exp $
 */

class testReadAlg : public Algorithm
{	
public:
    
    testReadAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:
    StatusCode readEvtHeader();

    StatusCode readLsfData();

    StatusCode readMcData();

    StatusCode readDigiData();
    StatusCode readGemData();

    StatusCode readReconData();
    
    StatusCode readOnboardFilter();
};

static const AlgFactory<testReadAlg>  Factory;
const IAlgFactory& testReadAlgFactory = Factory;


testReadAlg::testReadAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode testReadAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();
    
    return sc;
    
}

StatusCode testReadAlg::execute()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    sc = readEvtHeader();

    sc = readLsfData();

    sc = readMcData();

    sc = readDigiData();
    sc = readGemData();

    sc = readReconData();

    sc = readOnboardFilter();

    return sc;
}
StatusCode testReadAlg::readEvtHeader() {
    StatusCode sc = StatusCode::SUCCESS;
    SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);

    if (!evtTds) return sc;

    int evtId = evtTds->event();
    int runId = evtTds->run();

    SmartDataPtr<Event::MCEvent> mcEvt(eventSvc(), EventModel::MC::Event);
    if (!mcEvt) return sc;

    int sourceid = mcEvt->getSourceId();
    int seq = mcEvt->getSequence();
    double t = mcEvt->time();

    return sc;
}

StatusCode testReadAlg::readLsfData() {
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    SmartDataPtr<LsfEvent::MetaEvent> metaTds(eventSvc(), "/Event/MetaEvent");
    if (!metaTds)
       log << MSG::DEBUG << "No MetaEvent data" << endreq;
    else {
        log << MSG::DEBUG;
        metaTds->fillStream(log.stream());
        log << endreq;
     }

    SmartDataPtr<LsfEvent::LsfCcsds> ccsdsTds(eventSvc(), "/Event/Ccsds");
    if (!ccsdsTds)
       log << MSG::DEBUG << "No CCSDS data" << endreq;
    else {
        log << MSG::DEBUG;
        ccsdsTds->fillStream(log.stream());
        log << endreq;
     }

	return StatusCode::SUCCESS;
}

StatusCode testReadAlg::readMcData() {
    // Purpose and Method: Retrieve Monte Carlo data from the TDS.  This data 
    //    was put on the the TDS by the mcRootReaderAlg.
    // TDS Input:  EventModel::MC::McParticleCol, 
    //    EventModel::MC::McPositionHitCol, EventModel::MC::McIntegratingHitCol

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::McParticleCol> particles(eventSvc(), EventModel::MC::McParticleCol);
    if (particles) {
        log << MSG::DEBUG << "Retrieved McParticles from the TDS" << endreq;
        log << MSG::DEBUG << "Number of McParticles in the event = " << particles->size() << endreq;
        Event::McParticleCol::const_iterator p;
        for (p = particles->begin(); p != particles->end(); p++) {
            log << MSG::DEBUG;
            (*p)->fillStream(log.stream());
            log << endreq;
            log << MSG::DEBUG << "Number of daughters " << (*p)->daughterList().size() << endreq;
        }

    }

    SmartDataPtr<Event::McPositionHitVector> posHits(eventSvc(), EventModel::MC::McPositionHitCol);
    if (posHits) {
        log << MSG::DEBUG << "Retrieved McPositionHits from the TDS" << endreq;
        log << MSG::DEBUG << "Number of McPositionHits in the event = " << posHits->size() << endreq;
        Event::McPositionHitVector::const_iterator hit;
        for (hit = posHits->begin(); hit != posHits->end(); hit++ ) {   
            log << MSG::DEBUG;
            (*hit)->fillStream(log.stream());
            log << endreq;
        }
            
    }

    SmartDataPtr<Event::McIntegratingHitVector> intHits(eventSvc(), EventModel::MC::McIntegratingHitCol);
    if (intHits) {
        log << MSG::DEBUG << "Retrieved McIntegratingHits from the TDS" << endreq;
        log << MSG::DEBUG << "Number of McIntegratingHits in the event = " << intHits->size() << endreq;
        Event::McIntegratingHitVector::const_iterator hit;
        for (hit = intHits->begin(); hit != intHits->end(); hit++ ) {   
            log << MSG::DEBUG;
            (*hit)->fillStream(log.stream());
            log << endreq;
            Event::McIntegratingHit::energyDepositMap mcPartMap = (*hit)->itemizedEnergy();
            Event::McIntegratingHit::energyDepositMap::const_iterator partIt;
            log << MSG::DEBUG << "McIntegratingHit energy Map" << endreq;
            for (partIt = mcPartMap.begin(); partIt != mcPartMap.end(); partIt++) {
                log << MSG::DEBUG << "(McPartId, energy) = (" << partIt->first->particleProperty()
                    << ", " << partIt->second << ")" << endreq;
            }
        }
    }

    return sc;
}


StatusCode testReadAlg::readDigiData() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    SmartDataPtr<Event::CalDigiCol> calDigiColTds(eventSvc(), EventModel::Digi::CalDigiCol);
    if (!calDigiColTds) return sc;
    Event::CalDigiCol::const_iterator calDigiTds;

    for (calDigiTds = calDigiColTds->begin(); calDigiTds != calDigiColTds->end(); calDigiTds++) {
        log << MSG::DEBUG;
        (*calDigiTds)->fillStream(log.stream());
        log << endreq;
        idents::CalXtalId::CalTrigMode modeTds = (*calDigiTds)->getMode();
        idents::CalXtalId idTds = (*calDigiTds)->getPackedId();
        if (modeTds == idents::CalXtalId::BESTRANGE) {
            const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(0);
        } else {
            int range;
            for (range = idents::CalXtalId::LEX8; range <= idents::CalXtalId::HEX1; range++) {
                const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(range);
            }
        }
    }

    return sc;
}

StatusCode testReadAlg::readGemData() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<LdfEvent::Gem> gemTds(eventSvc(), "/Event/Gem");
    if (!gemTds) {
        log << MSG::INFO << "No GEM available" << endreq;
    } else {
       log << MSG::INFO;
       (*gemTds).fillStream(log.stream());
       log << endreq;
    }

    return sc;
}

StatusCode testReadAlg::readReconData() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<Event::TkrTrackCol> trackColTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
    if (!trackColTds)
        log << MSG::INFO << "No TKR track collection on the TDS" << endreq;
    else
        log << MSG::DEBUG << trackColTds->size() << " tracks in TDS" << endreq;

    SmartDataPtr<Event::TkrVertexCol> vertexColTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    if (!vertexColTds) 
        log << MSG::INFO << "No TKR vertex collection on TDS" << endreq;
    else 
        log << MSG::DEBUG << vertexColTds->size() << " TKR vertices on TDS" << endreq;


    SmartDataPtr<Event::CalXtalRecCol> xtalRecColTds(eventSvc(),EventModel::CalRecon::CalXtalRecCol);
    if (!xtalRecColTds)
        log << MSG::INFO << "No CAL recon xtal collection on the TDS" << endreq;
    else 
        log << MSG::DEBUG << xtalRecColTds->size() << " CAL recon xtals on TDS" << endreq;

    SmartDataPtr<Event::CalClusterCol> calClusterColTds(eventSvc(),EventModel::CalRecon::CalClusterCol);
    if (!calClusterColTds)
        log << MSG::INFO << "No CAL recon cluster collection on the TDS" << endreq;
    else 
        log << MSG::DEBUG << calClusterColTds->size() << " CAL clusters on TDS" << endreq;

    SmartDataPtr<Event::AcdRecon> acdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (!acdRecTds) 
        log << MSG::INFO << "No ACD recon data on TDS" << endreq;

    return sc;
}


StatusCode testReadAlg::readOnboardFilter() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    SmartDataPtr<OnboardFilterTds::FilterStatus> obfTds(eventSvc(), "/Event/Filter/FilterStatus"); 

    if (!obfTds) {
        log << MSG::INFO << "No OBF on TDS" << endreq;
        return StatusCode::SUCCESS;
    }

    log << MSG::DEBUG << "OnboardFilter: " << endreq;
    log << MSG::DEBUG << "Status: " << obfTds->get() << " StageEnergy: "
        << obfTds->getStageEnergy() << " TCIDS: " 
        << obfTds->getTcids() << endreq;
  
    log << MSG::DEBUG << "  GEM: " << endreq;
    log << " ThrTkr: " << obfTds->getGemThrTkr() << " CalHiLo: " 
        << obfTds->getGemCalHiLo() << " Condsumcno: " 
        << obfTds->getGemCondsumCno() << endreq;

    log << MSG::DEBUG << "Separation: " << obfTds->getSeparation() << endreq;

    return sc;

}

StatusCode testReadAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

