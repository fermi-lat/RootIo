#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/Digi/CalDigi.h"
#include "idents/CalXtalId.h"

#include <map>

/** @class testReadAlg
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/testReadAlg.cxx,v 1.3 2002/05/10 23:13:02 burnett Exp $
 */

class testReadAlg : public Algorithm
{	
public:
    
    testReadAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:

    StatusCode readMcData();

    StatusCode readDigiData();

    
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

    sc = readMcData();

    sc = readDigiData();

    return sc;
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
        idents::CalXtalId::CalTrigMode modeTds = (*calDigiTds)->getMode();
        idents::CalXtalId idTds = (*calDigiTds)->getPackedId();
        if (modeTds == idents::CalXtalId::BESTRANGE) {
            const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(0);
            log << MSG::DEBUG << "Using BestRange" << endreq;
        } else {
            log << MSG::DEBUG << "Using AllRange" << endreq;
            int range;
            for (range = idents::CalXtalId::LEX8; range <= idents::CalXtalId::HEX1; range++) {
                const Event::CalDigi::CalXtalReadout *readoutTds = (*calDigiTds)->getXtalReadout(range);
            }
        }
    }

    return sc;
}



StatusCode testReadAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

