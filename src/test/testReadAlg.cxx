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

#include <map>

/** @class testReadAlg
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/testReadAlg.cxx,v 1.2 2002/05/01 23:35:00 heather Exp $
 */

class testReadAlg : public Algorithm
{	
public:
    
    testReadAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:

    
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
            log << MSG::DEBUG << (*p)->fillStream(log.stream()) << endreq;
        }

    }

    SmartDataPtr<Event::McPositionHitVector> posHits(eventSvc(), EventModel::MC::McPositionHitCol);
    if (posHits) {
        log << MSG::DEBUG << "Retrieved McPositionHits from the TDS" << endreq;
        log << MSG::DEBUG << "Number of McPositionHits in the event = " << posHits->size() << endreq;
        Event::McPositionHitVector::const_iterator hit;
        for (hit = posHits->begin(); hit != posHits->end(); hit++ ) {   
            log << MSG::DEBUG << (*hit)->fillStream(log.stream()) << endreq;
        }
            
    }

    SmartDataPtr<Event::McIntegratingHitVector> intHits(eventSvc(), EventModel::MC::McIntegratingHitCol);
    if (intHits) {
        log << MSG::DEBUG << "Retrieved McIntegratingHits from the TDS" << endreq;
        log << MSG::DEBUG << "Number of McIntegratingHits in the event = " << intHits->size() << endreq;
        Event::McIntegratingHitVector::const_iterator hit;
        for (hit = intHits->begin(); hit != intHits->end(); hit++ ) {   
            log << MSG::DEBUG << (*hit)->fillStream(log.stream()) << endreq;
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


StatusCode testReadAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

