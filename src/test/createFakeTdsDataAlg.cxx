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
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/AcdRecon/AcdRecon.h"

#include "LdfEvent/Gem.h"

#include "idents/CalXtalId.h"
#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

/** @class creatFakeTdsDataAlg
 * @brief Creates dummy data to put on the TDS without requiring any
 * Monte Carlo generator or running any of the standard algorithms.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/createFakeTdsDataAlg.cxx,v 1.5.6.2 2004/08/18 20:37:46 heather Exp $
 */

class createFakeTdsDataAlg : public Algorithm
{	
public:
    
    createFakeTdsDataAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:

    StatusCode storeMcData();

    StatusCode storeDigiData();

    StatusCode storeReconData();

    StatusCode storeGemData();

    unsigned int m_numMcParticles, m_numMcPosHits, m_numMcIntHits;
    unsigned int m_numAcdDigi, m_numCalDigi, m_numTkrDigi;
    unsigned int m_numCalXtals, m_numTkrTracks;
    
};

static const AlgFactory<createFakeTdsDataAlg>  Factory;
const IAlgFactory& createFakeTdsDataAlgFactory = Factory;


createFakeTdsDataAlg::createFakeTdsDataAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{

	declareProperty("numMcParticles", m_numMcParticles=20);
    declareProperty("numMcPosHits", m_numMcPosHits=15);
    declareProperty("numMcIntHits", m_numMcIntHits=5);
    declareProperty("numAcdDigi", m_numAcdDigi=4);
    declareProperty("numCalDigi", m_numCalDigi=12);
	declareProperty("numTkrDigi", m_numTkrDigi=30);
    declareProperty("numCalXtals", m_numCalXtals=12);
    declareProperty("numTkrTracks", m_numTkrTracks=5);
	
}

StatusCode createFakeTdsDataAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();
    
    return sc;
    
}

StatusCode createFakeTdsDataAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    static unsigned int ievent = 0;

    // Retrieve the Event data for this event to set the ids
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) return sc;
    evt->setEvent(ievent);
    evt->setRun(1);
    ievent++;

    sc = storeMcData();

    sc = storeGemData();
    sc = storeDigiData();

    sc = storeReconData();

    return sc;
}

StatusCode createFakeTdsDataAlg::storeMcData() {
    // Purpose and Method: Create fake Monte Carlo data for the TDS.  
    // TDS Output:  EventModel::MC::McParticleCol, 
    //    EventModel::MC::McPositionHitCol, EventModel::MC::McIntegratingHitCol

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	// Make sure MC branch is created on TDS
	SmartDataPtr<Event::MCEvent> mcEvt(eventSvc(), EventModel::MC::Event);
    if(!mcEvt) {
        log << MSG::WARNING << EventModel::MC::Event  <<" could not be registered on data store" << endreq;
        return sc;
    }
	unsigned int run = 4;
	int sourceId = 7;
	unsigned int sequence = 3;
	mcEvt->initialize(run, sourceId, sequence, 99.0);

    // create the TDS location for the McParticle Collection
    Event::McParticleCol* mcParticleTdsCol = new Event::McParticleCol;
    sc = eventSvc()->registerObject(EventModel::MC::McParticleCol, mcParticleTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McParticle Collection" << endreq;
        return sc;
    }

	Event::McParticle *daughterPart = 0;

    unsigned int ipart;
    for (ipart = 0; ipart < m_numMcParticles; ipart++) {
        Event::McParticle *mcPart = new Event::McParticle();
        Event::McParticle *mom = mcPart;
        int id = ipart;
        unsigned int statusBits = 0;

        HepLorentzVector initialMom(ipart, ipart, ipart, ipart);
        double rand = RandGauss::shoot();
        HepLorentzVector finalMom(rand*ipart, rand, ipart*2., ipart);
        HepPoint3D initPos(0.5, 0.25, 0.3);
        HepPoint3D finalPos(1.0, 5.5, 10.3);

        // Setup the TDS version of the McParticle
        mcPart->init(mom, id, statusBits, initialMom, 
            finalMom, initPos, finalPos);

		if (ipart != 0) mcPart->addDaughter(mcPart);

        // Add the TDS McParticle to the TDS collection of McParticles
        mcParticleTdsCol->push_back(mcPart);
		daughterPart = mcPart;
    }


    // create the TDS location for the McParticle Collection
    Event::McPositionHitVector* mcPositionHitTdsCol = new Event::McPositionHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McPositionHitCol, mcPositionHitTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register McPositionHit Collection" << endreq;
        return sc;
    }

    unsigned int iposHit;
    for (iposHit = 0; iposHit < m_numMcPosHits; iposHit++) {
        Event::McPositionHit *posHitTds = new Event::McPositionHit();
        double edep = RandFlat::shoot();
        HepPoint3D entry(RandGauss::shoot(), RandGauss::shoot()*iposHit, iposHit);
        HepPoint3D exit(iposHit*2., iposHit*3., iposHit*4.);
        idents::VolumeIdentifier volId;
        volId.append(1);
        // setup the TDS McPositionHit
        posHitTds->init(edep, volId, entry, exit);
        // add the McPositionHit to the TDS collection of McPositionHits
        mcPositionHitTdsCol->push_back(posHitTds);
    }

    // create the TDS location for the McParticle Collection
    Event::McIntegratingHitVector* mcIntHitCol = new Event::McIntegratingHitVector;
    sc = eventSvc()->registerObject(EventModel::MC::McIntegratingHitCol, mcIntHitCol);
    if (sc.isFailure()) {
        log << "Failed to register McIntegratingHit" << endreq;
        return sc;
    }

    unsigned int intHit;
    for (intHit = 0; intHit < m_numMcIntHits; intHit++) {
        Event::McIntegratingHit *intHitTds = new Event::McIntegratingHit();
        idents::VolumeIdentifier volId;
        volId.append(1);
        intHitTds->setVolumeID(volId);
        mcIntHitCol->push_back(intHitTds);
    }

    return sc;
}


StatusCode createFakeTdsDataAlg::storeGemData() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
   
    // create the TDS location for the Gem
    LdfEvent::Gem *gemTds= new LdfEvent::Gem;
    LdfEvent::GemTileList tileListTds(1,2, 3, 4, 5, 6, 7);
    gemTds->initTrigger(8, 9, 10, 11, 12, 13, tileListTds);
    LdfEvent::GemOnePpsTime ppsTimeTds(14, 15);
    gemTds->initSummary(16, 17, 18, 19, 20, ppsTimeTds, 21);

    sc = eventSvc()->registerObject("/Event/Gem", gemTds);
    if (sc.isFailure()) {
        log << "Failed to register Gem" << endreq;
        return StatusCode::FAILURE;
    }
    return sc;

}
StatusCode createFakeTdsDataAlg::storeDigiData() {
    // Purpose and Method:  Create fake digitization data for ACD, CAL,
    //   TKR.  The data is stored on the TDS.
    // TDS Output: EventModel::Digi::AcdDigiCol, EventModel::Digi::CalDigiCol,
    //    EventModel::Digi::TkrDigiCol

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
   
    // create the TDS location for the AcdDigi Collection
    Event::AcdDigiCol* acdDigiTdsCol = new Event::AcdDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, acdDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register AcdDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    unsigned int iacdDigi;
    for (iacdDigi = 0; iacdDigi < m_numAcdDigi; iacdDigi++) {
        idents::AcdId id(0, 0, 3, 4);
        idents::VolumeIdentifier volId;
        volId.append(2);
        double energy = RandFlat::shoot();
        unsigned short pha[2] = {1, 2096};
        bool veto[2] = {true, true};
        bool low[2] = {true, true};
        bool high[2] = {false, true};
        Event::AcdDigi *acdDigiTds = new Event::AcdDigi(id, volId, 
            energy, pha, veto, low, high);

        acdDigiTdsCol->push_back(acdDigiTds);
    }

    // create the TDS location for the CalDigi Collection
    Event::CalDigiCol* calDigiTdsCol = new Event::CalDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::CalDigiCol, calDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register CalDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    unsigned int icalDigi;
    for (icalDigi = 0; icalDigi < m_numCalDigi; icalDigi++) {
        Event::CalDigi *calDigiTds = new Event::CalDigi();
        idents::CalXtalId idTds(4, 2, 1);
        calDigiTds->initialize(idents::CalXtalId::BESTRANGE, idTds);
        Event::CalDigi::CalXtalReadout r(idents::CalXtalId::LEX8, 1, 
            idents::CalXtalId::HEX8, 4095);
        calDigiTds->addReadout(r);
        calDigiTdsCol->push_back(calDigiTds);
    }

    // create the TDS location for the TkrDigi Collection
    Event::TkrDigiCol* tkrDigiTdsCol = new Event::TkrDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::TkrDigiCol, tkrDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register TkrDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    unsigned int itkrDigi;
    for (itkrDigi = 0; itkrDigi < m_numTkrDigi; itkrDigi++) {
        idents::TowerId towerTds(2, 3);
        int totTds[2] = {RandFlat::shoot(), RandFlat::shoot()};
        Event::TkrDigi *tkrDigiTds = 
            new Event::TkrDigi(5, idents::GlastAxis::X, towerTds, totTds);
        tkrDigiTds->addC0Hit(itkrDigi);
        tkrDigiTds->addC1Hit(itkrDigi+1);
        tkrDigiTdsCol->push_back(tkrDigiTds);
    }

    return sc;
}


StatusCode createFakeTdsDataAlg::storeReconData() {
    // Purpose and Method: Create fake reconstruction data for ACD, CAL, and
    //   TKR.
    // TDS Outputs:

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

StatusCode createFakeTdsDataAlg::finalize()
{    
    return StatusCode::SUCCESS;
}

