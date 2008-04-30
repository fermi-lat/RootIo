#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"
#include "idents/CalXtalId.h"
#include "idents/TowerId.h"

#include "Trigger/TriRowBits.h"

#include "OnboardFilterTds/FilterStatus.h"

#include "LdfEvent/DiagnosticData.h"
#include "LdfEvent/EventSummaryData.h"
#include "LdfEvent/LdfTime.h"
#include "LdfEvent/Gem.h"
#include "LdfEvent/ErrorData.h"

#include "AncillaryDataEvent/Digi.h"
#include "LdfEvent/LsfMetaEvent.h"
#include "LdfEvent/LsfCcsds.h"

#include "digiRootData/DigiEvent.h"

#include "facilities/Util.h"

#include "commonData.h"

#include "RootConvert/Digi/LsfDigiConvert.h"
#include "RootConvert/Digi/OnboardFilterConvert.h"
#include "RootConvert/Digi/AdfDigiConvert.h"

#include "RootIo/IRootIoSvc.h"
#include "RootConvert/Utilities/RootReaderUtil.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

/** @class digiRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootReaderAlg.cxx,v 1.94 2008/03/24 15:25:43 heather Exp $
 */

class digiRootReaderAlg : public Algorithm
{	
public:
    
    digiRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in read mode and creating a new TTree
    StatusCode initialize();
   
    /// Orchastrates reading from ROOT file and storing the data on the TDS for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
            
private:

    /// Reads top-level DigiEvent
    StatusCode readDigiEvent();

    /// Reads in EM event summary word
    StatusCode readEventSummary();

    StatusCode readGem();
 
    /// Reads in LDF Error data
    StatusCode readError();

    /// Reads in the EM Diagnostic trigger primitive data
    StatusCode readDiagnostic();

    /// Reads ACD digi data from ROOT and puts data on TDS
    StatusCode readAcdDigi();

    /// Reads CAL digi data from ROOT and puts data on TDS
    StatusCode readCalDigi();

    /// Reads TKR digi data from ROOT and puts it on the TDS
    StatusCode readTkrDigi();

    /// Reads OBF FilterStatus from ROOT and puts it on the TDS
    StatusCode readFilterStatus();

    /// Reads the Meta Event from ROOT and puts it on the TDS
    StatusCode readMetaEvent();

    /// Reads CCSDS froM ROOT and puts it on the TDS
    StatusCode readCcsds();

    /// Reads Ancillary Beamtest Digi data and puts it on the TDS
    StatusCode readAdf();

    /// Closes the ROOT file
    void close();

    /// Converts from ROOT's VolumeIdentifier to idents::VolumeIdentifier 
    void convertVolumeId(VolumeIdentifier rootVolId, 
        idents::VolumeIdentifier &tdsVolId);
   
    /// Top-level Monte Carlo ROOT object
    DigiEvent *m_digiEvt;
//    /// name of the input ROOT file
    std::string m_fileName;
    /// Array of input file names
    StringArrayProperty m_fileList;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Stores number of events available in the input ROOT TTree
    /// Branch name for events
    std::string m_branchName;
    /// Option string which will be passed to McEvent::Clear
    std::string m_clearOption;

    Long64_t m_numEvents;
  
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};

static const AlgFactory<digiRootReaderAlg>  Factory;
const IAlgFactory& digiRootReaderAlgFactory = Factory;


digiRootReaderAlg::digiRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator), m_digiEvt(0)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name  provided for backward compatibility, digiRootFileList is preferred
    // The files will be ignored if RootIoSvc is provided a meta (event collection) ROOT file
    declareProperty("digiRootFile",m_fileName="");
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initList.setValue(initVec);
    declareProperty("digiRootFileList",m_fileList=initList);
    // Input TTree name
    initVec.clear();
    declareProperty("digiTreeName", m_treeName="Digi");
    declareProperty("digiBranchName", m_branchName="DigiEvent");
    declareProperty("clearOption", m_clearOption="");
}

StatusCode digiRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    StatusCode headersSc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
   
    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "Event loop will not terminate gracefully" << endreq;
        m_rootIoSvc = 0;
        // We truly depend on the RootIoSvc now to read/write files, so we cannot just
        // continue if this service is not found.
        return StatusCode::FAILURE;
    } 

    if ( (m_fileList.value().size() > 0) && ( !m_fileName.empty() )) {
        log << MSG::WARNING << "Both digiRootFile and digiRootFileList have "
            << "been specified, digiRootFile is deprecated, please use "
            << "digiRootFileList" << endreq;
         return StatusCode::FAILURE;
    } else if ( (m_fileList.value().size() == 0) && ( !m_fileName.empty() ) )
        m_rootIoSvc->appendFileList(m_fileList, m_fileName);
    else if (m_fileList.value().size() == 0)
        m_rootIoSvc->appendFileList(m_fileList, "digi.root");



    // Set up new school system...
    // Use the TTree name as the key type 
    m_rootIoSvc->prepareRootInput("digi", m_treeName, m_branchName, m_fileList);

    return sc;
    
}

StatusCode digiRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    if (m_digiEvt) m_digiEvt->Clear(m_clearOption.c_str());
    m_digiEvt = 0;

    // Try reading the event this way... 
    // using treename as the key
    m_digiEvt = dynamic_cast<DigiEvent*>(m_rootIoSvc->getNextEvent("digi"));

    if (!m_digiEvt) return StatusCode::FAILURE;

    // Clear the digi common maps
    m_common.m_rootTkrDigiMap.clear();
    m_common.m_rootCalDigiMap.clear();

    // Clear the digi common maps
    m_common.m_tkrDigiMap.clear();
    m_common.m_calDigiMap.clear();

    sc = readDigiEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read top level DigiEvent" << endreq;
        return sc;
    }

    sc = readEventSummary();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read in eventsummary data" << endreq;
        return sc;
    }

    sc = readGem();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read in GEM" << endreq;
        return sc;
    }
    
    sc = readError();
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to read in error data" << endreq;
        // Do not return failure - error data may not be in all data
    }

    sc = readDiagnostic();
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to read in diagnostic data" << endreq;
        // Do not return failure - diagnostic data may not be in all data
    }


    sc = readAcdDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load AcdDigi" << endreq;
        return sc;
    }

    sc = readCalDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load CalDigi" << endreq;
        return sc;
    }

    sc = readTkrDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load TkrDigi" << endreq;
        return sc;
    }

    sc = readFilterStatus();
    // do not terminate job due to missing FilterStatus in ROOT file
    // older runs will not have this branch filled.
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to load FilterStatus" << endreq;
        // reset status code to success so we don't propagate the error
        sc = StatusCode::SUCCESS;
    }


    sc = readMetaEvent();
    // do not terminate job due to missing MetaEvent in ROOT file
    // MC events and older runs will not have this branch filled.
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to load MetaEvent" << endreq;
        // reset status code to success so we don't propagate the error
        sc = StatusCode::SUCCESS;
    }

    sc = readCcsds();
    // do not terminate job due to missing CCSDS in ROOT file
    // MC events and older runs will not have this branch filled.
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to load CCSDS" << endreq;
        // reset status code to success so we don't propagate the error
        sc = StatusCode::SUCCESS;
    }

    sc = readAdf();
    // do not terminate job due to missing ancillary data in ROOT file
    // MC events and older runs will not have this branch filled.
    if (sc.isFailure()) {
        log << MSG::DEBUG << "Failed to load ancillary data" << endreq;
        // reset status code to success so we don't propagate the error
        sc = StatusCode::SUCCESS;
    }
    else
        log << MSG::VERBOSE << "loaded ancillary data" << endreq;

    
    return sc;
}


StatusCode digiRootReaderAlg::readDigiEvent() {

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evt(eventSvc(), EventModel::EventHeader);
    if (!evt) {
        log << MSG::ERROR << "Failed to retrieve Event" << endreq;
        return StatusCode::FAILURE;
    }

    unsigned int eventIdTds = evt->event();
    unsigned int runIdTds = evt->run();

    unsigned int eventIdRoot = m_digiEvt->getEventId();
    unsigned int runIdRoot = m_digiEvt->getRunId();

    log << MSG::DEBUG << "Reading Event (run, event): (" << runIdRoot
        << ", " << eventIdRoot << ")" << endreq;

    // Check to see if the event and run ids have already been set.
    if (eventIdTds != eventIdRoot) evt->setEvent(eventIdRoot);
    if (runIdTds != runIdRoot) evt->setRun(runIdRoot);

    log << MSG::DEBUG << "Reading Event (run, event): (" << runIdRoot
        << ", " << eventIdRoot << ")" << endreq;

    TimeStamp timeObj(m_digiEvt->getTimeStamp());
    evt->setTime(timeObj);

    evt->setLivetime(m_digiEvt->getLiveTime());

    evt->setTrigger(m_digiEvt->getL1T().getTriggerWord());
    evt->setTriggerWordTwo(m_digiEvt->getL1T().getTriggerWordTwo());
    evt->setGemPrescale(m_digiEvt->getL1T().getGemPrescale());
    evt->setGltPrescale(m_digiEvt->getL1T().getGltPrescale());
    evt->setPrescaleExpired(m_digiEvt->getL1T().getPrescaleExpired());

    Event::DigiEvent* digiEventTds = 
        SmartDataPtr<Event::DigiEvent>(eventSvc(), EventModel::Digi::Event);
    if (!digiEventTds) {
        sc = eventSvc()->registerObject(EventModel::Digi::Event /*"/Event/Digi"*/,new DataObject);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register " << EventModel::Digi::Event /*<< /Event/Digi "*/ << endreq;
            return sc;
        }
    } else {
        bool fromMc = m_digiEvt->getFromMc();
        digiEventTds->initialize(fromMc);
    }

    TriRowBitsTds::TriRowBits *rowbits= new TriRowBitsTds::TriRowBits;
    sc = eventSvc()->registerObject("/Event/TriRowBits", rowbits);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Could not register TriRowBits" << endreq;
        return sc;
    }
    unsigned int iTower = 0;
    for (iTower = 0; iTower < 16; iTower++) {
        rowbits->setDigiTriRowBits(iTower, m_digiEvt->getL1T().getDigiTriRowBits(iTower));
        rowbits->setTrgReqTriRowBits(iTower, m_digiEvt->getL1T().getTrgReqTriRowBits(iTower));
    }
    
    LdfEvent::LdfTime *ldfTimeTds = new LdfEvent::LdfTime();
    if (ldfTimeTds) {
        ldfTimeTds->initialize(m_digiEvt->getEbfTimeSec(), m_digiEvt->getEbfTimeNanoSec(), m_digiEvt->getEbfUpperPpcTimeBase(), m_digiEvt->getEbfLowerPpcTimeBase());
        sc = eventSvc()->registerObject("/Event/Time", ldfTimeTds);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register /Event/Time " << endreq;
            return sc;
        }
    }
    return sc;
}

StatusCode digiRootReaderAlg::readEventSummary() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    EventSummaryData &evtSummary = m_digiEvt->getEventSummaryData();
    unsigned summaryWord = evtSummary.summary();
    unsigned eventFlags = evtSummary.eventFlags();
    //unsigned int evtSeq = 0; //evtSummary.eventSequence();
    unsigned long evtSizeInBytes = evtSummary.eventSizeInBytes();

    // Only update the eventflags on the TDS if the /Event/EventSummary
    // does not yet exist (digiRootReader may fill this for us)
    SmartDataPtr<LdfEvent::EventSummaryData> evtSumTds(eventSvc(), "/Event/EventSummary");
    if (!evtSumTds) {
      evtSumTds = new LdfEvent::EventSummaryData();
      sc = eventSvc()->registerObject("/Event/EventSummary", evtSumTds);
      if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register /Event/EventSummary " << endreq
;
        return sc;
      }
    }

    evtSumTds->initialize(summaryWord);
    evtSumTds->initEventFlags(eventFlags);
    evtSumTds->initEventSizeInBytes(evtSizeInBytes);

    evtSumTds->initContribLen((unsigned int*)evtSummary.temLength(), evtSummary.gemLength(), 
                              evtSummary.oswLength(), (unsigned int*)evtSummary.errLength(), 
                              (unsigned int*)evtSummary.diagLength(), evtSummary.aemLength());
    return sc;
}

StatusCode digiRootReaderAlg::readGem() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    // Skip GEM if this is simulated data
    if (m_digiEvt->getFromMc()) return sc;

    const Gem &gemRoot = m_digiEvt->getGem();
    GemTileList tileListRoot = gemRoot.getTileList();
    LdfEvent::Gem *gemTds = new LdfEvent::Gem();
    LdfEvent::GemTileList tileListTds(tileListRoot.getXzm(), tileListRoot.getXzp(), 
              tileListRoot.getYzm(), tileListRoot.getYzp(), tileListRoot.getXy(), 
              tileListRoot.getRbn(), tileListRoot.getNa());

    gemTds->initTrigger(gemRoot.getTkrVector(), gemRoot.getRoiVector(),
            gemRoot.getCalLeVector(), gemRoot.getCalHeVector(),
            gemRoot.getCnoVector(), gemRoot.getConditionSummary(),
            gemRoot.getMissed(), tileListTds);

    LdfEvent::GemOnePpsTime ppsTimeTds(gemRoot.getOnePpsTime().getTimebase(),
                            gemRoot.getOnePpsTime().getSeconds());
    LdfEvent::GemDataCondArrivalTime gemCondTimeTds;
    gemCondTimeTds.init(gemRoot.getCondArrTime().condArr());
    gemTds->initSummary(gemRoot.getLiveTime(), gemRoot.getPrescaled(),
                        gemRoot.getDiscarded(), gemCondTimeTds,
                        gemRoot.getTriggerTime(), ppsTimeTds, 
                        gemRoot.getDeltaEventTime(), 
                        gemRoot.getDeltaWindowOpenTime());

    sc = eventSvc()->registerObject("/Event/Gem", gemTds);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register /Event/Gem " << endreq;
        return sc;
    }
    return sc;
}

StatusCode digiRootReaderAlg::readError() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    LdfEvent::ErrorData *errCol = new LdfEvent::ErrorData();

    const TClonesArray *temRootCol = m_digiEvt->getTemCol();
    TIter temRootIt(temRootCol);
    Tem* temCur = 0;
    while ((temCur = (Tem*)temRootIt.Next())) {

        const ErrorData& errRoot = temCur->getError();
        LdfEvent::TowerErrorData err(temCur->getTowerId(), errRoot.getCal(), 
                   errRoot.getTkr(), errRoot.getPhs(), errRoot.getTmo());
        const UChar_t *tkrFifoCol = errRoot.getTkrFifoFullCol();
        unsigned int igtcc;
        for (igtcc=0;igtcc<enums::numGtcc;igtcc++)
            err.setTkrFifoFull(igtcc,tkrFifoCol[igtcc]);
        errCol->addTowerError(err);
    }

    sc = eventSvc()->registerObject("/Event/Error", errCol);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register " << "/Event/Error" << endreq;
        return sc;
    }
    return sc;
}

StatusCode digiRootReaderAlg::readDiagnostic() {

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;

    const TClonesArray *calCol = m_digiEvt->getCalDiagnosticCol();
    const TClonesArray *tkrCol = m_digiEvt->getTkrDiagnosticCol();
    LdfEvent::DiagnosticData *diagTds = new LdfEvent::DiagnosticData();
    TIter calIt(calCol);
    CalDiagnosticData *cDiagRoot;
    while ((cDiagRoot = (CalDiagnosticData*)calIt.Next())!=0) {
        LdfEvent::CalDiagnosticData cDiagTds(cDiagRoot->getDataWord(),cDiagRoot->tower(),cDiagRoot->layer());
        diagTds->addCalDiagnostic(cDiagTds);
    }

    TIter tkrIt(tkrCol);
    TkrDiagnosticData *tDiagRoot;
    while((tDiagRoot = (TkrDiagnosticData*)tkrIt.Next())!=0) {
        LdfEvent::TkrDiagnosticData tDiagTds(tDiagRoot->getDataWord(),tDiagRoot->tower(),tDiagRoot->gtcc());
        diagTds->addTkrDiagnostic(tDiagTds);
    }

    sc = eventSvc()->registerObject("/Event/Diagnostic", diagTds);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register " << "/Event/Diagnostic" << endreq;
        return sc;
    }
    

    return sc;
}


StatusCode digiRootReaderAlg::readAcdDigi() {
    // Purpose and Method:  Read in ACD digi collection from ROOT file
    //  and insert values into ACD digi collection on the TDS.

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const TClonesArray *acdDigiRootCol = m_digiEvt->getAcdDigiCol();
    if (!acdDigiRootCol) return sc;
    TIter acdDigiIter(acdDigiRootCol);

    // create the TDS location for the AcdDigi Collection
    Event::AcdDigiCol* acdDigiTdsCol = new Event::AcdDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, acdDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register AcdDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    AcdDigi *acdDigiRoot = 0;
    while ((acdDigiRoot = (AcdDigi*)acdDigiIter.Next())!=0) {
        float energyTds = acdDigiRoot->getEnergy();
        
        AcdId idRoot = acdDigiRoot->getId();
        idents::AcdId idTds(idRoot.getId(2));

        // Rather than using the volId stored in the ROOT file - recompute using the 
        // routine in the idents package.
        idents::VolumeIdentifier volIdTds = idTds.volId();

        unsigned short phaTds[2] = { acdDigiRoot->getPulseHeight(AcdDigi::A),
            acdDigiRoot->getPulseHeight(AcdDigi::B) };
        bool vetoTds[2] = { acdDigiRoot->getVeto(AcdDigi::A),
            acdDigiRoot->getVeto(AcdDigi::B)};
        bool lowTds[2] = { acdDigiRoot->getLowDiscrim(AcdDigi::A),
            acdDigiRoot->getLowDiscrim(AcdDigi::B) };
        bool highTds[2] = { acdDigiRoot->getHighDiscrim(AcdDigi::A),
            acdDigiRoot->getHighDiscrim(AcdDigi::B) };

        Event::AcdDigi *acdDigiTds = new Event::AcdDigi(idTds, volIdTds, 
            energyTds,phaTds, vetoTds, lowTds, highTds);

        Event::AcdDigi::Range range[2];
        range[0] = (acdDigiRoot->getRange(AcdDigi::A) == AcdDigi::LOW) ? Event::AcdDigi::LOW : Event::AcdDigi::HIGH;
        range[1] = (acdDigiRoot->getRange(AcdDigi::B) == AcdDigi::LOW) ? Event::AcdDigi::LOW : Event::AcdDigi::HIGH;

        Event::AcdDigi::ParityError err[4];
        err[0] = (acdDigiRoot->getOddParityError(AcdDigi::A) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        err[1] = (acdDigiRoot->getOddParityError(AcdDigi::B) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        err[2] = (acdDigiRoot->getHeaderParityError(AcdDigi::A) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        err[3] = (acdDigiRoot->getHeaderParityError(AcdDigi::B) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;

        acdDigiTds->initLdfParameters(acdDigiRoot->getTileName(), 
                             acdDigiRoot->getTileNumber(), range, err);

        acdDigiTds->initGem(acdDigiRoot->isNinja(), acdDigiRoot->getGemFlag());

        acdDigiTdsCol->push_back(acdDigiTds);
    }

    return sc;
}

StatusCode digiRootReaderAlg::readCalDigi() {
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    const TClonesArray *calDigiRootCol = m_digiEvt->getCalDigiCol();
    if (!calDigiRootCol) return sc;
    TIter calDigiIter(calDigiRootCol);

    // create the TDS location for the CalDigi Collection
    Event::CalDigiCol* calDigiTdsCol = new Event::CalDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::CalDigiCol, calDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register CalDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    CalDigi *calDigiRoot = 0;
    while ((calDigiRoot = (CalDigi*)calDigiIter.Next())!=0) {
        Event::CalDigi *calDigiTds = new Event::CalDigi();
        CalXtalId::CalTrigMode modeRoot = calDigiRoot->getMode();
        CalXtalId idRoot = calDigiRoot->getPackedId();
        idents::CalXtalId idTds(idRoot.getTower(), idRoot.getLayer(), idRoot.getColumn());
        idents::CalXtalId::CalTrigMode modeTds;
        if (modeRoot == CalXtalId::BESTRANGE) {
            modeTds = idents::CalXtalId::BESTRANGE;
            const CalXtalReadout *readoutRoot = calDigiRoot->getXtalReadout(0);
            Char_t rangePlusRoot = readoutRoot->getRange(CalXtalId::POS);
            UInt_t adcPlusRoot = readoutRoot->getAdc(CalXtalId::POS);
            Char_t rangeMinRoot = readoutRoot->getRange(CalXtalId::NEG);
            UInt_t adcMinRoot = readoutRoot->getAdc(CalXtalId::NEG);
            Event::CalDigi::CalXtalReadout r(rangePlusRoot, adcPlusRoot, 
                rangeMinRoot, adcMinRoot);
            calDigiTds->initialize(modeTds, idTds);
            calDigiTds->addReadout(r);
        } else {
            modeTds = idents::CalXtalId::ALLRANGE;
            calDigiTds->initialize(modeTds, idTds);
            int range;
            for (range = CalXtalId::LEX8; range <= CalXtalId::HEX1; range++) {
                const CalXtalReadout *readoutRoot = calDigiRoot->getXtalReadout(range);
                if (!readoutRoot) continue;
                Char_t rangePlusRoot = readoutRoot->getRange(CalXtalId::POS);
                UInt_t adcPlusRoot = readoutRoot->getAdc(CalXtalId::POS);
                Char_t rangeMinRoot = readoutRoot->getRange(CalXtalId::NEG);
                UInt_t adcMinRoot = readoutRoot->getAdc(CalXtalId::NEG);
                Event::CalDigi::CalXtalReadout r(rangePlusRoot, adcPlusRoot, 
                    rangeMinRoot, adcMinRoot);
                calDigiTds->addReadout(r);
            }
        }

        calDigiTdsCol->push_back(calDigiTds);

        m_common.m_rootCalDigiMap[calDigiRoot] = calDigiTds;
        m_common.m_calDigiMap[calDigiTds]      = calDigiRoot;
    }
 
    return sc;
}

StatusCode digiRootReaderAlg::readTkrDigi() {
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    const TObjArray *tkrDigiRootCol = m_digiEvt->getTkrDigiCol();
    if (!tkrDigiRootCol) return sc;
    TIter tkrDigiIter(tkrDigiRootCol);

    // create the TDS location for the CalDigi Collection
    Event::TkrDigiCol* tkrDigiTdsCol = new Event::TkrDigiCol;
    sc = eventSvc()->registerObject(EventModel::Digi::TkrDigiCol, tkrDigiTdsCol);
    if (sc.isFailure()) {
        log << "Failed to register TkrDigi Collection" << endreq;
        return StatusCode::FAILURE;
    }

    TkrDigi *tkrDigiRoot = 0;
    while ((tkrDigiRoot = (TkrDigi*)tkrDigiIter.Next())!=0) {
        TowerId towerRoot = tkrDigiRoot->getTower();
        idents::TowerId towerTds(towerRoot.ix(), towerRoot.iy());
        GlastAxis::axis axisRoot = tkrDigiRoot->getView();
        idents::GlastAxis::axis axisTds = 
            (axisRoot == GlastAxis::X) ? idents::GlastAxis::X : idents::GlastAxis::Y;
        int totTds[2] = { tkrDigiRoot->getToT(0), tkrDigiRoot->getToT(1)};
        Event::TkrDigi *tkrDigiTds = 
            new Event::TkrDigi(tkrDigiRoot->getBilayer(),
            axisTds, towerTds, totTds);
        int lastController0Strip = tkrDigiRoot->getLastController0Strip();
        unsigned int numStrips = tkrDigiRoot->getNumHits();
        unsigned int iHit;
        for (iHit = 0; iHit < numStrips; iHit++) {
            int strip = tkrDigiRoot->getHit(iHit);
            if (strip <= lastController0Strip) {
                tkrDigiTds->addC0Hit(strip);
            } else {
                tkrDigiTds->addC1Hit(strip);
            }
        }

        tkrDigiTdsCol->push_back(tkrDigiTds);

        m_common.m_rootTkrDigiMap[tkrDigiRoot] = tkrDigiTds;
        m_common.m_tkrDigiMap[tkrDigiTds]      = tkrDigiRoot;
    }

    return sc;
}


StatusCode digiRootReaderAlg::readFilterStatus() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const FilterStatus& filterStatusRoot = m_digiEvt->getFilterStatus();

    // create TDS location
    DataObject *pObj = new DataObject();
    sc = eventSvc()->registerObject("/Event/Filter", pObj);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register /Event/Filter" << endreq;
        return StatusCode::FAILURE;
    }

    OnboardFilterTds::FilterStatus *obfTds = new OnboardFilterTds::FilterStatus;

    sc = eventSvc()->registerObject("/Event/Filter/FilterStatus", obfTds);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register FilterStatus" << endreq;
        return StatusCode::FAILURE;
    }

    RootPersistence::convert(filterStatusRoot, *obfTds);

    // Now do the new style obf filter status class
    const ObfFilterStatus& obfFilterStatusRoot = m_digiEvt->getObfFilterStatus();

    OnboardFilterTds::ObfFilterStatus *obfFilterStatusTds = new OnboardFilterTds::ObfFilterStatus;

    RootPersistence::convert(obfFilterStatusRoot, *obfFilterStatusTds);

    sc = eventSvc()->registerObject("/Event/Filter/ObfFilterStatus", obfFilterStatusTds);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register ObfFilterStatus" << endreq;
        return StatusCode::FAILURE;
    }

    RootPersistence::convert(obfFilterStatusRoot, *obfFilterStatusTds);

    return sc;

}

StatusCode digiRootReaderAlg::readMetaEvent() {
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    const MetaEvent& metaEventRoot = m_digiEvt->getMetaEvent();

    // create the TDS location for the MetaEvent Collection
    LsfEvent::MetaEvent* metaEventTds = new LsfEvent::MetaEvent;
    sc = eventSvc()->registerObject("/Event/MetaEvent", metaEventTds);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register MetaEvent" << endreq;
        return StatusCode::FAILURE;
    }

    RootPersistence::convert(metaEventRoot,*metaEventTds);
    return sc;
}

StatusCode digiRootReaderAlg::readCcsds() {
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    const Ccsds& ccsdsRoot = m_digiEvt->getCcsds();

    // create the TDS location for the CCSDS data 
    LsfEvent::LsfCcsds* ccsdsTds = new LsfEvent::LsfCcsds;
    sc = eventSvc()->registerObject("/Event/Ccsds", ccsdsTds);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register CCSDS" << endreq;
        return StatusCode::FAILURE;
    }

    RootPersistence::convert(ccsdsRoot,*ccsdsTds);
    return sc;
}

StatusCode digiRootReaderAlg::readAdf() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    const AdfDigi& adfRoot = m_digiEvt->getAdfDigi();
    // create TDS location
    DataObject *pObj = new DataObject();
    sc = eventSvc()->registerObject("/Event/AncillaryEvent", pObj);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register /Event/AncillaryEvent" << endreq;
        return StatusCode::FAILURE;
    }

    AncillaryData::Digi* adfTds = new AncillaryData::Digi;
    sc = eventSvc()->registerObject("/Event/AncillaryEvent/Digi", adfTds);
    if (sc.isFailure()) {
        log << MSG::INFO << "Failed to register Adf" << endreq;
        return StatusCode::FAILURE;
    }
    RootPersistence::convert(adfRoot, *adfTds);
    return sc;
}

void digiRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.



}

StatusCode digiRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();
    return sc;
}

void digiRootReaderAlg::convertVolumeId(VolumeIdentifier rootVolId, 
                                      idents::VolumeIdentifier& tdsVolId) 
{
    // Purpose and Method:  We must store the volume ids as two 32 bit UInt_t
    //     in the ROOT class.  The idents::VolumeIdentifier class stores the
    //     data in one 64 bit word.  We must convert from the two 32 bit words
    //     into the 64 bit word.  We perform the conversion by iterating over
    //     all of the ids in the ROOT VolumeIdentifier and appending them to
    //     the TDS idents::VolumeIdentifier.
    // Input:  ROOT VolumeIdentifier
    // Ouput:  idents::VolumeIdentifier

    int index;
    for (index = 0; index < rootVolId.size(); index++) {
        tdsVolId.append(rootVolId.operator [](index));
    }

}
