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

#include "LdfEvent/DiagnosticData.h"
#include "LdfEvent/EventSummaryData.h"
#include "LdfEvent/LdfTime.h"
#include "LdfEvent/Gem.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter

#include "digiRootData/DigiEvent.h"

#include "facilities/Util.h"

#include "commonData.h"

#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

/** @class digiRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootReaderAlg.cxx,v 1.43.2.3 2005/01/25 19:44:30 heather Exp $
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

    /// Reads in the EM Diagnostic trigger primitive data
    StatusCode readDiagnostic();

    /// Reads ACD digi data from ROOT and puts data on TDS
    StatusCode readAcdDigi();

    /// Reads CAL digi data from ROOT and puts data on TDS
    StatusCode readCalDigi();

    /// Reads TKR digi data from ROOT and puts it on the TDS
    StatusCode readTkrDigi();

    /// Closes the ROOT file
    void close();

    /// Converts from ROOT's VolumeIdentifier to idents::VolumeIdentifier 
    void convertVolumeId(VolumeIdentifier rootVolId, 
        idents::VolumeIdentifier &tdsVolId);
   
    /// ROOT file pointer
    TFile *m_digiFile;
    /// ROOT tree pointer
    TChain *m_digiTree;
    /// Top-level Monte Carlo ROOT object
    DigiEvent *m_digiEvt;
    /// name of the input ROOT file
    std::string m_fileName;
    /// Array of input file names
    StringArrayProperty m_fileList;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;
    /// Stores number of events available in the input ROOT TTree
    Long64_t m_numEvents;
  
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};

static const AlgFactory<digiRootReaderAlg>  Factory;
const IAlgFactory& digiRootReaderAlgFactory = Factory;


digiRootReaderAlg::digiRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
    declareProperty("digiRootFile",m_fileName="");
    StringArrayProperty initList;
    std::vector<std::string> initVec;
    initVec.push_back("digi.root");
    initList.setValue(initVec);
    declareProperty("digiRootFileList",m_fileList=initList);
    // Input TTree name
    initVec.clear();
    declareProperty("digiTreeName", m_treeName="Digi");

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
        //return StatusCode::FAILURE;
    } 

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   

    m_digiTree = new TChain(m_treeName.c_str());

    std::string emptyStr("");
    if (m_fileName.compare(emptyStr) != 0) {
	  TFile f(m_fileName.c_str());
      if (!f.IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName.c_str()
            << " could not be opened for reading." << endreq;
        return StatusCode::FAILURE;
      }
	  f.Close();
	  m_digiTree->Add(m_fileName.c_str());
          log << MSG::INFO << "Opened file: " << m_fileName.c_str() << endreq;
    } else {
      const std::vector<std::string> fileList = m_fileList.value( );
      std::vector<std::string>::const_iterator it;
      std::vector<std::string>::const_iterator itend = fileList.end( );
      for (it = fileList.begin(); it != itend; it++) {
        std::string theFile = (*it);
	    TFile f(theFile.c_str());
        if (!f.IsOpen()) {
          log << MSG::ERROR << "ROOT file " << theFile.c_str()
              << " could not be opened for reading." << endreq;
          return StatusCode::FAILURE;
        }
	  f.Close();
	  m_digiTree->Add(theFile.c_str());
          log << MSG::INFO << "Opened file: " << theFile.c_str() << endreq;
	  }
    }


    m_digiEvt = 0;
    m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvt);
    m_common.m_digiEvt = m_digiEvt;

    m_numEvents = m_digiTree->GetEntries();
    
    if (m_rootIoSvc) {
      m_rootIoSvc->setRootEvtMax(m_numEvents);
      if (!m_digiTree->GetIndex()) m_digiTree->BuildIndex("m_runId", "m_eventId");
      m_rootIoSvc->registerRootTree(m_digiTree);
    }


    saveDir->cd();
    return sc;
    
}

StatusCode digiRootReaderAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the ROOT file and store
    //   data on the TDS.

    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    
    if (m_digiEvt) m_digiEvt->Clear();

    static Long64_t evtId = 0;
    Long64_t readInd;
    int numBytes;
    std::pair<int,int> runEventPair = (m_rootIoSvc) ? m_rootIoSvc->runEventPair() : std::pair<int,int>(-1,-1);
	
    if ((m_rootIoSvc) && (m_rootIoSvc->useIndex())) {
        readInd = m_rootIoSvc->index();
    } else if ((m_rootIoSvc) && (m_rootIoSvc->useRunEventPair())) {
        int run = runEventPair.first;
        int evt = runEventPair.second;
        readInd = m_digiTree->GetEntryNumberWithIndex(run, evt);
    } else {
        readInd = evtId;
    }

    if (readInd >= m_numEvents) {
        log << MSG::WARNING << "Requested index is out of bounds - no digi data loaded" << endreq;
        return StatusCode::SUCCESS;
    }

    if (m_rootIoSvc) m_rootIoSvc->setActualIndex(readInd);

    // ADDED FOR THE FILE HEADERS DEMO
    m_digiTree->LoadTree(readInd);
    m_headersTool->readConstDigiHeader(m_digiTree->GetFile()) ;
    
    numBytes = m_digiTree->GetEvent(readInd);
	
    if ((numBytes <= 0) || (!m_digiEvt)) {
        log << MSG::WARNING << "Failed to load digi event" << endreq;
        return StatusCode::SUCCESS;
    }


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

	evtId = readInd+1;
    
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

    // Check to see if the event and run ids have already been set.
    if (eventIdTds != eventIdRoot) evt->setEvent(eventIdRoot);
    if (runIdTds != runIdRoot) evt->setRun(runIdRoot);

    TimeStamp timeObj(m_digiEvt->getTimeStamp());
    evt->setTime(timeObj);

    evt->setLivetime(m_digiEvt->getLiveTime());

    evt->setTrigger(m_digiEvt->getL1T().getTriggerWord());

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
    sc = eventSvc()->registerObject("/Event/Digi/TriRowBits", rowbits);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Could not register TriRowBits" << endreq;
        return sc;
    }
    unsigned int iTower = 0;
    for (iTower = 0; iTower < 16; iTower++) {
        rowbits->setTriRowBits(iTower, m_digiEvt->getL1T().getTriRowBits(iTower));
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
    unsigned summaryWord = m_digiEvt->getEventSummaryData().summary();
    unsigned eventFlags = m_digiEvt->getEventSummaryData().eventFlags();
    EventSummaryData evtSummary = m_digiEvt->getEventSummaryData();

    LdfEvent::EventSummaryData *evtSumTds = new LdfEvent::EventSummaryData();
    evtSumTds->initialize(summaryWord);
    evtSumTds->initEventFlags(eventFlags);
    const unsigned int nTem = 16;
    unsigned int iTem;
    unsigned int tem[nTem];
    for (iTem = 0; iTem < nTem; iTem++) 
        tem[iTem] = evtSummary.temLength(iTem);
    evtSumTds->initContribLen(tem, evtSummary.gemLength(), evtSummary.oswLength(),
        evtSummary.errLength(), evtSummary.diagLength(), evtSummary.aemLength());
    sc = eventSvc()->registerObject("/Event/EventSummary", evtSumTds);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "could not register /Event/EventSummary " << endreq;
        return sc;
    }
    return sc;
}

StatusCode digiRootReaderAlg::readGem() {

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    const Gem &gemRoot = m_digiEvt->getGem();
    GemTileList tileListRoot = gemRoot.getTileList();
    LdfEvent::Gem *gemTds = new LdfEvent::Gem();
    LdfEvent::GemTileList tileListTds(tileListRoot.getXzm(), tileListRoot.getXzp(), 
              tileListRoot.getYzm(), tileListRoot.getYzp(), tileListRoot.getXy(), 
              tileListRoot.getRbn(), tileListRoot.getNa());
    gemTds->initTrigger(gemRoot.getTkrVector(), gemRoot.getRoiVector(),
            gemRoot.getCalLeVector(), gemRoot.getCalHeVector(),
            gemRoot.getCnoVector(), gemRoot.getConditionSummary(),
            tileListTds);
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
        log << MSG::ERROR << "could not register /Event/Gem " << endreq
;
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
        Event::AcdDigi::ParityError err[2];
        err[0] = (acdDigiRoot->getParityError(AcdDigi::A) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        err[1] = (acdDigiRoot->getParityError(AcdDigi::B) == AcdDigi::NOERROR) ? Event::AcdDigi::NOERROR : Event::AcdDigi::ERROR;
        acdDigiTds->initLdfParameters(acdDigiRoot->getTileName(), acdDigiRoot->getTileNumber(),
                    range, err);
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
    }

    return sc;
}

void digiRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    //TDirectory *saveDir = gDirectory;
    //m_digiFile->cd();
    //m_digiFile->Close();
    //saveDir->cd();
    if (m_digiTree) {
        delete m_digiTree ;
        m_digiTree = 0 ;
    }
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
