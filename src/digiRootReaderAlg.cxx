#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/CalDigi.h"
#include "idents/CalXtalId.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TCollection.h"  // Declares TIter

#include "digiRootData/DigiEvent.h"

#include "facilities/Util.h"


/** @class digiRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header$
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

    /// Reads CAL digi data from ROOT and put data on TDS
    StatusCode readCalDigi();

    /// Closes the ROOT file
    void close();
   
    /// ROOT file pointer
    TFile *m_digiFile;
    /// ROOT tree pointer
    TTree *m_digiTree;
    /// Top-level Monte Carlo ROOT object
    DigiEvent *m_digiEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the Monte Carlo TTree stored in the ROOT file
    std::string m_treeName;

};

static const AlgFactory<digiRootReaderAlg>  Factory;
const IAlgFactory& digiRootReaderAlgFactory = Factory;


digiRootReaderAlg::digiRootReaderAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input pararmeters that may be set via the jobOptions file
    // Input ROOT file name
    declareProperty("digiRootFile",m_fileName="digi.root");
    // Input TTree name
    declareProperty("digiTreeName", m_treeName="Digi");

}

StatusCode digiRootReaderAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for reading.

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    m_digiFile = new TFile(m_fileName.c_str(), "READ");
    if (!m_digiFile->IsOpen()) sc = StatusCode::FAILURE;
    m_digiFile->cd();
    m_digiTree = (TTree*)m_digiFile->Get(m_treeName.c_str());
    m_digiEvt = 0;
    m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvt);
    
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

    static UInt_t evtId = 0;
    m_digiTree->GetEvent(evtId);

    sc = readDigiEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to read top level DigiEvent" << endreq;
        return sc;
    }

    sc = readCalDigi();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to load CalDigi" << endreq;
        return sc;
    }

    m_digiEvt->Clear();
    evtId++;
    
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
    if (runIdTds != runIdTds) evt->setRun(runIdRoot);

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

    return sc;
}

StatusCode digiRootReaderAlg::readCalDigi() {
    MsgStream log(msgSvc(), name());

    StatusCode sc = StatusCode::SUCCESS;
    const TObjArray *calDigiRootCol = m_digiEvt->getCalDigiCol();
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
    while (calDigiRoot = (CalDigi*)calDigiIter.Next()) {
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

void digiRootReaderAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kOverWrite parameter is used in the Write method
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.

    TDirectory *saveDir = gDirectory;
    m_digiFile->cd();
    m_digiFile->Close();
    saveDir->cd();
}

StatusCode digiRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}