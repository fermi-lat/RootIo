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

/** @class digiRootReaderAlg
 * @brief Reads Digitization data from a persistent ROOT file and stores the
 * the data in the TDS.
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootReaderAlg.cxx,v 1.16 2003/08/25 18:47:38 heather Exp $
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
    int m_numEvents;
  
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

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
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
   
    if ( service("RootIoSvc", m_rootIoSvc).isFailure() ){
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
	  }
    }


    m_digiEvt = 0;
    m_digiTree->SetBranchAddress("DigiEvent", &m_digiEvt);
    m_common.m_digiEvt = m_digiEvt;

    m_numEvents = m_digiTree->GetEntries();
    
	if (m_rootIoSvc) {
		m_rootIoSvc->setRootEvtMax(m_numEvents);
		//m_digiTree->BuildIndex("m_runId", "m_eventId");
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

    static Int_t evtId = 0;
	int readInd, numBytes;
	std::pair<int,int> runEventPair = (m_rootIoSvc) ? m_rootIoSvc->runEventPair() : std::pair<int,int>(-1,-1);
	
	if ((m_rootIoSvc) && (m_rootIoSvc->index() >= 0)) {
		readInd = m_rootIoSvc->index();
	} else if ((m_rootIoSvc) && (runEventPair.first != -1) && (runEventPair.second != -1)) {
		int run = runEventPair.first;
		int evt = runEventPair.second;
		readInd = m_digiTree->GetEntryWithIndex(run, evt);
	} else {
		readInd = evtId;
	}

    if (readInd >= m_numEvents) {
        log << MSG::WARNING << "Requested index is out of bounds" << endreq;
        return StatusCode::FAILURE;
    }

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
    while (acdDigiRoot = (AcdDigi*)acdDigiIter.Next()) {
        float energyTds = acdDigiRoot->getEnergy();
        
        AcdId idRoot = acdDigiRoot->getId();
        idents::AcdId idTds(idRoot.getLayer(), idRoot.getFace(), 
            idRoot.getRow(), idRoot.getColumn());

        //const VolumeIdentifier volIdRoot = acdDigiRoot->getVolId();
        //convertVolumeId(volIdRoot, volIdTds);
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
    while (tkrDigiRoot = (TkrDigi*)tkrDigiIter.Next()) {
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
    if (m_digiTree) delete m_digiTree;
}

StatusCode digiRootReaderAlg::finalize()
{
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
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
