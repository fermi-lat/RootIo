#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "LdfEvent/EventSummaryData.h"

#include "facilities/Util.h"
#include "commonData.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TVector3.h"
#include "TMatrixD.h"

#include "gcrSelectRootData/GcrSelectEvent.h"
#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"

// low level converters

#include "RootConvert/GcrSelect/GcrSelectedXtalConvert.h"
#include "RootConvert/GcrSelect/GcrSelectValsConvert.h"

#include <cstdlib>

/** @class gcrSelectRootWriterAlg
* @brief Writes Recon TDS data to a persistent ROOT file.
*
* @author Heather Kelly and Tracy Usher
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/gcrSelectRootWriterAlg.cxx,v 1.8 2008/04/23 16:41:17 cohen Exp $
*/

class gcrSelectRootWriterAlg : public Algorithm
{	
public:
    
    gcrSelectRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
    
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
    
private:
    
    
    StatusCode writeGcrSelectEvent();

    StatusCode writeGcrSelect();
    
    /// These are the methods specific to filling the pieces of the CalRecon stuff    
    //CL: 08/22/06:
    void fillGcrSelectedXtal(GcrSelect *gcrSelect, Event::GcrSelectedXtalsCol* gcrSelectedXtalColTds); 
    void fillGcrSelectVals(GcrSelect *gcrSelect, Event::GcrSelectVals* gcrSelectValsTds); 
    
    /// Calls TTree::Fill for each event and clears m_mcEvt
    void writeEvent();
    
    /// Performs the final write to the ROOT file and closes
    void close();
    
    /// ROOT file pointer
    TFile *m_gcrSelectFile;
    /// ROOT tree pointer
    TTree *m_gcrTree;
    /// Top-level Monte Carlo ROOT object
    GcrSelectEvent *m_gcrSelEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the TTree in the ROOT file
    std::string m_treeName;
    /// ROOT split mode
    int m_splitMode;
    /// Buffer Size for the ROOT file
    int m_bufSize;
    /// Compression level for the ROOT file
    int m_compressionLevel;

    /// Keep track of relation between TDS objs and ROOT counterparts
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};


static const AlgFactory<gcrSelectRootWriterAlg>  Factory;
const IAlgFactory& gcrSelectRootWriterAlgFactory = Factory;

gcrSelectRootWriterAlg::gcrSelectRootWriterAlg(const std::string& name, 
                                       ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("gcrSelectRootFile",m_fileName="gcrSelect.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="GcrSelect");    
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool = 0 ;
}

StatusCode gcrSelectRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    StatusCode headersSc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }
    headersSc = m_headersTool->newGcrHeader() ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to create a new Gcr FileHeader" << endreq;
    }
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        m_rootIoSvc = 0;
        // RootIoSvc is now required for reading/writing, cannot continue without it
        return StatusCode::FAILURE;
    } 


    m_gcrTree = m_rootIoSvc->prepareRootOutput("gcr", m_fileName, m_treeName, 
        m_compressionLevel, "GLAST GcrSelect Data");

    m_gcrSelEvt = new GcrSelectEvent();
    m_rootIoSvc->setupBranch("gcr", "GcrSelectEvent", "GcrSelectEvent", &m_gcrSelEvt, m_bufSize, m_splitMode);
    
    return sc;
    
}

StatusCode gcrSelectRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log << MSG::DEBUG << "gcrSelectRootWriterAlg::execute BEGIN" << endreq;
        
    if(!m_gcrTree->GetCurrentFile()->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    
    if (m_gcrSelEvt) m_gcrSelEvt->Clear();

    sc = writeGcrSelectEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write GcrSelectEvent" << endreq;
        return sc;
    }
    
    sc = writeGcrSelect();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write GcrSelect" << endreq;
        return sc;
    }
       
    writeEvent();
    
        log << MSG::DEBUG << "gcrSelectRootWriterAlg::execute BEGIN" << endreq;


    return sc;
}


StatusCode gcrSelectRootWriterAlg::writeGcrSelectEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the ReconEvent ROOT object
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);
    
    if (!evtTds) return sc;
    
    UInt_t evtId = evtTds->event();
    UInt_t runId = evtTds->run();
    
    log << MSG::DEBUG;
    if( log.isActive())evtTds->fillStream(log.stream());
    log << endreq;
    
    m_gcrSelEvt->initialize(evtId, runId, new GcrSelect());

    // For simulated data - this may not exist on the TDS and that is ok
    // no need to fail for that
    SmartDataPtr<LdfEvent::EventSummaryData> summaryTds(eventSvc(), "/Event/EventSummary");
    if (summaryTds) m_gcrSelEvt->initEventFlags(summaryTds->eventFlags());
    
    return sc;
}



StatusCode gcrSelectRootWriterAlg::writeGcrSelect() {
    // Purpose and Method:  Retrieve the Cal Recon data from the TDS 
    //  calls the two helper methods that do the work of filling the ROOT    
    //  version of the CalRecon object for this event. 
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
   // log << MSG::INFO << "gcrSelectRootWriterAlg::writeGcrSelect BEGIN" << endreq;

    
    // Get pointer to CalRecon part of ReconEvent   
    GcrSelect* gcrSelect = m_gcrSelEvt->getGcrSelect();   
    if (!gcrSelect) return StatusCode::FAILURE;   
    gcrSelect->initialize();
    
    //C.L: 08/22/06: Retrieve GcrSelectedXtal collection
       //log << MSG::INFO << "just before retrieving GcrSelectedXtal collection" << endreq;
    SmartDataPtr<Event::GcrSelectedXtalsCol> gcrSelectedXtalColTds(eventSvc(), EventModel::CalRecon::GcrSelectedXtalsCol);   
    
    if (gcrSelectedXtalColTds) 
      fillGcrSelectedXtal(gcrSelect, gcrSelectedXtalColTds); 

    SmartDataPtr<Event::GcrSelectVals> gcrSelectValsTds(eventSvc(), EventModel::CalRecon::GcrSelectVals);   
    if (gcrSelectValsTds) 
      fillGcrSelectVals(gcrSelect, gcrSelectValsTds); 
    else
      log << MSG::DEBUG << "no gcrSelectVals found in TDS, do not fill GcrSelectVals" << endreq;
    
   //log << MSG::INFO << "gcrSelectRootWriterAlg::writeGcrSelect END" << endreq;
    
    return sc;
}


void gcrSelectRootWriterAlg::fillGcrSelectedXtal(GcrSelect *gcrSelect, Event::GcrSelectedXtalsCol* gcrSelectedXtalColTds) {   
    // Purpose and Method:  Given the GcrXtal collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
    MsgStream log(msgSvc(), name());
    Event::GcrSelectedXtalsCol::const_iterator gcrSelXtalColIterTds;
    
	log << MSG::DEBUG << "gcrSelectRootWriterAlg::fillGcrSelectedXtal BEGIN, gcrSelectedXtalColTds->size()=" << gcrSelectedXtalColTds->size()<< endreq;   
    
    for (gcrSelXtalColIterTds = gcrSelectedXtalColTds->begin(); gcrSelXtalColIterTds != gcrSelectedXtalColTds->end(); gcrSelXtalColIterTds++) {
	
	GcrSelectedXtal* gcrSelectedXtalRoot = gcrSelect->addGcrSelectedXtal();
	Event::GcrSelectedXtal* gcrSelXtalTds = *gcrSelXtalColIterTds;
	RootPersistence::convert(**gcrSelXtalColIterTds,*gcrSelectedXtalRoot) ; 
	
	log << MSG::DEBUG << "ROOT@@@@@ \n gcrSelectedXtalRoot->getXtalId()=" 	
	<< ", rowEnergy=" << gcrSelectedXtalRoot-> getRawEnergy()
	<< ", pathLength=" << gcrSelectedXtalRoot-> getPathLength()
	<< ", corrEnergy=" << gcrSelectedXtalRoot-> getCorrEnergy()
	<< ", crossedFaces=" << gcrSelectedXtalRoot-> getCrossedFaces()
	<< "\n" << endreq;
	}
	
    log << MSG::DEBUG << "gcrSelectRootWriterAlg::fillGcrSelectedXtal END" << endreq;   

    
    return;   
}


void gcrSelectRootWriterAlg::fillGcrSelectVals(GcrSelect *gcrSelect, Event::GcrSelectVals* gcrSelectValsTds) {

    // Purpose and Method:  Given the GcrXtal collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
     MsgStream log(msgSvc(), name());
 
    log << MSG::VERBOSE << "gcrSelectRootWriterAlg::fillGcrSelectVals BEGIN" << endreq;   
    
    
    GcrSelectVals* gcrSelectValsRoot = new GcrSelectVals();
    
    RootPersistence::convert(*gcrSelectValsTds,*gcrSelectValsRoot);
    
    //log << MSG::INFO << "gcrSelectValsTds->getInferedZ()" << gcrSelectValsTds->getInferedZ()<< endreq;   
    
    //log << MSG::INFO << "gcrSelectValsRoot->getInferedZ()" << gcrSelectValsRoot->getInferedZ()<< endreq;   
    
    log << MSG::DEBUG << "gcrOBFStatusWord=" << std::hex << gcrSelectValsRoot->getGcrOBFStatusWord() << std::dec << endreq;
   
    gcrSelect->addGcrSelectVals(gcrSelectValsRoot);
    
  
   log << MSG::VERBOSE << "gcrSelectRootWriterAlg::fillGcrSelectVals END" << endreq;   
    
    return;   
 
}




void gcrSelectRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the GcrSelectEvent data for this event in the ROOT
    //    tree.  
    
    m_rootIoSvc->fillTree("gcr");
    
    return;
}

void gcrSelectRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
   //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.
    
    m_rootIoSvc->closeFile("gcr");

    return;
}

StatusCode gcrSelectRootWriterAlg::finalize()
{
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool->writeGcrHeader(m_gcrTree->GetCurrentFile()) ;
    
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();

    log << MSG::DEBUG;
    if( log.isActive()) log.stream() << "Finalized";
    log << endreq;

    return sc;
}

