/** 
 * @file ntupleReaderSvc.cxx
 * @brief declare, implement the class ntupleReaderSvc
 *
 * Special service that directly writes ROOT tuples
 * It also allows multiple TTree's in the root file: see the addItem (by pointer) member function.
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/Attic/ntupleReaderSvc.cxx,v 1.1.2.1 2008/03/22 04:20:15 heather Exp $
 */

#include "GaudiKernel/Service.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Property.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/DataStoreItem.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"


#include "RootIo/INtupleReaderSvc.h"
#include "facilities/Util.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

// root includes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TLeaf.h"
#include "TVirtualIndex.h"

//#include <iomanip>
//#include <list>
#include <string>


class ntupleReaderSvc :  public Service, virtual public IIncidentListener,
        virtual public INtupleReaderSvc

{  
public:

    /// perform initializations for this service - required by Gaudi
    virtual StatusCode initialize ();

    /// clean up after processing all events - required by Gaudi
    virtual StatusCode finalize ();

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    

    /// Query interface - required of all Gaudi services
    virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvUnknown );

    // Load the next event in the chain
    virtual bool nextEvent(bool checkTds=true);

    // grab the requested event based on the input index
    virtual bool getEvent(Long64_t idx, bool checkTds=true);

    // Load the requested event based on run and event ids
    virtual bool getEvent(unsigned int run, unsigned int evt, bool checkTds=true);

    virtual long long numEvents() const { return (m_numEvents); }
    
    virtual std::string getItemType(const std::string& itemName) const; 

  //! acccess to an item (interface to the leaf)
   virtual std::string getItem(const std::string& itemName, void*& pointer)const;

   virtual bool ntupleReaderSvc::addFile(const std::string & fileName, const std::string & treeName);

   virtual bool clearChain();

private:
    /// Allow only SvcFactory to instantiate the service.
    friend class SvcFactory<ntupleReaderSvc>;

    ntupleReaderSvc ( const std::string& name, ISvcLocator* al );    

    bool setLeaves();

    // Looks at the run/event ids on the TDS and will compare them to the loaded event 
    // in the chain
    int compareTdsIds();

    /// routine to be called at the beginning of an event
    void beginEvent();
    /// routine that is called when we reach the end of an event
    StatusCode endEvent();

    /// List of input ntuples files
    StringArrayProperty m_fileList;
    /// Name of tree
    StringProperty m_treeName;

    TChain *m_chain;
    TVirtualIndex *m_runEvtIndex;    /// Save the RunId/EventId Index in case other indices are in use such as CompositeEventList


    // Index of currently loaded event
    Long64_t m_curEvtInd;

    // Total number of events in the chain
    Long64_t m_numEvents;

    // Ptr to the EvtRunId in the TTree
    TLeaf * m_runLeaf;
    // Ptr to the EvtEventId in the TTree
    TLeaf * m_evtLeaf;

    IDataProviderSvc* m_eventSvc;

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declare the service factories for the ntupleWriterSvc
static SvcFactory<ntupleReaderSvc> a_factory;
const ISvcFactory& ntupleReaderSvcFactory = a_factory;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//         Implementation of ntupleReaderSvc methods
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ntupleReaderSvc::ntupleReaderSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc), m_chain(0), m_runEvtIndex(0), m_curEvtInd(-1), m_numEvents(0), m_eventSvc(0)
{
    // declare the properties and set defaults
    declareProperty("fileList",  m_fileList);
    declareProperty("treename", m_treeName="MeritTuple");

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode ntupleReaderSvc::initialize () 
{
    StatusCode  status =  Service::initialize ();

    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigIllegalInstruction);
    gSystem->ResetSignal(kSigFloatingException); 

    // bind all of the properties for this service
    setProperties ();

    // open the message log
    MsgStream log( msgSvc(), name() );

    // use the incident service to register "end" events
    IIncidentSvc* incsvc = 0;
    status = service("IncidentSvc", incsvc, true);
    if( status.isFailure() ) {
        log << MSG::ERROR << "Failed to load IncidentSvc" << endreq;
        return status;
    }

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    if( serviceLocator() ) {
        status = serviceLocator()->service("EventDataSvc", m_eventSvc, true );
    }
    if(status.isFailure())
    {
        log << MSG::INFO << "Could not find EventDataSvc" << endreq;
        m_eventSvc = 0;
    }

    // -- create primary root file---
    TDirectory* curdir = gDirectory;
    m_chain = new TChain(m_treeName.value().c_str());
    std::vector<std::string>::const_iterator listIt;
    for (listIt = m_fileList.value().begin(); listIt != m_fileList.value().end(); listIt++) {
        m_chain->AddFile((*listIt).c_str());
    }

    m_numEvents = m_chain->GetEntries();
    if (!setLeaves())
        log << MSG::INFO << "Failed to locate either EvtRun or EvtEventId in the tuple" << endreq;

    curdir->cd(); // restore previous directory

    return status;
}

bool ntupleReaderSvc::setLeaves() {
    if( (m_runLeaf = m_chain->GetLeaf("EvtRun")) == 0) {
        m_runLeaf = 0;
        return false;
    }
    if( (m_evtLeaf = m_chain->GetLeaf("EvtEventId")) == 0) {
        m_evtLeaf = 0;
        return false;
    }
    return true;
}

bool ntupleReaderSvc::nextEvent(bool checkTds) {
    if(m_curEvtInd < m_numEvents) {
        int numBytes = m_chain->GetEntry(m_curEvtInd++);
        if (numBytes <= 0) return false;
        if (checkTds) {
            if (compareTdsIds() == -1) return false;
        }
        return true;
    }
    return false;
}

bool ntupleReaderSvc::getEvent(Long64_t idx, bool checkTds) {

    if (idx >= 0 && idx < m_numEvents)
    {
        int numBytes = m_chain->GetEntry(idx);
        if (numBytes <= 0) return false;
        if(checkTds) {
            if (compareTdsIds() == -1)
                return false;
        }
        return true;
    }
    return false;
}

bool ntupleReaderSvc::getEvent(unsigned int run, unsigned int evt, bool checkTds) {
 if (!m_runEvtIndex) 
   {
    m_chain->BuildIndex("EvtRunId", "EvtEventId") ;
    m_runEvtIndex = m_chain->GetTreeIndex();
   }
   m_chain->SetTreeIndex(m_runEvtIndex);
   int numBytes = m_chain->GetEntryWithIndex(run,evt) ;
   if (numBytes <= 0) return false;
   if (checkTds) {
       if (compareTdsIds() == -1)
           return false;
   }
   return true;
}

int ntupleReaderSvc::compareTdsIds() {
    if ( (!m_runLeaf) || (!m_evtLeaf)) return -4;
     unsigned int runId = static_cast<unsigned int>(m_runLeaf->GetValue());
     unsigned int evtId = static_cast<unsigned int>(m_evtLeaf->GetValue());
    // Retun -2, if no event Data Svc available
    if (!m_eventSvc) return -2;
    SmartDataPtr<Event::EventHeader> evt(m_eventSvc, EventModel::EventHeader);
    // return -3 if the event header is not available on the TDS
    if (!evt) return -3;

    unsigned int eventIdTds = evt->event();
    unsigned int runIdTds = evt->run();

    // return -1 if the run or event ids do not match
    if ((eventIdTds != evtId) || (runIdTds != runId)) {
        MsgStream log(msgSvc(),name());
        log << MSG::INFO << "run or event id in tuple and TDS do not match TDS: "
            << runIdTds << " " << eventIdTds << " Tuple: " << runId << " " << evtId << endreq;
        return -1;
    }

    // success
    return 0;
}

std::string ntupleReaderSvc::getItemType(const std::string& itemName) const {
    TDirectory *saveDir = gDirectory;
    std::string type_name("");
    if (!m_chain) return(type_name);
    TLeaf *leaf = m_chain->GetLeaf(itemName.c_str());
    if (!leaf) return(type_name);
    type_name = leaf->GetTypeName();
    saveDir->cd();
    return type_name;
}


std::string ntupleReaderSvc::getItem(const std::string& itemName, void*& pval)const
{
    MsgStream log(msgSvc(),name());
    StatusCode status = StatusCode::SUCCESS;
    TDirectory *saveDir = gDirectory;

    m_chain->GetCurrentFile()->cd();

    if( itemName.empty()){
        // assume this is a request for the tree
        pval = (void *)m_chain;
        saveDir->cd();
        return itemName;
    }
    TLeaf* leaf = m_chain->GetLeaf(itemName.c_str());
    if( leaf==0){
        throw std::invalid_argument(std::string("ntupleReaderSvc::getItem: did not find tuple or leaf: ")+ itemName);
    }
    pval = leaf->GetValuePointer();
    std::string type_name(leaf->GetTypeName());
    saveDir->cd();
    return type_name;
}

bool ntupleReaderSvc::addFile(const std::string & fileName, const std::string & treeName) {
    int status = 0;
    if (m_chain)
        status = m_chain->AddFile(fileName.c_str(),TChain::kBigNumber,treeName.c_str());
    else {
        m_chain = new TChain(treeName.c_str());
        if (!m_chain) return false;
        status = m_chain->AddFile(fileName.c_str());
    }
    // TChain::AddFile returns 1 if success, and 0 if failure
    if (status == 0) return false;
    return true;
}

bool ntupleReaderSvc::clearChain() {
    if (m_chain) {
        delete m_chain;
        m_chain = 0;
        m_runEvtIndex = 0;
        m_numEvents = 0;
        m_curEvtInd = -1;
    }
    return true;
}


void ntupleReaderSvc::handle(const Incident &inc)
{
    // Purpose and Method:  This routine is called when an "incident"
    //   occurs.  This method determines what action the ntupleReaderSvc
    //   will take in response to a particular event.  Currently, we handle
    //   BeginEvent and EndEvent events.

    if(inc.type()=="BeginEvent") beginEvent();
    if(inc.type()=="EndEvent") endEvent();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ntupleReaderSvc::beginEvent()
{

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode ntupleReaderSvc::endEvent()
    // must be called at the end of an event to update, allow pause
{         
//    MsgStream log(msgSvc(),name());
    StatusCode sc = SUCCESS;
//    TDirectory *saveDir = gDirectory;

//    saveDir->cd();
    return sc;

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StatusCode ntupleReaderSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IID_INtupleReaderSvc.versionMatch(riid) )  {
        *ppvInterface = (INtupleReaderSvc*)this;
    }
    else  {
        return Service::queryInterface(riid, ppvInterface);
    }
    addRef();
    return SUCCESS;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

StatusCode ntupleReaderSvc::finalize ()
{
    // open the message log
    MsgStream log( msgSvc(), name() );

    if (m_chain) delete m_chain;
    m_chain = 0;

    return StatusCode::SUCCESS;
}

