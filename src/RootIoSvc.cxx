/** 
* @file RootIoSvc.cxx
* @brief definition of the class RootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.27 2007/03/19 01:14:35 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/IObjManager.h"
#include "GaudiKernel/IToolFactory.h"
#include "GaudiKernel/IAlgManager.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"

#include "commonData.h"


#include <vector>
#include <algorithm>

/** 
* \class RootIoSvc
*
* \brief Service that implements the IRunable interface, to control the event loop.
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
*/

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IRunable.h"
#include "GaudiKernel/Property.h"
#include "RootIo/IRootIoSvc.h"
#include "TSystem.h"
#include "TFile.h"
#include "facilities/Util.h"

#include "RootConvert/Utilities/RootReaderUtil.h"

#include "RootInputDesc.h"

//forward declarations
template <class TYPE> class SvcFactory;
class IAppMgrUI;


class RootIoSvc : 
    virtual public Service, 
    virtual public IIncidentListener,
    virtual public IRootIoSvc,
    virtual public IRunable
{  
public:
    

    /// for the IRunnable interfce
    virtual StatusCode run();

    //------------------------------------------------------------------
    //  stuff required by a Service
    
    /// perform initializations for this service. 
    virtual StatusCode initialize ();
    
    /// perform the finalization, as required for a service.
    virtual StatusCode finalize ();
    
    /// Query interface
    virtual StatusCode queryInterface( const InterfaceID& riid, void** ppvUnknown );

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    
	
    // Old way, to be deleted
    virtual bool setRootFile(const char* mc, const char* digi, const char* rec, const char* gcr="");
    virtual std::string getMcFile() const { return m_mcFile; };
    virtual std::string getDigiFile() const { return m_digiFile; };
    virtual std::string getReconFile() const { return m_reconFile; };
    virtual std::string getGcrFile() const { return m_gcrFile; };
    virtual bool fileChange() const { return m_fileChange; };

    // New way
    virtual bool setFiles(const std::string& type, const StringArrayProperty& fileName);
    virtual StringArrayProperty getFileNameList(const std::string& type) const;
    virtual StatusCode registerIoAlgorithm(const std::string& type, 
                                           const std::string& tree, 
                                           const std::string& branch,
                                           const StringArrayProperty& fileList);
    
    virtual TObject* getNextEvent(const std::string& type);

    virtual Long64_t getEvtMax() { return m_evtMax; };

    virtual void setRootEvtMax(Long64_t max);

    virtual void setRootTimeMax(unsigned int max);

    virtual void registerRootTree(TChain *ch);

    virtual bool setIndex(Long64_t i);
    virtual void setActualIndex(Long64_t i) { m_index = i; };
    virtual Long64_t index() { return m_index; };

    virtual bool setRunEventPair(std::pair<int,int> ids);
    virtual std::pair<int,int> runEventPair() { return m_runEventPair; };

    virtual bool useIndex() { return m_useIndex; };
    virtual bool useRunEventPair() { return m_useRunEventPair; };

    virtual int getAutoSaveInterval() { return m_autoSaveInterval; };


protected: 
    
    /// Standard Constructor
    RootIoSvc ( const std::string& name, ISvcLocator* al );
    
    /// destructor
    virtual ~RootIoSvc();
    
private:

    bool                 fileChange(const std::string& type)    const;
    const RootInputDesc* getRootInputDesc(const std::string& type) const;
    RootInputDesc*       getRootInputDesc(const std::string& type); 

    void beginEvent();
    void endEvent();
    void loopStatus(int eventNumber, double currentTime, MsgStream& log);	
    /// Allow SvcFactory to instantiate the service.
    friend class SvcFactory<RootIoSvc>;
   
    /// Reference to application manager UI
    IAppMgrUI*    m_appMgrUI;
    IntegerProperty m_evtMax;
    IntegerProperty m_autoSaveInterval;
    UnsignedLongProperty m_treeSize;

    // starting and ending times for orbital simulation
    DoubleProperty m_startTime;
    DoubleProperty m_endTime;

    unsigned int m_rootEvtMax;
    Long64_t m_index;
    UnsignedLongProperty m_startIndex;
    std::pair<int, int> m_runEventPair;
    std::vector<TChain *> m_chainCol;

    commonData m_common;
    UInt_t m_objectNumber;

    bool m_useIndex, m_useRunEventPair;

    std::string m_mcFile, m_digiFile, m_reconFile, m_gcrFile;
    bool m_fileChange;
    int m_updated; // counter to check that algs have updated input root file

    std::map<std::string, RootInputDesc*> m_rootIoMap;
};

// declare the service factories for the RootIoSvc
static SvcFactory<RootIoSvc> a_factory;
const ISvcFactory& RootIoSvcFactory = a_factory;

// ------------------------------------------------
// Implementation of the RootIoSvc class
// ------------------------------------------------
/// Standard Constructor
RootIoSvc::RootIoSvc(const std::string& name,ISvcLocator* svc)
: Service(name,svc)
{
    
    declareProperty("EvtMax"     , m_evtMax=0);

    // disable for now
    //declareProperty("StartTime"   , m_startTime=0);
    //declareProperty("EndTime",      m_endTime=0);
    m_endTime = 0;
    m_startTime = 0;

    declareProperty("AutoSaveInterval", m_autoSaveInterval=1000);
    declareProperty("StartingIndex", m_startIndex=0);
    // limited by the size of an unsigned int, expecting MB
    declareProperty("MaxTreeSize", m_treeSize=0);

    m_index = 0;
    m_rootEvtMax = 0;
    m_runEventPair = std::pair<int,int>(-1,-1);
    m_chainCol.clear();
    m_useIndex = false;
    m_useRunEventPair = false;
    m_fileChange = false;
    m_updated = 0;

    m_rootIoMap.clear();
}


/// Standard Destructor
RootIoSvc::~RootIoSvc()  
{
    m_chainCol.clear();
}

StatusCode RootIoSvc::initialize () 
{   
    StatusCode  status =  Service::initialize ();
    
    // bind all of the properties for this service
    setProperties ();
    
    m_appMgrUI = 0 ;
    status = serviceLocator()->queryInterface(IID_IAppMgrUI, (void**)&m_appMgrUI);
    
    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);

    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    // Tell ROOT to reset signals to their default behavior
    gSystem->ResetSignal(kSigBus); 
    gSystem->ResetSignal(kSigSegmentationViolation); 
    gSystem->ResetSignal(kSigIllegalInstruction); 
    gSystem->ResetSignal(kSigFloatingException);  

    if (m_treeSize > 0) {
        // JO parameter in MB - need to convert to bytes
        Long64_t maxTreeSize = m_treeSize * 1000000;
        TTree::SetMaxTreeSize(maxTreeSize);
    } else if (m_treeSize == 0) {
        // 500 GB default
        //Long64_t maxTreeSize = 25000000000;
        Long64_t maxTreeSize = 500000000000;
        TTree::SetMaxTreeSize(maxTreeSize);
    }


    m_index = m_startIndex;
    if (m_startIndex >= 0) m_useIndex = true;

    return StatusCode::SUCCESS;
}

RootInputDesc*  RootIoSvc::getRootInputDesc(const std::string& type)
{
    RootInputDesc* rootIo = 0;

    // Look up this type in the map
    std::map<std::string, RootInputDesc*>::iterator typeItr = m_rootIoMap.find(type);

    // If one exists then get the pointer to the object
    if (typeItr != m_rootIoMap.end()) rootIo = typeItr->second;

    return rootIo;
}

const RootInputDesc* RootIoSvc::getRootInputDesc(const std::string& type) const
{
    const RootInputDesc* rootIo = 0;

    // Look up this type in the map
    std::map<std::string, RootInputDesc*>::const_iterator typeItr = m_rootIoMap.find(type);

    // If one exists then get the pointer to the object
    if (typeItr != m_rootIoMap.end()) rootIo = typeItr->second;

    return rootIo;
}

bool RootIoSvc::setFiles(const std::string& type, const StringArrayProperty& fileName)
{
    bool fileSet = false;

    RootInputDesc* rootIo = getRootInputDesc(type);

    if (rootIo)
    {
        rootIo->setFileList(fileName);

        fileSet = true;
    }
    else
    {
        // unregistered type so is an error
        fileSet = false;
    }

    return fileSet;
}

StringArrayProperty RootIoSvc::getFileNameList(const std::string& type) const
{
    StringArrayProperty result;

    const RootInputDesc* rootIo = getRootInputDesc(type);

    if (rootIo) result = rootIo->getFileList();

    return result;
}

StatusCode RootIoSvc::registerIoAlgorithm(const std::string& type, 
                                          const std::string& tree, 
                                          const std::string& branch,
                                          const StringArrayProperty& fileList)
{
    StatusCode sc = StatusCode::SUCCESS;

    // Create a new RootInputDesc object for this algorithm
    RootInputDesc* rootIo = new RootInputDesc(fileList, tree, branch);

    // Store the pointer to this in our map
    m_rootIoMap[type] = rootIo;

    // Register with RootIoSvc
    setRootEvtMax(rootIo->getNumEvents());
    registerRootTree(rootIo->getTChain());

    return sc;
}

TObject* RootIoSvc::getNextEvent(const std::string& type)
{
    MsgStream log( msgSvc(), name() );

    RootInputDesc* rootIo = getRootInputDesc(type);
    TObject*    pData  = 0;

    if (rootIo)
    {
        // Clear the event first
        //rootIo->clearEvent();

        // Read the event depending upon the mode
        if (useIndex())
        {
            pData = rootIo->getEvent(index());
        }
        // Otherwise we get by run/event number
        else
        {
            std::pair<int,int> runEvtPair = runEventPair();
            int                runNum     = runEvtPair.first;
            int                evtNum     = runEvtPair.second;

            pData = rootIo->getEvent(runNum, evtNum);
        }
    }

    return pData;
}


// finalize
StatusCode RootIoSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    return status;
}

/// Query interface
StatusCode RootIoSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IID_IRootIoSvc.versionMatch(riid) )  {
        *ppvInterface = (IRootIoSvc*)this;
    }else if (IID_IRunable.versionMatch(riid) ) {
      *ppvInterface = (IRunable*)this;
	} else if (IID_IIncidentListener.versionMatch(riid) ) {
		*ppvInterface = (IIncidentListener*)this;
	} else  {
        return Service::queryInterface(riid, ppvInterface);
    }

    addRef();
    return SUCCESS;
}


bool RootIoSvc::setRootFile(const char *mc, const char *digi, const char *rec, const char *gcr) 
{
    bool success = false;

    // Make local copies of the requested strings
    std::string mcFile    = mc;
    std::string digiFile  = digi;
    std::string reconFile = rec;
    std::string gcrFile   = gcr;

    // Expand to get rid of any environment variables
    facilities::Util::expandEnvVar(&mcFile);
    facilities::Util::expandEnvVar(&digiFile);
    facilities::Util::expandEnvVar(&reconFile);
    facilities::Util::expandEnvVar(&gcrFile);

    // at least one string must be non-null
    if (mcFile.empty() && digiFile.empty() && reconFile.empty() && gcrFile.empty())
        return success;

    // Check that these files exist
    // blank, skip - since that means we just won't read from that type of file
    if (!mcFile.empty()) 
    {
        if (!RootPersistence::fileExists(mcFile)) return success;

        // Change files in our RootInputDesc
        StringArrayProperty fileList;
        std::vector<std::string> fileVec;
        std::string type = "MC";
        fileVec.push_back("mc.root");
        fileList.setValue(fileVec);

        RootInputDesc* rootIo = getRootInputDesc(type);

        rootIo->setFileList(fileList);
    }

    if (!digiFile.empty()) 
    {
        if (!RootPersistence::fileExists(digiFile)) return success;
    }

    if (!reconFile.empty()) 
    {
        if (!RootPersistence::fileExists(reconFile)) return success;
        else 
        {
            TFile f(reconFile.c_str());
            if (f.IsOpen()) AcdRecon::fixAcdStreamer(f.GetVersion());
        }
    }

    if (!gcrFile.empty()) 
    {
        if (!RootPersistence::fileExists(gcrFile)) return success;
    }

    m_chainCol.clear(); // clear out TTrees from old files

    m_mcFile     = mc;
    m_digiFile   = digi;
    m_reconFile  = rec;
    m_gcrFile    = gcr;
    m_fileChange = true;

    success = true;
    return success;
}

void RootIoSvc::setRootEvtMax(Long64_t max) {
    // Purpose and Method:  Allow users of the RootIoSvc to specify the number
    //  of events found in their ROOT files
    if (m_rootEvtMax == 0) {
        m_rootEvtMax = max;
        return;
    } 
    
    if (m_rootEvtMax > max) m_rootEvtMax = max;
}

void RootIoSvc::setRootTimeMax(unsigned int max) {
    // Not yet used
    return;
}


void RootIoSvc::registerRootTree(TChain *ch) {
    m_chainCol.push_back(ch);
    ++m_updated;
}

bool RootIoSvc::setIndex(Long64_t i) {
     if (i < 0) return false;
     std::vector<TChain*>::iterator it;
     for(it = m_chainCol.begin(); it != m_chainCol.end(); it++) {
       if (i >= (*it)->GetEntries()) return false;
     }
     m_index = i;
     m_runEventPair = std::pair<int, int>(-1,-1);
     m_useIndex = true;
     m_useRunEventPair = false;
     return true;
}


bool RootIoSvc::setRunEventPair(std::pair<int, int> ids) {
    // If we just changed ROOT files, we may not have run the reader algs
    // yet, so the TTrees are not set to read, so we temporarily load the
    // the TTrees so we can check to see if the requested run/event pair
    // exists.  
    TFile *mc=0, *digi=0, *rec=0, *gcr=0;
    TChain *mcChain, *digiChain, *recChain, *gcrChain;

    if ((m_fileChange) && (m_chainCol.size() == 0)) {
        if (!m_mcFile.empty()) {
           mc = new TFile(m_mcFile.c_str(), "READ");
           if (!mc->IsOpen()) return false;
           mcChain = (TChain*)mc->Get("Mc");
           if (!mcChain) return false;
           if (!mcChain->GetTreeIndex()) 
               mcChain->BuildIndex("m_runId", "m_eventId");
           registerRootTree(mcChain);
        }
        if (!m_digiFile.empty()) {
           digi = new TFile(m_digiFile.c_str(), "READ");
           if (!digi->IsOpen()) return false;
           digiChain = (TChain*)digi->Get("Digi");
           if (!digiChain) return false;
           if (!digiChain->GetTreeIndex()) 
               digiChain->BuildIndex("m_runId", "m_eventId");
           registerRootTree(digiChain);
        }
        if (!m_reconFile.empty()) {
           rec = new TFile(m_reconFile.c_str(), "READ");
           if (!rec->IsOpen()) return false;
           recChain = (TChain*)rec->Get("Recon");
           if (!recChain) return false;
           if (!recChain->GetTreeIndex()) 
               recChain->BuildIndex("m_runId", "m_eventId");
           registerRootTree(recChain);
        }

	if (!m_gcrFile.empty()) {
           gcr = new TFile(m_gcrFile.c_str(), "READ");
           if (!rec->IsOpen()) return false;
           gcrChain = (TChain*)gcr->Get("GcrSelect");
           if (!gcrChain) return false;
           if (!gcrChain->GetTreeIndex()) 
               gcrChain->BuildIndex("m_runId", "m_eventId");
           registerRootTree(gcrChain);
        }


    }

    std::vector<TChain*>::iterator it;
    for(it = m_chainCol.begin(); it != m_chainCol.end(); it++) {
        Long64_t readInd = (*it)->GetEntryNumberWithIndex(ids.first, ids.second);
        if ( (readInd < 0) || (readInd >= (*it)->GetEntries()) ) {
          if (m_fileChange) {
              m_chainCol.clear();
              if (mc) delete mc;
              if (digi) delete digi;
              if (rec) delete rec;
              if (gcr) delete gcr;
          }
          return false;
        }
    }
    m_runEventPair = ids;
    m_useIndex = false;
    m_useRunEventPair = true;
    if (m_fileChange) {
        m_chainCol.clear();
        if (mc) delete mc;
        if (digi) delete digi;
        if (rec) delete rec;
        if (gcr) delete gcr;
    }
    return true;
}

// handle "incidents"
void RootIoSvc::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
}


void RootIoSvc::beginEvent() // should be called at the beginning of an event
{ 
    m_objectNumber = TProcessID::GetObjectCount();
}

void RootIoSvc::endEvent()  // must be called at the end of an event to update, allow pause
{        
    m_useIndex = false;
    m_useRunEventPair = false;

    //temp
    m_useIndex = true;
    m_index++;

    m_runEventPair = std::pair<int, int>(-1,-1);

    // clear out the static maps
    m_common.clear();

    // reset object in order to avoid memleak
    TProcessID::SetObjectCount(m_objectNumber);

    // assuming all algs have gotten the new files by now
    if (m_updated > 0) {
        m_fileChange = false;
        m_updated = 0;
    }
}

StatusCode RootIoSvc::run(){
    // Purpose and Method:  Control the event loop

    StatusCode status = StatusCode::FAILURE;
    MsgStream log( msgSvc(), name() );

    if ( 0 == m_appMgrUI )  return status; 

    IProperty* propMgr=0;
    status = serviceLocator()->service("ApplicationMgr", propMgr );
    if( status.isFailure()) {
        log << MSG::ERROR << "Unable to locate PropertyManager Service" << endreq;
        return status;
    }
    
    IntegerProperty evtMax("EvtMax",0);
    status = propMgr->getProperty( &evtMax );
    if (status.isFailure()) return status;
 
    // Determine if the min number of ROOT events is less than the
    // requested number of events in the jobOptions file
    IntegerProperty rootEvtMax("EvtMax", m_rootEvtMax);
    if ( static_cast<int>(rootEvtMax-m_startIndex) < evtMax) {
       evtMax = rootEvtMax - m_startIndex;
       setProperty(evtMax);
    } else setProperty(evtMax);

    // now find the top alg so we can monitor its error count
    //
    IAlgManager* theAlgMgr =0 ;
    status = serviceLocator( )->getService( "ApplicationMgr",
        IID_IAlgManager,
        (IInterface*&)theAlgMgr );
    IAlgorithm* theIAlg;
    Algorithm*  theAlgorithm=0;
    IntegerProperty errorProperty("ErrorCount",0);
    
    status = theAlgMgr->getAlgorithm( "Top", theIAlg );
    if ( status.isSuccess( ) ) {
        try{
            theAlgorithm = dynamic_cast<Algorithm*>(theIAlg);
        } catch(...){
            status = StatusCode::FAILURE;
        }
    }
    if ( status.isFailure( ) ) {
        log << MSG::WARNING << "Could not find algorithm 'Top'; will not monitor errors" << endreq;
    }
    
    
    // loop over the events

    int eventNumber= 0;
    double currentTime=m_startTime;
    
    { bool noend=true;
    log << MSG::INFO << "Runable interface starting event loop as :" ; 
    if( m_evtMax>0)  { log << " MaxEvt = " << m_evtMax; noend=false;  }
    if( m_endTime>0) { log << " EndTime= " << m_endTime; noend=false; }
    log << endreq;
    
    if(noend) { 
        log << MSG::WARNING << "No end condition specified: will not process any events!" << endreq; 
    }
    }
    // Not yet using time as a control on the event loop for ROOT
    while( m_evtMax>0 && eventNumber < m_evtMax
        || m_endTime>0 && currentTime< m_endTime ) {
        
        loopStatus( eventNumber, currentTime, log);

        status =  m_appMgrUI->nextEvent(1); // currently, always success
        
        // the single event may have created a failure. Check the ErrorCount propery of the Top alg.
        if( theAlgorithm !=0) theAlgorithm->getProperty(&errorProperty);
        if( status.isFailure() || errorProperty.value() > 0){
            status = StatusCode::FAILURE;
        }
        
        if( status.isFailure()) break;
        //if(flux!=0){
         //   currentTime = flux->gpsTime();
       // }
        eventNumber++;
    }
    if( status.isFailure()){
        log << MSG::ERROR << "Terminating RootIoSvc loop due to error" << endreq;
        
    }else if( m_endTime>0 && currentTime >= m_endTime ) {
        log << MSG::INFO << "Loop terminated by time " << endreq;
    }else {
        log << MSG::INFO << "Processing loop terminated by event count" << endreq;
    }
    return status;
}

void RootIoSvc::loopStatus(int eventNumber, double currentTime, MsgStream& log)
// purpose: print periodic messages
{
    static int last_fraction=0;

    double efrac =   (m_evtMax>0? 100.*eventNumber/m_evtMax: 0.0), 
        tfrac =   (m_endTime>0? 100.*(currentTime-m_startTime)/(m_endTime-m_startTime) : 0.0) ;

    int percent_complete= static_cast<int>(  std::max( efrac, tfrac)  );
    if( percent_complete!=last_fraction){
        last_fraction=percent_complete;
        if( percent_complete<10 || percent_complete%10 ==0){
            log << MSG::INFO <<   percent_complete << "% complete: "
                << " event "<< eventNumber<<",  time "<< currentTime << endreq;
        }

    }

}
