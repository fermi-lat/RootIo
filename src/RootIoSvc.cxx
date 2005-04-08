/** 
* @file RootIoSvc.cxx
* @brief definition of the class RootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.16 2005/04/06 20:07:43 heather Exp $
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

//#include "CLHEP/Random/Random.h"

#include <vector>
#include <algorithm>

/** 
* \class RootIoSvc
*
* \brief Service that implements the IRunable interface, to control the event loop.
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.16 2005/04/06 20:07:43 heather Exp $
*/

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IRunable.h"
#include "GaudiKernel/Property.h"
#include "RootIo/IRootIoSvc.h"
#include "TSystem.h"

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
    virtual StatusCode queryInterface( const IID& riid, void** ppvUnknown );

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    
	
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

    void beginEvent();
    void endEvent();
	
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

    bool m_useIndex, m_useRunEventPair;

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
        // 25 GB default
        Long64_t maxTreeSize = 25000000000;
        TTree::SetMaxTreeSize(maxTreeSize);
    }


    m_index = m_startIndex;
    if (m_startIndex > 0) m_useIndex = true;

    return StatusCode::SUCCESS;
}


// finalize
StatusCode RootIoSvc::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    return status;
}

/// Query interface
StatusCode RootIoSvc::queryInterface(const IID& riid, void** ppvInterface)  {
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
    std::vector<TChain*>::iterator it;
    for(it = m_chainCol.begin(); it != m_chainCol.end(); it++) {
        Long64_t readInd = (*it)->GetEntryNumberWithIndex(ids.first, ids.second);
        if ( (readInd < 0) || (readInd >= (*it)->GetEntries()) ) return false;
    }
    m_runEventPair = ids;
    m_useIndex = false;
    m_useRunEventPair = true;
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
}

void RootIoSvc::endEvent()  // must be called at the end of an event to update, allow pause
{        
    m_useIndex = false;
    m_useRunEventPair = false;

    m_runEventPair = std::pair<int, int>(-1,-1);
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
    if ((rootEvtMax-m_startIndex) < evtMax) {
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

