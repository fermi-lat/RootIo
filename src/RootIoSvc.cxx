/** 
* @file RootIoSvc.cxx
* @brief definition of the class RootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.1 2003/06/02 01:53:19 heather Exp $
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

#include "CLHEP/Random/Random.h"


#include <algorithm>
/** 
* \class RootIoSvc
*
* \brief Service that implements the IRunable interface, to control the event loop.
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.1 2003/06/02 01:53:19 heather Exp $
*/

// includes
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IRunable.h"
#include "GaudiKernel/Property.h"
#include "RootIo/IRootIoSvc.h"

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
	
	virtual int getEvtMax() { return m_evtMax; };

    virtual void setRootEvtMax(unsigned int max);

    virtual void setRootTimeMax(unsigned int max);

	virtual void setIndex(int i);
	virtual int index() { return m_index; };

	virtual void setRunEventPair(std::pair<int,int> ids);
	virtual std::pair<int,int> runEventPair() { return m_runEventPair; };

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

    // starting and ending times for orbital simulation
    DoubleProperty m_startTime;
    DoubleProperty m_endTime;

    unsigned int m_rootEvtMax;
	int m_index;
	std::pair<int, int> m_runEventPair;

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
    declareProperty("StartTime"   , m_startTime=0);
    declareProperty("EndTime",      m_endTime=0);
    m_rootEvtMax = 0;
	m_index = -1;
	m_runEventPair = std::pair<int,int>(-1,-1);
}


/// Standard Destructor
RootIoSvc::~RootIoSvc()  
{
}


// initialize
StatusCode RootIoSvc::initialize () 
{   
    StatusCode  status =  Service::initialize ();
    
    // bind all of the properties for this service
    setProperties ();
    
    // open the message log
    MsgStream log( msgSvc(), name() );
    
    status = serviceLocator()->queryInterface(IID_IAppMgrUI, (void**)&m_appMgrUI);
    
    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);

    if( status.isFailure() ) return status;

    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

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


void RootIoSvc::setRootEvtMax(unsigned int max) {
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

void RootIoSvc::setIndex(int i) {
	m_index = i;
	m_runEventPair = std::pair<int, int>(-1,-1);
}

void RootIoSvc::setRunEventPair(std::pair<int, int> ids) {
	m_runEventPair = ids;
	m_index=-1;
}
// handle "incidents"
void RootIoSvc::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
}


void RootIoSvc::beginEvent() // should be called at the beginning of an event
{ 
	setIndex(-1);
	setRunEventPair(std::pair<int,int>(-1,-1));
}

void RootIoSvc::endEvent()  // must be called at the end of an event to update, allow pause
{        
	setIndex(-1);
	setRunEventPair(std::pair<int,int>(-1,-1));
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
    if (rootEvtMax < evtMax) setProperty(rootEvtMax);
    else setProperty(evtMax);

    // now find the top alg so we can monitor its error count
    //
    IAlgManager* theAlgMgr;
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
        eventNumber ++;
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

