#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "RootIo/IRootIoSvc.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

/** @class checkTds
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/checkTds.cxx,v 1.2.656.1 2010/08/31 03:08:51 heather Exp $
 */

class checkTds : public Algorithm
{	
public:
    
    checkTds(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:
 
};

//static const AlgFactory<checkTds>  Factory;
//const IAlgFactory& checkTdsFactory = Factory;
DECLARE_ALGORITHM_FACTORY(checkTds);


checkTds::checkTds(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode checkTds::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();

    return sc;
    
}

StatusCode checkTds::execute()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    log << MSG::DEBUG << "Current Event" << endreq;
    SmartDataPtr<Event::EventHeader> header(eventSvc(), EventModel::EventHeader);
    log << MSG::DEBUG << "run: " << header->run() << " event: " 
        << header->event() << endreq;

	return sc;
}


StatusCode checkTds::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

