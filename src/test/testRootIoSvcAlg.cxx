#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "RootIo/IRootIoSvc.h"

/** @class testRootIoSvcAlg
 * @brief Takes data from the TDS to test reading from ROOT files
 *
 * @author Heather Kelly
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/testRootIoSvcAlg.cxx,v 1.2 2004/01/12 19:28:54 heather Exp $
 */

class testRootIoSvcAlg : public Algorithm
{	
public:
    
    testRootIoSvcAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
        
private:
    IRootIoSvc*   m_rootIoSvc;
 
};

static const AlgFactory<testRootIoSvcAlg>  Factory;
const IAlgFactory& testRootIoSvcAlgFactory = Factory;


testRootIoSvcAlg::testRootIoSvcAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode testRootIoSvcAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();

    if ( service("RootIoSvc", m_rootIoSvc).isFailure() ){
        log << MSG::ERROR << "Couldn't find the RootIoSvc!" << endreq;
        return StatusCode::FAILURE;
    }    
	
	return sc;
    
}

StatusCode testRootIoSvcAlg::execute()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	static int flag = 0;
	
	if (flag % 2)
		m_rootIoSvc->setIndex(flag);
    else {
		m_rootIoSvc->setRunEventPair(std::pair<int,int>(10,48));
        log << MSG::INFO << "Requesting run/event (10,48) randomly" << endreq;
    }
        //m_rootIoSvc->setIndex(flag);

	flag++;


	return sc;
}


StatusCode testRootIoSvcAlg::finalize()
{    
    StatusCode sc = StatusCode::SUCCESS;
    return sc;
}

