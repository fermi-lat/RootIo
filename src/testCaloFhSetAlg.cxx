
#include "src/testCaloFhSetAlg.h"
#include "RootIo/FhTool.h"

#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"

#include <GaudiKernel/MsgStream.h>
#include <GaudiKernel/AlgFactory.h>
#include <GaudiKernel/IService.h>
#include <GaudiKernel/ISvcLocator.h>
#include <GaudiKernel/IDataProviderSvc.h>
#include <GaudiKernel/SmartDataPtr.h>

#include <cstdio>

static const AlgFactory<testCaloFhSetAlg> Factory ;
const IAlgFactory& testCaloFhSetAlgFactory = Factory ;

StatusCode testCaloFhSetAlg::finalize() {
    
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFhTool * headersTool ;
    sc = toolSvc()->retrieveTool("FhTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }
    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    
    // cal string
    if (digiHeader) {
        digiHeader->setString("CaloString","Demo Cal String") ;
    }
    
    // end
	return StatusCode::SUCCESS ;
}
    
testCaloFhSetAlg::testCaloFhSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode testCaloFhSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode testCaloFhSetAlg::execute() {
	return StatusCode::SUCCESS ;
}
    
