
#include "src/testCaloFileHeadersSetAlg.h"
#include "src/FileHeadersTool.h"

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

static const AlgFactory<testCaloFileHeadersSetAlg> Factory ;
const IAlgFactory& testCaloFileHeadersSetAlgFactory = Factory ;

StatusCode testCaloFileHeadersSetAlg::finalize() {
    
    MsgStream log(msgSvc(),name()) ;
    StatusCode sc ;
    
    // retrieve the headers
    IFileHeadersTool * headersTool ;
    sc = toolSvc()->retrieveTool("FileHeadersTool",headersTool) ;
    if (sc.isFailure()) {
        log<<MSG::WARNING << "Failed to retrieve headers tool" << endreq ;
        return StatusCode::FAILURE ;
    }
    FileHeader * mcHeader = headersTool->mcHeader() ;
    FileHeader * digiHeader = headersTool->digiHeader() ;
    FileHeader * reconHeader = headersTool->reconHeader() ;
    
    // cal string
    if (digiHeader) {
        digiHeader->setString("caloString","Demo Cal String") ;
    }
    
    // end
	return StatusCode::SUCCESS ;
}
    
testCaloFileHeadersSetAlg::testCaloFileHeadersSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode testCaloFileHeadersSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode testCaloFileHeadersSetAlg::execute() {
	return StatusCode::SUCCESS ;
}
    
