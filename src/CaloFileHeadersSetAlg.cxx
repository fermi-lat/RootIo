
#include "src/CaloFileHeadersSetAlg.h"
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

static const AlgFactory<CaloFileHeadersSetAlg> Factory ;
const IAlgFactory& CaloFileHeadersSetAlgFactory = Factory ;

StatusCode CaloFileHeadersSetAlg::finalize() {
    
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
    
CaloFileHeadersSetAlg::CaloFileHeadersSetAlg
( const std::string& name, ISvcLocator* pSvcLocator)
:  Algorithm(name, pSvcLocator) {
}

StatusCode CaloFileHeadersSetAlg::initialize() {
	return StatusCode::SUCCESS ;
}

StatusCode CaloFileHeadersSetAlg::execute() {
	return StatusCode::SUCCESS ;
}
    
