
#include "RootIo/FhTool.h"
#include "RootIo/FhSystemEnv.h"
#include "RootIo/FhCmtConfig.h"
#include "RootIo/FhJobOptions.h"
#include "RootIo/FhChronoStatTable.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include <sstream>
#include <string>
#include "TOrdCollection.h"

/** @class testFhGetAlg
 * @brief Takes and display few headers attributes
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/testFhAlg.cxx,v 1.2 2004/09/17 12:38:43 chamont Exp $
 */

class testFhGetAlg : public Algorithm
{	
public:
    
    testFhGetAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
    
private:

    StatusCode finalize_common( MsgStream &,
      const std::string &, const FileHeader * );   
    
    IFhTool * m_headersTool ;
};

static const AlgFactory<testFhGetAlg>  Factory;
const IAlgFactory& testFhGetAlgFactory = Factory;


testFhGetAlg::testFhGetAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode testFhGetAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();
    
    StatusCode headersSc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }

    return sc;
    
}

StatusCode testFhGetAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

StatusCode testFhGetAlg::finalize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // common part of each header
    if (finalize_common(log,"Mc",m_headersTool->mcHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Digi",m_headersTool->digiHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Recon",m_headersTool->reconHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Mc Const",m_headersTool->constMcHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Digi Const",m_headersTool->constDigiHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Recon Const",m_headersTool->constReconHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }

    // digi specific test
    const FileHeader * digiHeader = m_headersTool->digiHeader() ;
    if (digiHeader) {
        log << MSG::INFO << "Digi CaloString: " << digiHeader->getString("CaloString") << endreq ;
    }
    digiHeader = m_headersTool->constDigiHeader() ;
    if (digiHeader) {
        log << MSG::INFO << "Digi Const CaloString: " << digiHeader->getString("CaloString") << endreq ;
    }

    return sc;
}

StatusCode testFhGetAlg::finalize_common(
    MsgStream & log, const std::string & prefix, const FileHeader * header ) {
    
    if (!header) return StatusCode::SUCCESS ;

    std::string word, word2 ;
    
    FhSystemEnv systemEnv ;
    systemEnv.init(header) ;
    log << MSG::INFO << prefix << " SystemEnv.RootIoShr: " << systemEnv.getValue("RootIoShr") << endreq ;

    FhCmtConfig cmtConfig ;
    cmtConfig.init(header) ;
    log << MSG::INFO << prefix << " CmtUses.Event: " << cmtConfig.getPackageSelection("Event") << endreq ;
    log << MSG::INFO << prefix << " CmtMacros.package: " << cmtConfig.getMacroValue("package") << endreq ;
    log << MSG::INFO << prefix << " CmtMacros.ROOTIOVERSION: " << cmtConfig.getMacroValue("ROOTIOVERSION") << endreq ;

    // all job options
    FhJobOptions jobOptions ;
    jobOptions.init(header) ;
    TOrdCollection jobOptionNames ;
    jobOptions.getNames(".*",jobOptionNames) ;
    jobOptionNames.Sort() ;
    TIterator * iter = jobOptionNames.MakeIterator() ;
    TObjString * jobOptionName ;
    while (( jobOptionName = (TObjString *)(iter->Next()) )) {
        log<<MSG::INFO<<prefix<<" "
          <<"JobOptions."<<jobOptionName->String()<<": "
          <<jobOptions.getValue(*jobOptionName)
          <<endreq ;
    }
    delete iter ;
    
    // all chrono stat table
    FhChronoStatTable chronoStatTable ;
    chronoStatTable.init(header) ;
    TOrdCollection tagNames ;
    chronoStatTable.getTagNames(".*",tagNames) ;
    tagNames.Sort() ;
    iter = tagNames.MakeIterator() ;
    TObjString * tagName ;
    while (( tagName = (TObjString *)(iter->Next()) )) {
        log<<MSG::INFO<<prefix<<" "
          <<"ChronoStatTable."<<tagName->String()<<": "
          <<chronoStatTable.getStringTime(*tagName)
          <<endreq ;
    }
    delete iter ;
    
    return StatusCode::SUCCESS ;
}
