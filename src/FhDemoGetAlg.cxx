
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

#include <string>
#include "TOrdCollection.h"

#ifdef DEFECT_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

/** @class FhDemoGetAlg
 * @brief Takes and display few headers attributes
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/FhDemoGetAlg.cxx,v 1.4.672.1 2010/08/31 03:08:50 heather Exp $
 */

class FhDemoGetAlg : public Algorithm
{	
public:
    
    FhDemoGetAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
    
private:

    StatusCode finalize_common( MsgStream &,
      const std::string &, const FileHeader * );   
    
    IFhTool * m_headersTool ;
};

//static const AlgFactory<FhDemoGetAlg>  Factory;
//const IAlgFactory& FhDemoGetAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(FhDemoGetAlg);


FhDemoGetAlg::FhDemoGetAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode FhDemoGetAlg::initialize()
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

StatusCode FhDemoGetAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

StatusCode FhDemoGetAlg::finalize()
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

StatusCode FhDemoGetAlg::finalize_common(
    MsgStream & log, const std::string & prefix, const FileHeader * header ) {
    
    if (!header) return StatusCode::SUCCESS ;

    std::string word, word2 ;
    
    FhSystemEnv systemEnv ;
    systemEnv.init(header) ;
    TString cmtPath(systemEnv.getValue("CMTPATH")) ;
    TString rootIoShr(systemEnv.getValue("RootIoShr")) ;
    rootIoShr.Remove(0,cmtPath.Length()) ;
    log << MSG::INFO << prefix << " SystemEnv.RootIoShr: $CMTPATH" << rootIoShr << endreq ;

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
          <<"ChronoStatTable."<<tagName->String()
// DavidC: I do not display the time so to make the output stable for diff
//          <<": "<<chronoStatTable.getStringTime(*tagName)
          <<endreq ;
    }
    delete iter ;
    
    return StatusCode::SUCCESS ;
}
