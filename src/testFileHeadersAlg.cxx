
#include "src/FileHeadersTool.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include <sstream>
#include <string>
#include <map>

/** @class testFileHeadersSetAlg
 * @brief Takes and display few headers attributes
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/test/testReadAlg.cxx,v 1.9 2003/08/21 18:29:49 heather Exp $
 */

class testFileHeadersSetAlg : public Algorithm
{	
public:
    
    testFileHeadersSetAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    StatusCode initialize();
   
    StatusCode execute();
    
    StatusCode finalize();   
    
private:

    StatusCode finalize_common( MsgStream &,
      const std::string &, const FileHeader * );   
    
    IFileHeadersTool * m_headersTool ;
};

static const AlgFactory<testFileHeadersSetAlg>  Factory;
const IAlgFactory& testFileHeadersSetAlgFactory = Factory;


testFileHeadersSetAlg::testFileHeadersSetAlg(const std::string& name, ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
}

StatusCode testFileHeadersSetAlg::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    setProperties();
    
    StatusCode headersSc = toolSvc()->retrieveTool("FileHeadersTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }

    return sc;
    
}

StatusCode testFileHeadersSetAlg::execute()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

StatusCode testFileHeadersSetAlg::finalize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // common part of each header
    if (finalize_common(log,"Mc FileHeader",m_headersTool->mcHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Digi FileHeader",m_headersTool->digiHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Recon FileHeader",m_headersTool->reconHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Mc Const FileHeader",m_headersTool->constMcHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Digi Const FileHeader",m_headersTool->constDigiHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }
    if (finalize_common(log,"Recon Const FileHeader",m_headersTool->constReconHeader())==StatusCode::FAILURE)
     { sc = StatusCode::FAILURE ; }

    // digi specific test
    const FileHeader * digiHeader = m_headersTool->digiHeader() ;
    if (digiHeader) {
        log << MSG::INFO << "Digi FileHeader CaloString: " << digiHeader->getString("caloString") << endreq ;
    }
    digiHeader = m_headersTool->constDigiHeader() ;
    if (digiHeader) {
        log << MSG::INFO << "Digi Const FileHeader CaloString: " << digiHeader->getString("caloString") << endreq ;
    }

    return sc;
}

StatusCode testFileHeadersSetAlg::finalize_common(
    MsgStream & log, const std::string & prefix, const FileHeader * header ) {

    if (!header) return StatusCode::SUCCESS ;
    
//    better stored event by event  
//    log << MSG::INFO << prefix << " RunId: " << header->getString("runId") << endreq ;

    std::istringstream cmtUses(header->getString("cmtUses").Data()) ;
    std::string word, word2 ;
    while ((cmtUses>>word)&&(word!="Selection")) ;
    while (cmtUses>>word) {
    	if (word=="Event") {
            cmtUses>>word2 ;
            log << MSG::INFO << prefix << " Cmt Use: Event " << word2 << endreq ;
    	}
    }

    std::istringstream cmtMacros(header->getString("cmtMacros").Data()) ;
    while (cmtMacros>>word) {
    	if (word.find("package=")==0) {
            log << MSG::INFO << prefix << " Cmt Macro: " << word << endreq ;
    	}
    	else if (word.find("version=")==0) {
            log << MSG::INFO << prefix << " Cmt Macro: " << word << endreq ;
    	}
    }

    std::istringstream jobOptions(header->getString("jobOptions").Data()) ;
    while (jobOptions>>word>>word2) {
    	log << MSG::INFO << prefix << " Job Option: " << word << " " << word2 << endreq ;
    }
    
    std::istringstream chronoStatTable(header->getString("chronoStatTable").Data()) ;
    std::string words[7], value, unit, name ;
    int iword = 0 ;
    while (chronoStatTable>>words[iword]) {
        if (words[iword]=="Tot=") {
            name = words[(iword-5+7)%7] ;
            chronoStatTable>>value>>unit ;
            log << MSG::INFO << prefix << " Chrono Stat Table: "
              << name << " " << value << " " <<unit  << endreq ;
        }
        iword = (iword+1)%7 ;
    }
    
    return StatusCode::SUCCESS ;
}
