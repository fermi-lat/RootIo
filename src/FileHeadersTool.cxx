
#include "src/FileHeadersTool.h"
#include <GaudiKernel/ToolFactory.h>
#include <GaudiKernel/MsgStream.h>

static ToolFactory<FileHeadersTool> s_factory ;
const IToolFactory & FileHeadersToolFactory = s_factory ;

FileHeadersTool::FileHeadersTool( const std::string & type,
    const std::string & name, const IInterface * parent )
    : AlgTool(type,name,parent) {
    	
    declareInterface<IFileHeadersTool>(this) ;
    
    m_mcHeader = 0 ;
    m_digiHeader = 0 ;
    m_reconHeader = 0 ;
   
    m_constMcHeader = 0 ;
    m_constDigiHeader = 0 ;
    m_constReconHeader = 0 ;
   
}

FileHeadersTool::~FileHeadersTool() {
    	
    delete m_mcHeader ;
    delete m_digiHeader ;
    delete m_reconHeader ;
   
    delete m_constMcHeader ;
    delete m_constDigiHeader ;
    delete m_constReconHeader ;
   
}

FileHeader * FileHeadersTool::mcHeader() {
	return m_mcHeader ;
}

FileHeader * FileHeadersTool::digiHeader() {
	return m_digiHeader ;
}

FileHeader * FileHeadersTool::reconHeader() {
	return m_reconHeader ;
}

const FileHeader * FileHeadersTool::constMcHeader() {
	return m_constMcHeader ;
}

const FileHeader * FileHeadersTool::constDigiHeader() {
	return m_constDigiHeader ;
}

const FileHeader * FileHeadersTool::constReconHeader() {
	return m_constReconHeader ;
}

StatusCode FileHeadersTool::newMcHeader() {
    if (m_mcHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::ERROR<<"Mc FileHeader already existing"<<endreq ;
        return StatusCode::FAILURE ;
    }
    m_mcHeader = new FileHeader ;
	return StatusCode::SUCCESS ;
}

StatusCode FileHeadersTool::newDigiHeader() {
    if (m_digiHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::ERROR<<"Digi FileHeader already existing"<<endreq ;
        return StatusCode::FAILURE ;
    }
    m_digiHeader = new FileHeader ;
	return StatusCode::SUCCESS ;
}

StatusCode FileHeadersTool::newReconHeader() {
    if (m_reconHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::ERROR<<"Recon FileHeader already existing"<<endreq ;
        return StatusCode::FAILURE ;
    }
    m_reconHeader = new FileHeader ;
	return StatusCode::SUCCESS ;
}

StatusCode FileHeadersTool::writeMcHeader( TFile * file ) {
	return writeHeader(file,m_mcHeader) ;
}

StatusCode FileHeadersTool::writeDigiHeader( TFile * file ) {
	return writeHeader(file,m_digiHeader) ;
}

StatusCode FileHeadersTool::writeReconHeader( TFile * file ) {
	return writeHeader(file,m_reconHeader) ;
}

StatusCode FileHeadersTool::readConstMcHeader( TFile * file ) {
	return readHeader(file,m_constMcHeader) ;
}

StatusCode FileHeadersTool::readConstDigiHeader( TFile * file ) {
	return readHeader(file,m_constDigiHeader) ;
}

StatusCode FileHeadersTool::readConstReconHeader( TFile * file ) {
	return readHeader(file,m_constReconHeader) ;
}

StatusCode FileHeadersTool::writeHeader( TFile * file,
  FileHeader * header ) {
	
    static TFile * oldFile = 0 ;
    if ( file != oldFile ) {
    	
        oldFile = file ;

	    //header->Print() ;
	    
	    TDirectory * saveDir = gDirectory ;  
	    file->cd() ;
	    header->Write("header") ;
	    saveDir->cd() ;
	    
    }
    
    return StatusCode::SUCCESS ;
    
}

StatusCode FileHeadersTool::readHeader( TFile * file,
  const FileHeader * & header ) {

    static TFile * oldFile = 0 ;
    if ( file != oldFile ) {
        oldFile = file ;
	    // just get the header from the file
	    // and call TObject::Print()
	    TObject * object = file->Get("header") ;
	    if (!object) {
            MsgStream log(msgSvc(),name()) ;
	        log<<MSG::WARNING<<"Failed to get the file header"<<endreq ;
	    	return StatusCode::FAILURE ;
	    }
	    else { 
            delete header ;
	    	header = (const FileHeader *) object ;
	    	//header->Print() ;
        }
    }
	return StatusCode::SUCCESS ;

}

