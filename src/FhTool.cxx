
#include "RootIo/FhTool.h"
#include <GaudiKernel/ToolFactory.h>
#include <GaudiKernel/MsgStream.h>

//static ToolFactory<FhTool> s_factory ;
//const IToolFactory & FhToolFactory = s_factory ;
DECLARE_TOOL_FACTORY(FhTool);

FhTool::FhTool( const std::string & type,
    const std::string & name, const IInterface * parent )
    : AlgTool(type,name,parent) {
    	
    declareInterface<IFhTool>(this) ;
    
    m_mcHeader = 0 ;
    m_digiHeader = 0 ;
    m_reconHeader = 0 ;
    m_gcrHeader = 0 ;
    m_meritHeader = 0 ;
   
    m_constMcHeader = 0 ;
    m_constDigiHeader = 0 ;
    m_constReconHeader = 0 ;
    m_constGcrHeader = 0 ;
    m_constMeritHeader = 0 ;
   
}

FhTool::~FhTool() {
    	
    delete m_mcHeader ;
    delete m_digiHeader ;
    delete m_reconHeader ;
    delete m_gcrHeader ;
    delete m_meritHeader ;
   
    delete m_constMcHeader ;
    delete m_constDigiHeader ;
    delete m_constReconHeader ;
    delete m_constGcrHeader ;
    delete m_constMeritHeader ;
   
}

FileHeader * FhTool::mcHeader() {
	return m_mcHeader ;
}

FileHeader * FhTool::digiHeader() {
	return m_digiHeader ;
}

FileHeader * FhTool::reconHeader() {
	return m_reconHeader ;
}

FileHeader * FhTool::gcrHeader() {
	return m_gcrHeader ;
}

FileHeader * FhTool::meritHeader() {
	return m_meritHeader ;
}

const FileHeader * FhTool::constMcHeader() {
	return m_constMcHeader ;
}

const FileHeader * FhTool::constDigiHeader() {
	return m_constDigiHeader ;
}

const FileHeader * FhTool::constReconHeader() {
	return m_constReconHeader ;
}

const FileHeader * FhTool::constGcrHeader() {
	return m_constGcrHeader ;
}

const FileHeader * FhTool::constMeritHeader() {
	return m_constMeritHeader ;
}

StatusCode FhTool::newMcHeader() {
    if (m_mcHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::WARNING<<"Mc FileHeader already existing"<<endreq ;
    }
    else {
        m_mcHeader = new FileHeader ;
    }
	return StatusCode::SUCCESS ;
}

StatusCode FhTool::newDigiHeader() {
    if (m_digiHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::WARNING<<"Digi FileHeader already existing"<<endreq ;
    }
    else {
        m_digiHeader = new FileHeader ;
    }
	return StatusCode::SUCCESS ;
}

StatusCode FhTool::newReconHeader() {
    if (m_reconHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::WARNING<<"Recon FileHeader already existing"<<endreq ;
    }
    else {
        m_reconHeader = new FileHeader ;
    }
	return StatusCode::SUCCESS ;
}

StatusCode FhTool::newGcrHeader() {
    if (m_gcrHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::WARNING<<"Gcr FileHeader already existing"<<endreq ;
    }
    else {
        m_gcrHeader = new FileHeader ;
    }
	return StatusCode::SUCCESS ;
}

StatusCode FhTool::newMeritHeader() {
    if (m_meritHeader) {
        MsgStream log(msgSvc(),name()) ;
	    log<<MSG::WARNING<<"Merit FileHeader already existing"<<endreq ;
    }
    else {
        m_meritHeader = new FileHeader ;
    }
	return StatusCode::SUCCESS ;
}

StatusCode FhTool::writeMcHeader( TFile * file ) {
	return writeHeader(file,m_mcHeader) ;
}

StatusCode FhTool::writeDigiHeader( TFile * file ) {
	return writeHeader(file,m_digiHeader) ;
}

StatusCode FhTool::writeReconHeader( TFile * file ) {
	return writeHeader(file,m_reconHeader) ;
}

StatusCode FhTool::writeGcrHeader( TFile * file ) {
	return writeHeader(file,m_gcrHeader) ;
}

StatusCode FhTool::writeMeritHeader( TFile * file ) {
	return writeHeader(file,m_meritHeader) ;
}

StatusCode FhTool::readConstMcHeader( TFile * file ) {
	return readHeader(file,m_constMcHeader) ;
}

StatusCode FhTool::readConstDigiHeader( TFile * file ) {
	return readHeader(file,m_constDigiHeader) ;
}

StatusCode FhTool::readConstReconHeader( TFile * file ) {
	return readHeader(file,m_constReconHeader) ;
}

StatusCode FhTool::readConstGcrHeader( TFile * file ) {
	return readHeader(file,m_constGcrHeader) ;
}

StatusCode FhTool::readConstMeritHeader( TFile * file ) {
	return readHeader(file,m_constMeritHeader) ;
}

StatusCode FhTool::writeHeader( TFile * file,
  FileHeader * header ) {
      
    MsgStream log(msgSvc(),name()) ;
	
    static TFile * oldFile = 0 ;
	if (!header) {
	    log<<MSG::WARNING<<"No file header found"<<endreq ;
	}
    else if ( file != oldFile ) {
    	
        log << MSG::DEBUG << "BEGIN writeHeader" << endreq ;

        oldFile = file ;

	    //header->Print() ;
	    
	    TDirectory * saveDir = gDirectory ;  
	    file->cd() ;
	    header->Write("header") ;
	    saveDir->cd() ;
	    
        log << MSG::DEBUG << "END writeHeader" << endreq ;
    }
    
    return StatusCode::SUCCESS ;
    
}

StatusCode FhTool::readHeader( TFile * file,
  const FileHeader * & header ) {

    static TFile * oldFile = 0 ;
    if ( file != oldFile ) {
        oldFile = file ;
	    // just get the header from the file
	    // and call TObject::Print()
	    TObject * object = file->Get("header") ;
	    if (!object) {
            MsgStream log(msgSvc(),name()) ;
	        log<<MSG::WARNING<<"No file header found"<<endreq ;
	    }
	    else { 
            delete header ;
	    	header = (const FileHeader *) object ;
	    	//header->Print() ;
        }
    }
	return StatusCode::SUCCESS ;

}

