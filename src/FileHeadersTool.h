
#ifndef FileHeadersTool_H
#define FileHeadersTool_H

#include "commonRootData/FileHeader.h"
#include <GaudiKernel/AlgTool.h>
#include <GaudiKernel/StatusCode.h>
#include <TFile.h>
#include <iostream>

/** @class IFileHeadersTool
 * @brief Interface for FileHeadersTool.
 *
 * $Header: $
 */

static const InterfaceID IID_IFileHeadersTool("IFileHeadersTool",1,0) ;

class IFileHeadersTool : virtual public IAlgTool {

public:

    /// retreive interface id
    static const InterfaceID & interfaceID() { return IID_IFileHeadersTool ; }
    
    /// create a writable mc header
    virtual StatusCode newMcHeader() =0 ;
    
    /// access to the current writable mc file header
    virtual FileHeader * mcHeader() =0 ;
    
    /// write the mc header
    virtual StatusCode writeMcHeader( TFile * ) =0 ;
    
    /// if the file is new, extract its header
    virtual StatusCode readConstMcHeader( TFile * ) =0 ;

    /// access to the current read-only mc file header
    virtual const FileHeader * constMcHeader() =0 ;
    
    /// create a writable digi header
    virtual StatusCode newDigiHeader() =0 ;

    /// access to the current writable digi header
    virtual FileHeader * digiHeader() =0 ;
    
    /// write a digi header
    virtual StatusCode writeDigiHeader( TFile * ) =0 ;
    
    /// if the file is new, extract its header
    virtual StatusCode readConstDigiHeader( TFile * ) =0 ;

    /// access to the current read-only digi header
    virtual const FileHeader * constDigiHeader() =0 ;
    
    /// create a writable recon header
    virtual StatusCode newReconHeader() =0 ;

    /// access to the current writable recon header
    virtual FileHeader * reconHeader() =0 ;
    
    /// write the recon header
    virtual StatusCode writeReconHeader( TFile * ) =0 ;
    
    /// if the file is new, extract its header
    virtual StatusCode readConstReconHeader( TFile * ) =0 ;

    /// access to the current read-only recon header
    virtual const FileHeader * constReconHeader() =0 ;
    
} ; 
 
/** @class FileHeadersTool
 * @brief Tool which manage the files headers.
 *
 * $Header: $
 */

class FileHeadersTool : public AlgTool, virtual public IFileHeadersTool {

public:

    FileHeadersTool( const std::string & type,
        const std::string & name, const IInterface * parent ) ;
    virtual ~FileHeadersTool() ;
    
    StatusCode newMcHeader() ;
    FileHeader * mcHeader() ;
    StatusCode writeMcHeader( TFile * ) ;
    StatusCode readConstMcHeader( TFile * ) ;
    const FileHeader * constMcHeader() ;
    
    StatusCode newDigiHeader() ;
    FileHeader * digiHeader() ;
    StatusCode writeDigiHeader( TFile * ) ;
    StatusCode readConstDigiHeader( TFile * ) ;
    const FileHeader * constDigiHeader() ;
    
    StatusCode newReconHeader() ;
    FileHeader * reconHeader() ;
    StatusCode writeReconHeader( TFile * ) ;
    StatusCode readConstReconHeader( TFile * ) ;
    const FileHeader * constReconHeader() ;
    
private:

    StatusCode writeHeader( TFile *, FileHeader * ) ;
    StatusCode readHeader( TFile *, const FileHeader * & ) ;

    FileHeader * m_mcHeader ;
    FileHeader * m_digiHeader ;
    FileHeader * m_reconHeader ;

    const FileHeader * m_constMcHeader ;
    const FileHeader * m_constDigiHeader ;
    const FileHeader * m_constReconHeader ;

} ; 
 
#endif





