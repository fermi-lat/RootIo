
#ifndef FhTool_H
#define FhTool_H

#include "commonRootData/FileHeader.h"
#include <GaudiKernel/AlgTool.h>
#include <GaudiKernel/StatusCode.h>
#include <TFile.h>
#include <iostream>

/*!

 @class IFhTool        
 @brief Interface for the tool managing the instances of FileHeader

 Interface for FhTool, the tool which manage the instances
 of FileHeader. Using a Gaudi tool for such task permits any piece
 of client code to access the shared headers, provided it can retrieve
 the tool.
 
 At any given time, the tool can keep simultaneously six header
 objects in memory : for each kind of file (mc, digi and recon),
 there can be one non-const header (useful for the write jobs) and
 one const header (useful for the read jobs). This should be enough
 for any of our use-cases.
 
 The tool do not systematically create the headers, nor write or
 read them the ROOT files. This should be asked by the classes which
 are driving the IO tasks, such as mcRootReaderAlg or reconRootWriterAlg.
 This is why, when any other client ask for a given header, it must
 check of the result is not 0 before using it. A 0 is not necessarily
 an error. For example, it is perfectly valid that there is no recon header
 in a job which is reading mc files and producing digi files.

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/
/** @class IFhTool
 * @brief Interface for FhTool.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/FhTool.h,v 1.1 2004/10/07 12:15:16 chamont Exp $
 */

static const InterfaceID IID_IFhTool("IFhTool",1,0) ;

class IFhTool : virtual public IAlgTool {

public:

    /// retreive interface id
    static const InterfaceID & interfaceID() { return IID_IFhTool ; }
    
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

    /// create a writable gcr header
    virtual StatusCode newGcrHeader() =0 ;

    /// access to the current writable gcr header
    virtual FileHeader * gcrHeader() =0 ;
    
    /// write the recon header
    virtual StatusCode writeGcrHeader( TFile * ) =0 ;
    
    /// if the file is new, extract its header
    virtual StatusCode readConstGcrHeader( TFile * ) =0 ;

    /// access to the current read-only gcr header
    virtual const FileHeader * constGcrHeader() =0 ;

    
} ; 
 
/** @class FhTool
 * @brief Tool which manage the files headers.
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/FhTool.h,v 1.1 2004/10/07 12:15:16 chamont Exp $
 */

class FhTool : public AlgTool, virtual public IFhTool {

public:

    FhTool( const std::string & type,
        const std::string & name, const IInterface * parent ) ;
    virtual ~FhTool() ;
    
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
    
    StatusCode newGcrHeader() ;
    FileHeader * gcrHeader() ;
    StatusCode writeGcrHeader( TFile * ) ;
    StatusCode readConstGcrHeader( TFile * ) ;
    const FileHeader * constGcrHeader() ;

    
private:

    StatusCode writeHeader( TFile *, FileHeader * ) ;
    StatusCode readHeader( TFile *, const FileHeader * & ) ;

    FileHeader * m_mcHeader ;
    FileHeader * m_digiHeader ;
    FileHeader * m_reconHeader ;
    FileHeader * m_gcrHeader ;

    const FileHeader * m_constMcHeader ;
    const FileHeader * m_constDigiHeader ;
    const FileHeader * m_constReconHeader ;
    const FileHeader * m_constGcrHeader ;

} ; 
 
#endif





