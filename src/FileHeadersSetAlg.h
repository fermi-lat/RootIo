
#include "GaudiKernel/Algorithm.h"

/** @class FileHeadersSetAlg
 * @brief Prepare the common attributes and set them in available headers.
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootWriterAlg.cxx,v 1.24 2004/03/25 20:18:14 heather Exp $
 */

class FileHeadersSetAlg : public Algorithm {	

public:

    FileHeadersSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    ///
    StatusCode initialize();
   
    ///
    StatusCode execute();
    
    /// prepare the common attributes and set them in available headers
    StatusCode finalize();
    
} ;

