
#include "GaudiKernel/Algorithm.h"

/** @class CaloFileHeadersSetAlg
 * @brief Prepare the calo specific attributes and set them in the relevant headers.
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/digiRootWriterAlg.cxx,v 1.24 2004/03/25 20:18:14 heather Exp $
 */

class CaloFileHeadersSetAlg : public Algorithm  {
        
public:
    
    CaloFileHeadersSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    /// 
    StatusCode initialize();
   
    /// 
    StatusCode execute();
    
    /// prepare the calo specific attributes and set them in the relevant headers
    StatusCode finalize();
    
} ;

