
#ifndef FhSetAlg_H
#define FhSetAlg_H


#include "GaudiKernel/Algorithm.h"

/*!

 @class FhSetAlg        
 @brief Prepare the common attributes and set them in all available headers

 This Gaudi algorithm is expected to prepare and set the file header data which is
 common to all files. Its implementation implicitly define what is the header
 data common to all ROOT files. It should be declared in the job options top algorithms,
 each time that a new ROOT file should be produced by the job.
 
 As any algorithm which is preparing file headers data, it is using the tool
 FhTool so to get access to the current writable headers. On the other
 contrary of other algorithms, it is setting the same common data in all
 the available writable headers (mc, digi and recon).
 Worth to note, the task is done during the finalize() step.

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/
/** @class FhSetAlg
 * @brief Prepare the common attributes and set them in available headers.
 *
 * @author David Chamont
 * $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/FhSetAlg.h,v 1.3 2004/09/22 05:09:39 heather Exp $
 */

class FhSetAlg : public Algorithm {	

public:

    FhSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    ///
    StatusCode initialize();
   
    ///
    StatusCode execute();
    
    /// prepare the common attributes and set them in available headers
    StatusCode finalize();

private:
	/// provide a local popen method to handle both linux and windows
	FILE* popen(const char* command, const char *mode);
	/// provide a local pclose method to handle both linux and windows
	int pclose(FILE* fp);
    
} ;

#endif

