
#ifndef FhDemoCaloSetAlg_H
#define FhDemoCaloSetAlg_H

#include "GaudiKernel/Algorithm.h"

/*!

 @class FhDemoCaloSetAlg        
 @brief Prepare the calo specific attributes and set them in the relevant file headers

 This Gaudi algorithm is expected to prepare and set the file header data which is
 specific to calorimetry. It is a kind of demo class, to be reused as is and/or
 eventually transformed by calo developers. It should be declared in the job options
 top algorithms.
 
 As any algorithm which is preparing file headers data, it is using the tool
 FhTool so to get access to the current writable headers. Worth to note,
 the task is done during the finalize() step.

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

*/

class FhDemoCaloSetAlg : public Algorithm  {
        
public:
    
    FhDemoCaloSetAlg( const std::string& name, ISvcLocator* pSvcLocator ) ;
    
    /// 
    StatusCode initialize();
   
    /// 
    StatusCode execute();
    
    /// prepare the calo specific attributes and set them in the relevant headers
    StatusCode finalize();
    
} ;

#endif

