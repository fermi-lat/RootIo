#ifndef ROOTIO_COMMONDATA_H
#define ROOTIO_COMMONDATA_H 1


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"

#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"


#include "TFile.h"
#include "TTree.h"
#include "TRef.h"
#include "mcRootData/McEvent.h"
#include "digiRootData/DigiEvent.h"

#include "facilities/Util.h"

#include <map>

/** @class commonData
* @brief Provides global RootIo access to the ROOT data pointers.  This is 
* necessary when writing/reading relations which cross the boundaries of
* Monte Carlo, Digitization, and Reconstruction data.
*
* @author Heather Kelly
* $Header$
*/

class commonData 
{	
public:
    
    commonData() {  }; 
    ~commonData() {  };
    
    void clear();
    
    static std::map<const Event::McParticle*, TRef> m_mcPartMap;
    static std::map<const Event::McIntegratingHit*, TRef> m_mcIntHitMap;
    static std::map<const Event::McPositionHit*, TRef> m_mcPosHitMap;
    
    static std::map<const Event::TkrDigi*, TRef> m_tkrDigiMap;
    static std::map<const Event::CalDigi*, TRef> m_calDigiMap;
    
    static McEvent *m_mcEvt;
    static DigiEvent *m_digiEvt;
    
    
};




#endif