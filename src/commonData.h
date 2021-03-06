#ifndef ROOTIO_COMMONDATA_H
#define ROOTIO_COMMONDATA_H 1


//#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTrajectory.h"

#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/Digi/CalDigi.h"
#include "Event/Digi/TkrDigi.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/TkrRecon/TkrTree.h"

#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"

#include "TFile.h"
#include "TTree.h"
#include "TRef.h"
#include "mcRootData/McEvent.h"
#include "digiRootData/DigiEvent.h"
#include "reconRootData/ReconEvent.h"

#include "facilities/Util.h"

#include <map>

/** @class commonData
* @brief Provides global RootIo access to the ROOT data pointers.  This is 
* necessary when writing/reading relations which cross the boundaries of
* Monte Carlo, Digitization, and Reconstruction data.
*
* @author Heather Kelly
* $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/src/commonData.h,v 1.11 2013/02/05 15:52:28 usher Exp $
*/

class commonData 
{	
public:
    
    commonData() {  }; 
    ~commonData() {  };
    
    void clear();
    
    /// Create a set of maps between MC data in TDS and the TRefs in the ROOT file
    static std::map<const Event::McParticle*, TRef> m_mcPartMap;
    static std::map<const Event::McIntegratingHit*, TRef> m_mcIntHitMap;
    static std::map<const Event::McPositionHit*, TRef> m_mcPosHitMap;
    static std::map<const Event::McTrajectory*, TRef> m_mcTrajectoryMap;
    static std::map<const Event::McTrajectoryPoint*, TRef> m_mcTrajectoryPointMap;
    
    /// Create a set of maps between Digi data in the TDS and the TRefs in the ROOT file
    static std::map<const Event::TkrDigi*, TRef> m_tkrDigiMap;
    static std::map<const Event::CalDigi*, TRef> m_calDigiMap;
    
    /// Create a set of maps between Recon data in TDS and TRefs in ROOT
    static std::map<const Event::TkrCluster*,       TRef> m_tkrClusterMap;
    static std::map<const Event::TkrTrack*,         TRef> m_tkrTrackMap;
    static std::map<const Event::TkrVertex*,        TRef> m_tkrVertexMap;
    static std::map<const Event::TkrVecPoint*,      TRef> m_tkrVecPointMap;
    static std::map<const Event::TkrVecPointsLink*, TRef> m_tkrVecPointsLinkMap;
    static std::map<const Event::TkrVecNode*,       TRef> m_tkrVecNodeMap;
    static std::map<const Event::TkrTree*,          TRef> m_tkrTreeMap;
    static std::map<const Event::TkrFilterParams*,  TRef> m_tkrFilterParamsMap;

    static std::map<const Event::CalCluster*,       TRef> m_calClusterMap;
    static std::map<const Event::CalXtalRecData*,   TRef> m_calXtalRecDataMap;
    static std::map<const Event::CalEventEnergy*,   TRef> m_calEventEnergyMap;

    /// Create a set of maps between ROOT MC objects and the TDS MC data
    static std::map<const TObject*, const Event::McParticle*>        m_rootMcPartMap;
    static std::map<const TObject*, const Event::McIntegratingHit*>  m_rootMcIntHitMap;
    static std::map<const TObject*, const Event::McPositionHit*>     m_rootMcPosHitMap;
    static std::map<const TObject*, const Event::McTrajectory*>      m_rootMcTrajectoryMap;
    static std::map<const TObject*, const Event::McTrajectoryPoint*> m_rootMcTrajectoryPointMap;
    
    /// Create a set of maps between ROOT Digi objects and TDS Digi data
    static std::map<const TObject*, const Event::TkrDigi*>          m_rootTkrDigiMap;
    static std::map<const TObject*, const Event::CalDigi*>          m_rootCalDigiMap;

    /// Create a set of maps between ROOT Recon objects and TDS Recon data
    static std::map<const TObject*, const Event::TkrCluster*>       m_rootTkrClusterMap;
    static std::map<const TObject*, const Event::TkrTrack*>         m_rootTkrTrackMap;
    static std::map<const TObject*, const Event::TkrVertex*>        m_rootTkrVertexMap;
    static std::map<const TObject*, const Event::TkrVecPoint*>      m_rootTkrVecPointMap;
    static std::map<const TObject*, const Event::TkrVecPointsLink*> m_rootTkrVecPointsLinkMap;
    static std::map<const TObject*, const Event::TkrVecNode*>       m_rootTkrVecNodeMap;
    static std::map<const TObject*, const Event::TkrTree*>          m_rootTkrTreeMap;
    static std::map<const TObject*, const Event::TkrFilterParams*>  m_rootTkrFilterParamsMap;

    static std::map<const TObject*, const Event::CalCluster*>       m_rootCalClusterMap;
    static std::map<const TObject*, const Event::CalXtalRecData*>   m_rootCalXtalRecDataMap;
    static std::map<const TObject*, const Event::CalEventEnergy*>   m_rootCalEventEnergyMap;

    /// Provide access to the ROOT event pointers
    static McEvent *m_mcEvt;
    static DigiEvent *m_digiEvt;
    static ReconEvent *m_reconEvt;
    
    
};




#endif
