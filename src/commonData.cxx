#include "commonData.h"


std::map<const Event::McParticle*, TRef>                  commonData::m_mcPartMap;
std::map<const Event::McIntegratingHit*, TRef>            commonData::m_mcIntHitMap;
std::map<const Event::McPositionHit*, TRef>               commonData::m_mcPosHitMap;
std::map<const Event::McTrajectory*, TRef>                commonData::m_mcTrajectoryMap;
std::map<const Event::McTrajectoryPoint*, TRef>           commonData::m_mcTrajectoryPointMap;
std::map<const Event::TkrDigi*, TRef>                     commonData::m_tkrDigiMap;
std::map<const Event::CalDigi*, TRef>                     commonData::m_calDigiMap;
std::map<const Event::TkrCluster*, TRef>                  commonData::m_tkrClusterMap;
std::map<const Event::TkrTrack*, TRef>                    commonData::m_tkrTrackMap;
std::map<const Event::TkrVertex*, TRef>                   commonData::m_tkrVertexMap;
std::map<const Event::CalCluster*, TRef>                  commonData::m_calClusterMap;
std::map<const Event::CalXtalRecData*, TRef>              commonData::m_calXtalRecDataMap;

std::map<const TObject*, const Event::McParticle*>        commonData::m_rootMcPartMap;
std::map<const TObject*, const Event::McIntegratingHit*>  commonData::m_rootMcIntHitMap;
std::map<const TObject*, const Event::McPositionHit*>     commonData::m_rootMcPosHitMap;
std::map<const TObject*, const Event::McTrajectory*>      commonData::m_rootMcTrajectoryMap;
std::map<const TObject*, const Event::McTrajectoryPoint*> commonData::m_rootMcTrajectoryPointMap;
std::map<const TObject*, const Event::TkrDigi*>           commonData::m_rootTkrDigiMap;
std::map<const TObject*, const Event::CalDigi*>           commonData::m_rootCalDigiMap;
std::map<const TObject*, const Event::TkrCluster*>        commonData::m_rootTkrClusterMap;
std::map<const TObject*, const Event::TkrTrack*>          commonData::m_rootTkrTrackMap;
std::map<const TObject*, const Event::TkrVertex*>         commonData::m_rootTkrVertexMap;
std::map<const TObject*, const Event::CalCluster*>        commonData::m_rootCalClusterMap;
std::map<const TObject*, const Event::CalXtalRecData*>    commonData::m_rootCalXtalRecDataMap;

McEvent* commonData::m_mcEvt;
DigiEvent* commonData::m_digiEvt;
ReconEvent* commonData::m_reconEvt;


void commonData::clear() { 
    m_mcPartMap.clear(); 
    m_mcIntHitMap.clear();
    m_mcPosHitMap.clear();
    m_mcTrajectoryMap.clear();
    m_mcTrajectoryPointMap.clear();
    m_tkrDigiMap.clear();
    m_calDigiMap.clear(); 
    m_tkrClusterMap.clear();
    m_tkrTrackMap.clear();
    m_tkrVertexMap.clear();
    m_calClusterMap.clear();
    m_calXtalRecDataMap.clear();

    m_rootMcPartMap.clear(); 
    m_rootMcIntHitMap.clear();
    m_rootMcPosHitMap.clear();
    m_rootMcTrajectoryMap.clear();

    m_rootTkrDigiMap.clear();
    m_rootCalDigiMap.clear();
    
    m_rootTkrClusterMap.clear();
    m_rootTkrTrackMap.clear();
    m_rootTkrVertexMap.clear();
    m_rootCalClusterMap.clear();
    m_rootCalXtalRecDataMap.clear();
}
