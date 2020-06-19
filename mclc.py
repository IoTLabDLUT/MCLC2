#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import matlab.engine
import math
from typing import List, Optional, Tuple, Any

deltaT = 1
TIME_HEADWAY = 3.0
VEHICLE_LENGTH = 5.0
inf: float = float('Inf')
NUM_ENGINE = 3

# we need to import python modules from the $SUMO_HOME/tools directory
try:
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', "tools"))  # tutorial in tests
    sys.path.append(os.path.join(os.environ.get("SUMO_HOME", os.path.join(os.path.dirname(__file__), "..", "..",
                                                                          "..")), "tools"))  # tutorial in docs
    from sumolib import checkBinary
    import traci
    from traci import vehicle, lane, simulation
except ImportError:
    sys.exit("please declare environment variable 'SUMO_HOME' as the root directory of your sumo installation"
             "(it should contain folders 'bin', 'tools' and 'docs')")


def get_gcv(subject_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    lane_id = vehicle.getLaneID(subject_id)
    target_lane_id = lane_id[:-1] + '0'
    lane_position = vehicle.getLanePosition(subject_id)
    neighbors: List[Tuple[str, float]] = [
        (vehicle_id, vehicle.getLanePosition(vehicle_id)) for vehicle_id in lane.getLastStepVehicleIDs(target_lane_id)]
    neighbors = sorted(neighbors, key=lambda item: item[1])
    leaders = [x[0] for x in neighbors if x[1] >= lane_position]
    followers = [x[0] for x in neighbors if x[1] < lane_position]
    leader1 = leaders[0] if len(leaders) else None
    follower1 = followers[-1] if len(followers) else None
    follower2 = followers[-2] if len(followers) > 1 else None
    return leader1, follower1, follower2


def mclc(subject_id: str, N: int, tb: float) -> None:
    global eng

    lane_id = vehicle.getLaneID(subject_id)
    leader, follower1, follower2 = get_gcv(subject_id)

    ps = []
    if leader:
        prev = vehicle.getLeader(leader)
        ps.append((leader, prev[0] if prev else None))
        if not follower1:
            ps.append((None, leader))
    if follower1:
        prev = vehicle.getLeader(leader)
        ps.append((follower1, prev[0] if prev else None))
        if not follower2:
            ps.append((None, follower1))
    if follower2:
        prev = vehicle.getLeader(leader)
        ps.append((follower2, prev[0] if prev else None))

    rets = []
    x_sv = vehicle.getLanePosition(subject_id)
    v_sv = vehicle.getSpeed(subject_id)
    u_sv = vehicle.getAcceleration(subject_id)
    x_b = lane.getLength(lane_id)
    leader = vehicle.getLeader(subject_id)
    if leader:
        leader = leader[0]
        x_leader = vehicle.getLanePosition(leader)
        v_leader = vehicle.getSpeed(leader)
    else:
        x_leader = -1.0
        v_leader = -1.0

    for gcv, pv in ps:
        if gcv:
            x_gcv = vehicle.getLanePosition(gcv)
            v_gcv = vehicle.getSpeed(gcv)
        else:
            x_gcv = -1.0
            v_gcv = -1.0

        if pv:
            x_pv = vehicle.getLanePosition(pv)
            v_pv = vehicle.getSpeed(pv)
        else:
            x_pv = -1.0
            v_pv = -1.0

        n, u, agcv, cost = eng.mclc(x_sv, v_sv, u_sv, x_leader, v_leader, x_gcv, v_gcv, x_pv, v_pv, x_b, N, tb,
                                    nargout=4)
        rets.append((gcv if gcv else None, n, u, agcv, cost))

    if not len(rets):
        vehicle.setSpeed(subject_id, -1)
        return
    rets = sorted(rets, key=lambda x: (x[4], x[1]))
    gcv, n, u, agcv, c = rets[0]
    if c > 1e8:
        if lane.getLength(lane_id) - vehicle.getLanePosition(subject_id) < 50:
            vehicle.setSpeed(subject_id, min(v_sv, math.sqrt(
                2 * 2.5 * lane.getLength(lane_id) - vehicle.getLanePosition(subject_id))))
            vehicle.changeLaneRelative(subject_id, -1, 3)
        else:
            vehicle.setSpeed(subject_id, -1)
    else:
        vehicle.setSpeed(subject_id,
                         min(math.sqrt(2 * 1.5 * lane.getLength(lane_id) - vehicle.getLanePosition(subject_id)),
                             max(0, vehicle.getSpeed(subject_id) + u)))
        if gcv:
            vehicle.slowDown(gcv, max(0, vehicle.getSpeed(gcv) + agcv), 1)
        if n == 1:
            vehicle.changeLaneRelative(subject_id, -1, 1)


def get_gcv_in(subject_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    m = {'in1': 'round4_0', 'in2': 'round1_0', 'in3': 'round2_0', 'in4': 'round3_0'}
    lane_id = m[vehicle.getRoadID(subject_id)]
    subject_length = lane.getLength(vehicle.getLaneID(subject_id))
    target_length = lane.getLength(lane_id)
    lane_position = vehicle.getLanePosition(subject_id) - subject_length
    neighbors: List[Tuple[str, float]] = [(vehicle_id, vehicle.getLanePosition(vehicle_id) - target_length) for
                                          vehicle_id in lane.getLastStepVehicleIDs(lane_id) if
                                          vehicle.getRoute(vehicle_id)[-1][-1] != vehicle.getRoadID(subject_id)[-1]]
    neighbors = sorted(neighbors, key=lambda item: item[1])
    leaders = [x for x in neighbors if x[1] >= lane_position]
    followers = [x for x in neighbors if x[1] < lane_position]
    leader = leaders[0] if len(leaders) else None
    follower1 = followers[-1] if len(followers) else None
    follower2 = followers[-2] if len(followers) > 1 else None
    return leader, follower1, follower2


def mclc_in(subject_id: str, N: int, tb: float) -> None:
    global eng

    lane_id = vehicle.getLaneID(subject_id)
    leader, follower1, follower2 = get_gcv_in(subject_id)

    ps = []
    if leader:
        prev = vehicle.getLeader(leader)
        ps.append((leader, prev[0] if prev else None))
        if not follower1:
            ps.append((None, leader))
    if follower1:
        prev = vehicle.getLeader(leader)
        ps.append((follower1, prev[0] if prev else None))
        if not follower2:
            ps.append((None, follower1))
    if follower2:
        prev = vehicle.getLeader(leader)
        ps.append((follower2, prev[0] if prev else None))

    rets = []
    x_sv = vehicle.getLanePosition(subject_id) - lane.getLength(lane_id)
    v_sv = vehicle.getSpeed(subject_id)
    u_sv = vehicle.getAcceleration(subject_id)
    x_b = 0.
    leader = vehicle.getLeader(subject_id)
    if leader:
        leader = leader[0]
        x_leader = vehicle.getLanePosition(leader) - lane.getLength(lane_id)
        v_leader = vehicle.getSpeed(leader)
    else:
        x_leader = -1.0
        v_leader = -1.0

    for gcv, pv in ps:
        if gcv:
            x_gcv = vehicle.getLanePosition(gcv[0]) - lane.getLength(vehicle.getLaneID(gcv[0]))
            v_gcv = vehicle.getSpeed(gcv[0])
        else:
            x_gcv = -1.0
            v_gcv = -1.0

        if pv and vehicle.getRoute(pv[0])[-1][-1] != vehicle.getRoadID(subject_id)[-1]:
            x_pv = vehicle.getLanePosition(pv[0]) - lane.getLength(vehicle.getLaneID(pv[0]))
            v_pv = vehicle.getSpeed(pv[0])
        else:
            x_pv = -1.0
            v_pv = -1.0

        u, agcv, cost = eng.mclc_in(x_sv, v_sv, u_sv, x_leader, v_leader, x_gcv, v_gcv, x_pv, v_pv, x_b, N, tb,
                                    nargout=3)
        rets.append((gcv[0] if gcv else None, u, agcv, cost))
        
    if not len(rets):
        return
    rets = sorted(rets, key=lambda x: (x[3], x[1]))
    gcv, u, agcv, c = rets[0]
    if c > 1e8:
        if lane.getLength(lane_id) - vehicle.getLanePosition(subject_id) < 50:
            vehicle.setSpeed(subject_id, min(v_sv, math.sqrt(
                2 * 2.5 * (lane.getLength(lane_id) - vehicle.getLanePosition(subject_id)))))
        else:
            vehicle.setSpeed(subject_id, -1)
    else:
        vehicle.setSpeed(subject_id, max(0, vehicle.getSpeed(subject_id) + u))
        if False and gcv:
            vehicle.slowDown(gcv, max(0, vehicle.getSpeed(gcv) + agcv), 1)


def run(N: int = 20, tb: float = 3.0) -> None:
    """execute the TraCI control loop"""
    while simulation.getMinExpectedNumber() > 0:
        traci.simulationStep()
        if N == 0:
            continue
        for vehicle_id in vehicle.getIDList():
            vehicle.setSpeed(vehicle_id, -1)
        for in_lane in ('in1_0', 'in2_0', 'in3_0', 'in4_0'):
            for vehicle_id in lane.getLastStepVehicleIDs(in_lane):
                mclc_in(vehicle_id, N, tb)
        for current_lane, out_road in (
                ('round1_1', 'out2'), ('round2_1', 'out3'), ('round3_1', 'out4'), ('round4_1', 'out1')):
            for vehicle_id in lane.getLastStepVehicleIDs(current_lane):
                if vehicle.getRoute(vehicle_id)[-1] == out_road:
                    vehicle.setLaneChangeMode(vehicle_id, 512)
                    mclc(vehicle_id, N, tb)


eng = None


def main() -> None:
    global eng
    eng = matlab.engine.start_matlab()
    traci.start([checkBinary('sumo'), "-c", "road/road.sumocfg"])
    run(20, 3.0)
    traci.close()
    sys.stdout.flush()


# this is the main entry point of this script
if __name__ == "__main__":
    main()
