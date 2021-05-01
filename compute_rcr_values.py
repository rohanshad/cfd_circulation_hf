#!/usr/bin/python

import os, re
import numpy as np 
import matplotlib.pyplot as plt

debug = True

def compute_rcr_parameters(area, Q_goal, P_min, P_max, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor=1.0):

    tol = 1e-12

    # total resistances 
    R_total = []
    for Q in Q_goal:
        R_total.append(P_mean/Q)

    if debug:
        print(("R_total = ", R_total))

    # ratio of resistors is constant 
    # resistors sum to total resistance 
    R_d = []
    for R_t in R_total:
        R_d.append(R_t / (1.0 + ratio_prox_to_distal_resistors))
        
    R_p = []
    for R_dist, R_t in zip(R_d, R_total):
        R_p.append(R_t - R_dist) 

    for R_prox, R_dist, R_t in zip(R_p, R_d, R_total):
        assert abs(R_dist + R_prox - R_t) < tol 

    C = []
    for R_dist in R_d:
        C.append(-C_prefactor * decay_time / (R_dist * np.log(P_min/P_max)))

    return R_p, C, R_d, R_total



def read_flow_file(file_name): 
    # times in first column, flows in second 
    # no checks for uncompliant files listed 
    times = []
    flows = []
    with open(file_name, 'r') as f:
        for line in f:  
            times.append(float(line.split()[0]))
            flows.append(float(line.split()[1]))

    n_times = len(times)

    beat_time = times[n_times-1]

    if flows[0] != flows[n_times-1]:
        raise ValueError('Expecting periodic flow, did not get it')

    # remove last element, do not count it twice 
    times.pop()
    flows.pop()

    plots = False
    if plots: 
        print("should be plotting... ")
        plt.plot(times, flows)
        plt.xlabel('time (s)')
        plt.ylabel('flow (ml/s)')
        plt.show()

    t_at_min_flow = times[np.argmin(flows)]
    t_at_max_flow = times[np.argmax(flows)]

    assert t_at_max_flow > t_at_min_flow

    Q_min  = np.min(flows)  
    Q_max  = np.max(flows)  
    Q_mean = np.mean(flows)
    
    return beat_time, Q_min, Q_max, Q_mean, t_at_min_flow, t_at_max_flow


if __name__== "__main__":

    if debug:
        print("Debug info:")

    tol = 1e-12

    MMHG_TO_CGS = 1333.22368

    # four outlets 

    P_sys = 82.66666667 * MMHG_TO_CGS
    P_dia = 63.5 * MMHG_TO_CGS
    P_mean = 0.7*P_dia + 0.3*P_sys # tune this value 

    C_prefactor = 0.1 # tune this value, 
                      # 1 for no adjustment 
                      # lower for faster response 
                      # higher for slower response 

    P_min = P_dia
    P_max = P_sys

    file_name = '97_test.flow'
    beat_time, Q_min, Q_max, Q_mean, t_at_min_flow, t_at_max_flow = read_flow_file(file_name)

    decay_time = 0.6 * beat_time # tune this parameter, should be duration of diastole 

    heart_rate = 1/beat_time # beats per second 
    if debug: 
        print(("heart_rate = ", heart_rate, "beats per second, ", 60*heart_rate, " bpm"))

        print(("P_dia, P_sys, P_mean = ", P_dia, P_sys, P_mean))
        print(("Q_min, Q_max, Q_mean = ", Q_min, Q_max, Q_mean))

        print(("t_at_min_flow = ", t_at_min_flow, ", t_at_max_flow = ", t_at_max_flow))

    dt_transition = t_at_max_flow - t_at_min_flow
    timescale_rd_c = 2.0*dt_transition

    names = ['innominate', 'L_carotid', 'L_subclavian', 'aorta']
    area = np.array([2.03203, 0.327643, 0.633362, 2.17127])

    Q_goal = np.zeros(4)

    Q_goal[3] = .7  # https://doi.org/10.1161/JAHA.115.002657, fig 4 bottom right
                            # .7 of total flow to descending aorta 
    Q_goal[0] = .15 # right gets half of remaining flow to upper body 

    area_L = area[1] + area[2] # total area of L_carotid, L_subclavian
    Q_goal[1] = (area[1]/area_L) * .15 
    Q_goal[2] = (area[2]/area_L) * .15 

    # scale by total flow 
    Q_goal *= Q_mean

    if debug:
        print(("Q_goal = ", Q_goal))

    # quick sanity check 
    assert abs(Q_mean - np.sum(Q_goal)) < tol

    ratio_prox_to_distal_resistors = 77.0 / 1185.0 # constant from https://www.physiology.org/doi/10.1152/jappl.1990.69.1.112 

    R_p, C, R_d, R_total = compute_rcr_parameters(area, Q_goal, P_min, P_max, P_mean, Q_mean, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

    if debug:
        print("\n\n\n")

    print("name,\tr_p,\tC,\t,R_d,")
    for name, R_p_tmp, C_tmp, R_d_tmp in zip(names, R_p, C, R_d):
        print((name + "\t" + str('{:.6f}'.format(R_p_tmp)) + "," + str('{:.6f}'.format(C_tmp)) + "," + str('{:.6f}'.format(R_d_tmp))))



















