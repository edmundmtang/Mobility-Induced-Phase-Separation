'''3-Dimensional Random Walk With Self-Propulsion

A set of functions to simulate self-propelling microparticles
and record the simulations to a mySQL database. This module
also contains functions for calculating auxiliary parameters
and visualizing the results.

Changelog
2021/04/18 - base version
2021/04/24 - began work on 3d plotting function
2021/04/26 - redid 3d plotting function for trajectories
2021/04/30 - added regression function for 2d msd

Edmund Tang 2021/05/01
'''

import numpy as np
from numpy import random
import numpy.linalg as la
import mysql.connector
import plotly.express as px
import plotly.graph_objs as go
from scipy.optimize import curve_fit
import pandas as pd
from lmfit import Model

class RW_step:
    '''Iterator that yields positions and orientations in a
    Brownian, self-propelling particle assuming a non-
    dimensionalized diffusivity of 4/3.'''

    def __init__(self,dt,v):
        self.dt = dt # dt is the time step size
        self.v = v # V is the non-dimensional propulsive speed (scalar)

    def __iter__(self):
        self.pos = np.array([0, 0, 0]) # initial position is [0, 0, 0]
        q = random.rand(3)
        q_norm = q / np.sqrt(np.sum(q**2))
        self.ori = q_norm # initial orientation is random unit vector
        return self

    def __next__(self): # apply changes in position and orientation
        dt = self.dt
        # calculate translation
        vP = self.v * self.ori # propulsive velocity, in direction of orientation
        vB = np.array([(8/(3*dt))**(1/2)*random.normal() for _ in range(3)])
        vTot = vP + vB # total velocity
        self.pos = self.pos + vTot*dt
        # calculate rotation
        tB = np.array([(2/dt)**(1/2)*random.normal() for _ in range(3)])
        tTot = tB # total torque
        q = self.ori
        q = q + dt*np.cross(tTot,q) # determine effect on orientation
        self.ori = q / np.sqrt(np.sum(q**2)) # normalize orientation vector

def RW_sim(simLen, stepsPerObs, dt, v, cnx):
    '''Function that performs the 3-D random walk with self-
    propulsion simulation according to the specified
    parameters. The function also commits the results to a
    database and returns a simID.'''

    cursor = cnx.cursor() # get cursor from the mySQL connection
    # create experiment entry in database
    ins_stmt = """INSERT INTO `experiments` (`simLen`,`stepsPerObservation`, `stepSize`, `propulsiveSpeed`)
                  VALUES (%s, %s, %s, %s)"""
    ins_data = (simLen, stepsPerObs, dt, v)
    cursor.execute(ins_stmt, ins_data)
    # determine new simulation ID
    newSimID = cursor.lastrowid
    print("\n==========\nSimulation ID: {}".format(newSimID))
    print("Simulation Length: {}\nSteps per Observation: {}\nTime Step Size: {}\nPropulsive Speed: {}".format(simLen,stepsPerObs,dt,v))
    print("\nRunning simulation...")
        
    sim = RW_step(dt, v) # create iterator
    iter(sim)
    # create first entry in trajectory table for this experiment
    ins_stmt = """INSERT INTO `trajectories` (`simID`, `obsNum`, `xpos`, `xori`, `ypos`, `yori`, `zpos`, `zori`)
                  VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"""
    obsNum = 0
    ins_data = [newSimID, obsNum]
    for i in range(3):
        ins_data.append(float(sim.pos[i]))
        ins_data.append(float(sim.ori[i]))
    cursor.execute(ins_stmt, ins_data)
    for _ in range(simLen):
        for _ in range(stepsPerObs): # perform a certain number of steps between each observation
            next(sim)
        obsNum += 1
        ins_data = [newSimID, obsNum]
        for i in range(3):
            ins_data.append(float(sim.pos[i]))
            ins_data.append(float(sim.ori[i]))
        cursor.execute(ins_stmt, ins_data)
    cnx.commit()
    print("Simulation complete!")
    return newSimID

def calc_angles(simID, cnx):
    '''Calculates the angles from the orientation vectors
    of a simulation with a particular simID'''

    cursor = cnx.cursor() # get cursor from the mySQL connection
    
    # obtain list of all entries that will be updated
    cursor.execute("SELECT `entryID` FROM `trajectories` WHERE `simID` = {}".format(simID))
    entryIDs = [item for t in cursor.fetchall() for item in t]
    ins_stmt = ("UPDATE `trajectories` SET `theta` = %s, `psi` = %s WHERE `entryID` = %s")
    
    for entryID in entryIDs:
        # calculate theta and psi
        cursor.execute("SELECT `xori`, `yori`, `zori` FROM `trajectories` WHERE `entryID` = {}".format(entryID))
        [rx, ry, rz] = [item for t in cursor.fetchall() for item in t]
        theta = np.arccos(rz)
        psi = np.arctan2(rx,ry)
        # update values to database
        ins_data = (theta, psi, entryID)
        cursor.execute(ins_stmt, ins_data)

    cnx.commit()

def calc_MSD(simID, cnx):
    '''Calculate the mean squared displacement from the
    trajectory data of a simulation given a particular
    simulation ID. Do so for both position and orientation'''
    
    cursor = cnx.cursor()
    # get some experiment parameters
    cursor.execute("SELECT `simLen`, `stepsPerObservation`, `stepSize` FROM `experiments` WHERE `simID` = {}".format(simID))
    [simLen, stepsPerObs, stepSize] = [item for t in cursor.fetchall() for item in t]
    print("\n==========\nSimulation ID: {}\nSimulation Length: {}".format(simID, simLen))
    print("Steps per Observation: {}\nTime Step Size: {}".format(stepsPerObs, stepSize))
    # get list of entry IDs
    cursor.execute("SELECT `entryID` FROM `trajectories` WHERE `simID` = {}".format(simID))
    entryIDs = [item for t in cursor.fetchall() for item in t]
    if entryIDs == []:
        print("\nSimulation (simID = {}) does not exist.".format(simID))
        print("Skipping calculation of MSD.")
        return
    # check that calculation was not already performed
    cursor.execute("SELECT COUNT(*) FROM `MSDs` WHERE `simID` = {} LIMIT 1".format(simID))
    msdExist = cursor.fetchone()[0]
    if msdExist != 0:
        print("\nMSDs have already been calculated for simID = {}".format(simID))
    else:
        print("\nCalcuating mean squared displacements...")
        # extract position data
        cursor.execute("SELECT `xpos`, `ypos`, `zpos`, `xori`, `yori`, `zori` FROM `trajectories` WHERE `simID` = {}".format(simID))
        data = cursor.fetchall()
        r = [(item[0], item[1], item[2]) for item in data]
        q = [(item[3], item[4], item[5]) for item in data]
        ins_stmt = ("""INSERT INTO `MSDs` (`simID`, `dt`, `msd`, `msd2d`, `msad`)
                       VALUES (%s, %s, %s, %s, %s)""")
        for i, entryID in zip(range(simLen),entryIDs):
            dt = i*stepsPerObs*stepSize
            if i == 0:                
                ins_data = (simID, dt, 0, 0, 0) # first entry is always the same
            else:
                sdList = []
                sdList2d = []
                sadList = []
                for j in range(simLen-i):
                    # pull position data & calculate squared displacements
                    [rx1, ry1, rz1] = r[j]
                    [rx2, ry2, rz2] = r[j+i]
                    drx = rx2 - rx1
                    dry = ry2 - ry1
                    drz = rz2 - rz1
                    sdList.append(drx**2 + dry**2 + drz**2) # 3D
                    sdList2d.append(drx**2 + dry**2) # 2D
                    # pull orientation data & calculate squared angular displacements
                    q1 = np.array(q[j])
                    q2 = np.array(q[j+i])
                    ad = np.arccos((np.dot(q1,q2))/(la.norm(q1)*la.norm(q2)))
                    sadList.append(ad**2)
                # calculate msd and msad
                msd = np.mean(sdList)
                msd2d = np.mean(sdList2d)
                msad = np.mean(sadList)
                ins_data = (simID, dt, msd, msd2d, msad)
            # add results to database
            cursor.execute(ins_stmt, ins_data)
    cnx.commit()
    print("Calculations complete!")
    return simID





   
def plotTrajectory3D(df, steps, mode = "position",
                     fileName = "rw3d_trailing_scatterplot.html",
                     xrange = [-20, 20], yrange = [-20, 20], zrange = [-20, 20]):
    '''Create an animated 3D plot of a set of trajectories
    with trailing lines representing where each particle has
    gone. Takes a pandas data frame as input.'''

    # check mode
    if not (mode.lower() == 'position' or mode.lower() == 'orientation'):
        print("Error: Expected either 'position' or 'orientation' for mode.")
        return

    # test that columns are properly named
    if mode.lower() == 'position':
        xstr = "rx"; ystr = "ry"; zstr = "rz"
    elif mode.lower() == 'orientation':
        xstr = "qx"; ystr = "qy"; zstr = "qz"
    columnNames = ["Simulation ID", "obsNum", xstr, ystr, zstr]
    try:
        for name in columnNames:
            df[name]
    except:
        print("ValueError: Could not find column. Expected one of ['Simulation ID', 'obsNum', 'rx', 'ry', 'rz', 'qx', 'qy', 'qz'] in the data frame.")
        return
    
    simIDs = [] # build list of all unique simulation IDs
    for simID in list(df["Simulation ID"]):
        if simID not in simIDs:
            simIDs.append(simID)

    # make figure
    fig_dict = {
        "data": [],
        "layout": {},
        "frames": []
    }

    # fill in most of layout
    fig_dict["layout"]["hovermode"] = "closest"
    fig_dict["layout"]["width"] = 1300
    fig_dict["layout"]["height"] = 900
    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 0, "redraw": True},
                                    "fromcurrent": True, "transition": {"duration": 5,
                                                                        "easing": "quadratic-in-out"}}],
                    "label": "Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                      "mode": "immediate",
                                      "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": True,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "Step:",
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 50, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }

    # make data

    for i, simID in zip(range(len(simIDs)), simIDs):
        # generate each pair of traces by simulation ID
        data_by_simID = df[df["Simulation ID"] == simID]
        # trace for the marker
        trace0 = go.Scatter3d(x=[list(data_by_simID[xstr])[0]],
                              y=[list(data_by_simID[ystr])[0]],
                              z=[list(data_by_simID[zstr])[0]],
                              mode = "markers",
                              marker = dict(
                                  color = px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)], #https://plotly.com/python/discrete-color/
                                  size = 4
                                  ), 
                              name = "Particle {}".format(simID)
                              )
        # trace for the line
        trace1 = go.Scatter3d(x=trace0.x,
                              y=trace0.y,
                              z=trace0.z,
                              mode = "lines",
                              line = dict(
                                  color = trace0.marker.color,
                                  width = 3),
                              name = "Trajectory {}".format(simID)
                              )
        # add traces to data
        fig_dict["data"].append(trace0); fig_dict["data"].append(trace1)

    # make frames
    frames = []
    for step in steps:
        frame_data = []
        for simID in simIDs:
            data_by_simID = df[df["Simulation ID"] == simID]
            frame_mrkr_by_simID = dict(type = "scatter3d", # end marker
                                       x=[list(data_by_simID[xstr])[step]],
                                       y=[list(data_by_simID[ystr])[step]],
                                       z=[list(data_by_simID[zstr])[step]]
                                       )
            frame_data.append(frame_mrkr_by_simID)
            frame_traj_by_simID = dict(type = "scatter3d", # trajectory line
                                       x=list(data_by_simID[xstr])[:step+1],
                                       y=list(data_by_simID[ystr])[:step+1],
                                       z=list(data_by_simID[zstr])[:step+1]                                   
                                       )
            frame_data.append(frame_traj_by_simID)
        frame = dict(data= frame_data,
                     traces = list(range(0,len(simIDs)*2)),
                     name = step)
        frames.append(frame)
                               
    fig_dict["frames"] = frames

    # slider frames
    for step in steps:  
        slider_step = {"args": [
            [step],
            {"frame": {"duration": 50, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 50}}
        ],
            "label": step,
            "method": "animate"}
        sliders_dict["steps"].append(slider_step)


    fig_dict["layout"]["sliders"] = [sliders_dict]

    fig = go.Figure(fig_dict)

    fig.update_layout(scene_aspectmode="cube")
    if mode.lower() == 'orientation':
        xrange = [-1, 1]; yrange = [-1, 1]; zrange = [-1, 1]
    fig.update_scenes(
        xaxis_autorange=False,
        yaxis_autorange=False,
        zaxis_autorange=False,
        xaxis_range=xrange,
        yaxis_range=yrange,
        zaxis_range=zrange
        )

    fig.write_html(fileName, auto_open=True)
    print("Figure generated!")
    print(f"Filename: {fileName}")
    return fig 

def explinear2d(t, D, V, tau):
    '''Expolinear curve typically used to describe the MSD
    of a propulsive particle'''
    
    return 4*D*t + (V**2*tau**2)/3*(2*t/tau + np.exp(-2*t/tau)-1)

def msdFit2d(simIDs, p0, cnx, fitLim=1.5):
    '''Perform a fit on the 2D MSD from a set of
    trajectories and record the results to the database.'''

    cursor = cnx.cursor()
    allResults = []
    for simID in simIDs:
        cursor.execute(f"SELECT `simID`, `dt`, `msd2d` FROM `MSDs` WHERE (simID = {simID} AND dt <= {fitLim})")
        result = cursor.fetchall()
        for row in result:
            allResults.append(row)

    df = pd.DataFrame([[ij for ij in i] for i in allResults])
    df.rename(columns = {0:"Simulation ID", 1:"dt", 2:"msd"}, inplace = True)

    t_data = np.array(df["dt"])    
    msd_data = np.array(df["msd"])
    
    explinModel = Model(explinear2d)
    result = explinModel.fit(msd_data, t=t_data, D=p0[0], V=p0[1], tau=p0[2])

    cursor.execute(f"""INSERT INTO `msdfit` (`fitLim`, `fitD`, `fitV`, `fitTau`)
                    VALUES ({fitLim},{result.params['D'].value},{result.params['V'].value},{result.params['tau'].value})
                    """)
    newSetID = cursor.lastrowid    
    for simID in simIDs:
        cursor.execute(f"""INSERT INTO `msdfitset` (`setID`, `simID`)
                        VALUES ({newSetID},{simID})
                        """)

    cnx.commit()
    return result









    
