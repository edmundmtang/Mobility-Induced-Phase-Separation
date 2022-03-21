'''3-Dimensional Random Walk With Self-Propulsion

A set of functions to simulate self-propelling microparticles
and record the simulations to a mySQL database. This module
also contains functions for calculating auxiliary parameters
and visualizing the results.

Changelog
2022/02/15 - Updated RW_step() and RW_sim() to utilize fields
2022/02/09 - Added field class
2022/02/08 - Defined zone class for rectangles and circles
2022/02/05 - Imported code from 3DPropulsiveWalk

Edmund Tang 2022/02/05
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

class zone:
    '''Defines the speed modifier for a region of a
    simulation. A zone object checks if a position exists
    within the boundaries it defines, and outputs a speed
    modifier if it is true.'''

    def __init__(self,zoneType,parameters):
        validZones = ['rectangle', 'circle'] # list of defined zones
        try:
            if zoneType in validZones:
                self.zoneType = zoneType # parameters depend on the type
                if type(parameters) == list or type(parameters) == dict:
                    self.params = parameters # parameters are a list of values describing the zone
                else:
                    print("Error: Parameters are not a list or dictionary.")
            else:
                print("Error: Unexpected zoneType.")
        except:
            print("Error: Invalid input for zone class.")

    def check(self,pos): # determines the speed modifier for a valid position
        # new zone types need to be defined here
        # pos is the position, defined by a vector
        spdMod = 'undefined'
        if self.zoneType == 'rectangle':
            try:
                tmp_pos = np.delete(pos,self.params[0]) # 0th param defines orientation of the rectangular prism
                isContained = ((self.params[1] <= tmp_pos[0] <= self.params[2])
                               & (self.params[3] <= tmp_pos[1] <= self.params[4]))
                if isContained:
                    spdMod = self.params[5]
            except:
                print("Error: Attempted to execute a zone with invalid parameters.")
        elif self.zoneType == 'circle':
            try:
                tmp_pos = np.delete(pos,self.params[0]) # 0th param defines orientation of the circular prism
                isContained = ((tmp_pos[0] - self.params[1])**2 + (tmp_pos[1] - self.params[2])**2
                               <= self.params[3]**2)
                if isContained:
                    spdMod = self.params[4]
            except:
                print("Error: Attempted to execute a zone with invalid parameters.")
        else:
            print("Error: Attempted to execute an invalid zone")
        
        return isContained, spdMod
        
class field:
    '''Defines how the speed of a particle varies with
    location. A field is a list of zones contained in
    self.zones. Use list methods to build this list. 
    '''

    def __init__(self, bounds = [], spdDefault = 0, zones = []):
        # bounds = [xMin, xMax, yMin, yMax, zMin, zMax]
        self.bounds = bounds
        self.spdDefault = spdDefault
        isZones = all(isinstance(i, zone) for i in zones)
        if (zones == []) or isZones:
            self.zones = zones # zones is a list of zone objects
        else:
            print("Error: A valid list of zones was not provided.")

    def check(self,pos):
        # checks that current position is within bounds, corrects if not, then
        # evaluates the speed modifier corresponding to the input position
        if self.bounds != []:
            if len(self.bounds) != 6:
                print("Error: Incorrect number of coordinates provided for the bounds of the simulation.")
            if pos[0] < self.bounds[0]: pos[0] = self.bounds[0]
            elif pos[0] > self.bounds[1]: pos[0] = self.bounds[1]
            if pos[1] < self.bounds[2]: pos[1] = self.bounds[2]
            elif pos[1] > self.bounds[3]: pos[1] = self.bounds[3]
            if pos[2] < self.bounds[4]: pos[2] = self.bounds[4]
            elif pos[2] > self.bounds[5]: pos[2] = self.bounds[5]

        zoneID = 1 # zoneID = 0 represents no zone
        for zone in self.zones:
            isContained, spdMod = zone.check(pos)
            if isContained:
                return pos, zoneID, spdMod
            zoneID += 1
        return pos, 0, self.spdDefault

class RW_step:
    '''Iterator that yields positions and orientations in a
    Brownian, self-propelling particle assuming a non-
    dimensionalized diffusivity of 4/3.'''

    def __init__(self,dt,v0,fld):
        self.dt = dt # dt is the time step size
        self.v0 = v0 # V is the non-dimensional propulsive speed (scalar)
        self.fld = fld # field that the simulation exists in

    def __iter__(self):
        next_pos = np.array([0, 0, 0]) # initial position is [0, 0, 0]
        next_pos, zoneID, spdMod = self.fld.check([0,0,0])
        self.pos = next_pos

        self.zone = zoneID # initial zone
        self.v = spdMod*self.v0 # initial speed

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
        next_pos = self.pos + vTot*dt
        # calculate rotation
        tB = np.array([(2/dt)**(1/2)*random.normal() for _ in range(3)])
        tTot = tB # total torque
        q = self.ori
        q = q + dt*np.cross(tTot,q) # determine effect on orientation
        next_ori = q / np.sqrt(np.sum(q**2)) # normalize orientation vector
        # correct position and determine speed of next step
        next_pos, zoneID, spdMod = self.fld.check(next_pos)
        self.pos = next_pos
        self.zone = zoneID
        self.v = spdMod*self.v0
        self.ori = next_ori

def RW_sim(simLen, stepsPerObs, dt, v0, fld, cnx):
    '''Function that performs the 3-D random walk with self-
    propulsion simulation according to the specified
    parameters. The function also commits the results to a
    database and returns a simID.'''

    cursor = cnx.cursor() # get cursor from the mySQL connection
    # create experiment entry in database
    ins_stmt = """INSERT INTO `experiments` (`simLen`,`stepsPerObservation`,`stepSize`,`baseSpeed`,`spdModDefault`,`xMin`,`xMax`,`yMin`,`yMax`,`zMin`,`zMax`)
                  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""
    ins_data = [simLen, stepsPerObs, dt, v0, fld.spdDefault]
    for i in fld.bounds:
        ins_data.append(i) # add bounds to insertion data
    cursor.execute(ins_stmt, ins_data)
    # determine new simulation ID
    simID = cursor.lastrowid
    print("\n==========\nSimulation ID: {}".format(simID))
    print("Simulation Length: {}\nSteps per Observation: {}\nTime Step Size: {}\nBase Speed: {}".format(simLen,stepsPerObs,dt,v0))
    print("\nRunning simulation...")
        
    sim = RW_step(dt, v0, fld) # create iterator
    iter(sim)
    
    # create entry for the zones in this experiment
    ins_stmt = """INSERT INTO `zones` (`simID`, `zoneID`, `type`, `parameters`, `spdMod`)
                  VALUES (%s, %s, %s, %s, %s)"""
    zoneID = 1
    for zone in fld.zones:
        delim = ","
        zoneSpdMod = zone.params[-1] # spdMod, note that it's redundant in params
        zoneParams = delim.join(list(map(str,zone.params))) # convert params into str
        ins_data = [simID, zoneID, zone.zoneType, zoneParams, zoneSpdMod]
        cursor.execute(ins_stmt, ins_data)
        zoneID += 1
    # create first entry in trajectory table for this experiment
    ins_stmt = """INSERT INTO `trajectories` (`simID`, `obsNum`, `zoneID` , `xpos`, `xori`, `ypos`, `yori`, `zpos`, `zori`)
                  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"""
    obsNum = 0
    ins_data = [simID, obsNum, sim.zone]
    
    for i in range(3):
        ins_data.append(float(sim.pos[i]))
        ins_data.append(float(sim.ori[i]))
    
    cursor.execute(ins_stmt, ins_data)
    for _ in range(simLen):
        for _ in range(stepsPerObs): # perform a certain number of steps between each observation
            next(sim)
        obsNum += 1
        ins_data = [simID, obsNum, sim.zone]
        for i in range(3):
            ins_data.append(float(sim.pos[i]))
            ins_data.append(float(sim.ori[i]))
        cursor.execute(ins_stmt, ins_data)
    cnx.commit()
    print("Simulation complete!")
    return simID

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

def readField(expNum, cnx):
    '''Recreate a field object by reading it from the
    SQL database. Does so for a specific experiment.'''
    cursor = cnx.cursor() # get cursor from the mySQL connection
    cursor.execute(f"""SELECT `spdModDefault`,`xMin`,`xMax`,`yMin`,`yMax`,`zMin`,`zMax` from `experiments` WHERE
                    (`simID` = {expNum})""")
    results = cursor.fetchone()
    spdDefault = results[0]
    bounds = results[1:7]
    cursor.execute(f"""SELECT * FROM `zones` WHERE
                    (`simID` = {expNum})""")
    results = cursor.fetchall()
    zones = []
    for row in results:
        zoneType = row[2]
        parameters = list(map(int,row[3].split(","))) # read the parameters str and convert to list of ints
        zones.append(zone(zoneType,parameters))       
    fld = field(bounds,spdDefault,zones)
    return fld
    
def plotTrajectory2D(df, fld,
                     steps, excludeDim = 2,
                     filename = "scatter2D_animated.html",
                     xrange = [-1, 1], yrange = [-1, 1]):
    '''Create an animated 2D plot of a set of trajectories
    with lines indicating zones. Takes pandas data frame
    and a field object as inputs.'''

    if excludeDim == 0:
        r1 = "ry"; r2 = "rz"; xmin = fld.bounds[2]; xmax = fld.bounds[3]; ymin = fld.bounds[4]; ymax = fld.bounds[5]
    elif excludeDim == 1:
        r1 = "rx"; r2 = "rz"; xmin = fld.bounds[0]; xmax = fld.bounds[1]; ymin = fld.bounds[4]; ymax = fld.bounds[5]
    elif excludeDim == 2:
        r1 = "rx"; r2 = "ry"; xmin = fld.bounds[0]; xmax = fld.bounds[1]; ymin = fld.bounds[2]; ymax = fld.bounds[3]

    # test that columns in df are properly named
    columnNames = ["Simulation ID", "obsNum", r1, r2]
    try:
        for name in columnNames:
            df[name]
    except:
        print(f"""ValueError: Could no find column. Expected one of {columnNames} in the data frame.""")
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

    # make data (zones)
    for i, zone in zip(range(len(fld.zones)), fld.zones):
        if zone.zoneType == 'rectangle':
            p = zone.params[1:5]
            xList1 = [p[0],p[0],p[1]]
            yList1 = [p[2],p[3],p[3]]
            xList2 = [p[0],p[1],p[1]]
            yList2 = [p[2],p[2],p[3]]
        trace_z0 = go.Scatter(x=xList1,
                            y=yList1,
                            mode="lines",
                            name = "Zone {}".format(i),
                            line_color = px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)]
                            )
        trace_z1 = go.Scatter(x=xList2,
                            y=yList2,
                            mode="lines",
                            name = "Zone {}".format(i),
                            fill='tonexty',
                            line_color = px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)]
                            )
        fig_dict["data"].append(trace_z0)
        fig_dict["data"].append(trace_z1)
    # make data (bounds)
    xBoundary = [xmin, xmin, xmax, xmax, xmin]
    yBoundary = [ymin, ymax, ymax, ymin, ymin]
    trace_B = go.Scatter(x=xBoundary,
                         y=yBoundary,
                         mode="lines",
                         line_color="black",
                         name = "Boundary"
                         )
    fig_dict["data"].append(trace_B)
    # make data (particles)
    data_by_step = df[df["obsNum"] == 0]
    traceP = go.Scatter(x=list(data_by_step[r1]),
                        y=list(data_by_step[r2]),
                        mode = "markers",
                        marker = dict(
                            color = 'black',
                            size = 4
                            ), 
                        name = "Particle {}".format(simID)
                        )
    fig_dict["data"].append(traceP)

    # make frames
    frames = []
    for step in steps:
        frame_data = []
        # make frame for zones
        for zone in fld.zones:
            if zone.zoneType == 'rectangle':
                p = zone.params[1:5]
                xList1 = [p[0],p[0],p[1]]
                yList1 = [p[2],p[3],p[3]]
                xList2 = [p[0],p[1],p[1]]
                yList2 = [p[2],p[2],p[3]]
            frame_z0 = dict(type = "scatter",
                            x=xList1,
                            y=yList1
                            )
            frame_z1 = dict(type = "scatter",
                            x=xList2,
                            y=yList2
                            )
            frame_data.append(frame_z0)
            frame_data.append(frame_z1)
        # make frame for boundary
        frame_B = dict(type = "scatter",
                       x=xBoundary,
                       y=yBoundary
                       )
        frame_data.append(frame_B)

        # make frame for particles
        data_by_step = df[df["obsNum"] == step]
        frameP = dict(type = "scatter",
                      x=list(data_by_step[r1]),
                      y=list(data_by_step[r2]),
                      )
        frame_data.append(frameP)

        frame = dict(data = frame_data,
                     traces = list(range(0,2+len(fld.zones)*2)),
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

    fig.update_xaxes(range=xrange,
                     constrain="domain")
    fig.update_yaxes(range=yrange,
                     scaleanchor = "x",
                     scaleratio = 1
                     )
    

    fig.write_html(filename, auto_open=True)
    print("Figure generated!")
    print(f"Filename: {filename}")

    return fig, fig_dict
                     
def plotTrajectory3D(df, steps, mode = "position",
                     fileName = "trailing_scatterplot.html",
                     xrange = [-20, 20], yrange = [-20, 20], zrange = [-20, 20]):
    '''Create an animated 3D plot of a set of trajectories
    with trailing lines representing where each particle has
    gone. Takes a pandas data frame as input.'''

    # check mode
    if not (mode.lower() == 'position' or mode.lower() == 'orientation'):
        print("Error: Expected either 'position' or 'orientation' for mode.")
        return

    # test that columns in df are properly named
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









    
