import propulsiveSensing as ps
import mysql.connector
from db_login import *
# Connect to database
cnx = mysql.connector.connect(
    host = host,
    user = user,
    password = password,
    database = database
    )
# Simulation parameters
simLen = 2000 # simulation length
stepsPerObs = 10**3 # steps in the simulation between recorded observations
dt = 10**-5 # time between simulation steps
v0 = 5 # base propulsive speed
# Build the simulation field
z1 = ps.zone('rectangle',[2,-1,0,-1,1,0])
z2 = ps.zone('rectangle',[2, 0,1,-1,1,1])
fld = ps.field(bounds = [-1, 1,-1,1,0,0], zones = [z1,z2])
# Run simulations
for i in range(128):
    ps.RW_sim(simLen,stepsPerObs,dt,v0,fld,cnx)
cnx.close()
