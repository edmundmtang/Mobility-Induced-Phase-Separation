import propulsiveSensing as ps
import mysql.connector
import pandas as pd
from db_login import *
# Connect to database
cnx = mysql.connector.connect(
    host = host,
    user = user,
    password = password,
    database = database
    )
cursor = cnx.cursor()

obsLimit = 1000

fld = ps.readField(1,cnx)

# filter for experiments of interest
cursor.execute(f"""SELECT `simID` FROM `experiments` WHERE
                (`simID` <= 200)""") # change the number here to limit experiments in final plot 
res1 = cursor.fetchall()
simIDs = [x for t in res1 for x in t]

# pick out results for experiments of interest
allResults = []
for simID in simIDs:
    cursor.execute(f"""SELECT `simID`, `obsNum`, `xpos`, `ypos`, `zpos`
                    FROM `trajectories` WHERE (`simID` = {simID})""")
    results = cursor.fetchall()
    for row in results:
        allResults.append(row)
x = [[ij for ij in i] for i in allResults]
df = pd.DataFrame([[ij for ij in i] for i in allResults]) 
df.rename(columns = {0:"Simulation ID", 1:"obsNum", 2:"rx", 3:"ry", 4:"rz"}, inplace = True)

steps = [x for x in list(range(obsLimit+1)) if x % 1 == 0]

fig, fig_dict = ps.plotTrajectory2D(df,fld,steps,excludeDim = 2, filename = "test.html",
                          xrange=[-1.05,1.05],yrange=[-1.05,1.05])

cnx.close()
