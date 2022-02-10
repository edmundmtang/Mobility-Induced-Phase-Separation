import propulsiveSensing as ps
import mysql.connector
from db_login import *
'''
cnx = mysql.connector.connect(
    host = host,
    user = user,
    password = password,
    database = database
    )

cnx.close()
'''

z1 = ps.zone('rectangle',[2,-1,0,1,1,1])
z2 = ps.zone('rectangle',[2,0,1,-1,1,2])
fld = ps.field([z1,z2])
fld.check([0,0,0])
