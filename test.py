import propulsiveSensing as ps
import mysql.connector
from db_login import *

cnx = mysql.connector.connect(
    host = host,
    user = user,
    password = password,
    database = database
    )

cnx.close()
