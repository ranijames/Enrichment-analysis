#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import yaml
from sqlalchemy import *
import pymysql

def get_connection(settings):
    """
    Creates database connection.

    @param: settings: dictionary loaded from a yaml file with the necessary details for the connection with unidb.
    """
    host = settings["database_credentials"]["host"]
    user = settings["database_credentials"]["user"]
    passwd = settings["database_credentials"]["password"]
    database = settings["database_credentials"]["database"]  
    print("db:", host, "user:", user, "database:", database)

    # Database connection
    db=pymysql.connect(
            host=host,
            user=user,
            passwd=passwd,
            db=database,
            autocommit=True)

    engine = create_engine(f'mysql+pymysql://{user}:{passwd}@{host}/{database}')
    print(f'mysql+pymysql://{user}:{passwd}@{host}/{database}')

    cursor = db.cursor()
    conn = engine.connect()
    metadata = MetaData()
    return cursor, conn

