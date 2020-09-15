from collections import namedtuple
from time import time
import pandas as pd

Data = namedtuple("data",["type","head","content","appendix","file"])

class Cell:
    def __init__(self, name:str, data:Data=None, log:str=""):
        #name, data exposed
        #log system
        self.name = name
        if data == None:
            self.datas = {}
        else:
            self.add_data(data)
        if log == "":
            self.logs = []
        else:
            self.add_log(log)
    def add_log(self, msg):
        #time stamp
        self.logs.append(msg)
    def get_log(self):
        return "\n".join(self.logs)
    def add_data(self, data):
        self.datas[data.type] = data
    