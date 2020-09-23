from collections import namedtuple
from time import time
import pandas as pd

#Data = namedtuple("data",["type","head","content","appendix","file"])
Task = namedtuple("task",["data_type", "out_file", "num_thread"])
class Data:
    def __init__(self, type:str, head:str, content:str, appendix:str, file_name:str):
        self.type = type
        self.head = head
        self.content = content
        self.appendix = appendix
        self.file_name = file_name
class Cell:
    def __init__(self, name:str, data:Data=None, log:str=""):
        #name, data exposed
        #log system
        self.name = name
        self.datas = {}
        if data != None:
            self.add_data(data)
        self.logs = []
        if log != "":
            self.add_log(log)
    def add_log(self, msg):
        #time stamp
        self.logs.append(msg)
    def get_log(self):
        return "\n".join(self.logs)
    def add_data(self, data):
        self.datas[data.type] = data
    def get_data(self, data_type:str):
        return self.datas[data_type]
