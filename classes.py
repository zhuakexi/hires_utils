import pandas as pd

class Cell:
    def __init__(self, name:str, appendix, data:"dataframe", file_head:str="", log:str=[]):
        #name, data exposed
        #log system
        self.name = name
        self.data = data
        self.appendix = appendix
        self.log = log
        self.file_head = file_head
    def add_log(self, msg):
        #time stamp
        self.log.append(msg)
    def get_log(self):
        return "\n".join(self.log)