import sys
import json

class Pre_Processing:
    def __init__(self, input_settings):
        self.read_settings(input_settings)


    def read_settings(self, input_settings):
        with open(input_settings, "r") as file:
            self.settings = json.load(file)
    
