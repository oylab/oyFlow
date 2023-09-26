# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:58:53 2016

@author: Alonyan
"""
import os

from . import Utilities
from .Flow import Workspace
from . import Gating

os.environ["OMP_NUM_THREADS"] = "4"
