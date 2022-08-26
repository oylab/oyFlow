# Metadata module for flow cytometry data
# Currently supports FCS files

from os import walk, listdir, path
import os
from os.path import join, isdir
import sys

import pandas as pd
import numpy as np
import dill as pickle
from ast import literal_eval
import warnings
import logging
from natsort import natsorted, natsort_keygen
from FlowCytometryTools import FCMeasurement


md_logger = logging.getLogger(__name__)
md_logger.setLevel(logging.DEBUG)

usecolslist = []


class Metadata(object):
    def __init__(self, pth):
        if path.isdir(pth):
            self.base_pth = pth
        else:
            self.base_pth, self._md_name = path.split(pth)

        self.type='Flow'
        self._load_metadata_FCS_GUI(self.base_pth)



#functions for generic loading!

    def _load_metadata_FCS_GUI(self,pth=''):
        GlobText = self._getPatternsFromPathGUI(pth=pth)
        box1 = self._getMappingFromGUI(GlobText)


    def _getPatternsFromPathGUI(self, pth=''):
        from ipywidgets import Button, Text, widgets, HBox, Layout
        from tkinter import Tk, filedialog
        from IPython.display import clear_output, display
        from oyFlow.Processing.Utilities import findregexp, findstem, extractFieldsByRegex
        out1 = widgets.Output()
        out2 = widgets.Output()
        out3 = widgets.Output()

        FolderText = Text(value='',placeholder='Enter path to image files',description='Path:',layout=Layout(width='70%', height='30px'))

        GlobText = Text(value='',placeholder='*',description='Regular expression:',layout=Layout(width='70%', height='30px'), style={'description_width': 'initial'})
        GlobText.flag=0
        def select_files(b):
            with out1:
                out1.clear_output()
                root = Tk()
                root.withdraw() # Hide the main window.
                root.call('wm', 'attributes', '.', '-topmost', True) # Raise the root to the top of all windows.
                b.dir = filedialog.askdirectory() # List of selected folder will be set button's file attribute.
                FolderText.value=b.dir+os.path.sep+'*'
                print(b.dir) # Print the list of files selected.

        fileselect = Button(description="Browse directory")
        fileselect.on_click(select_files)

        left_box = HBox([fileselect, FolderText])
        display(left_box)

        def on_value_change_path(change):
            with out1:
                out1.clear_output()
                print(change['new'])
            with out2:
                out2.clear_output()
                GlobText.fnames = fnamesFromPath(change['new'])
                if len(GlobText.fnames):
                    GlobText.value = findregexp(GlobText.fnames)
                    GlobText.globExp = GlobText.value
                    GlobText.patterns = extractFieldsByRegex(GlobText.globExp, GlobText.fnames)
                    GlobText.path = FolderText.value
                    GlobText.flag=1
                else:
                    GlobText.value = 'Empty File List!'
                    GlobText.flag=0
                print(*GlobText.fnames[0:min(5,len(GlobText.fnames))],sep = '\n')

        def fnamesFromPath(pth):
            import glob
            import os
            fnames = glob.glob(pth+'**'+os.path.sep+'*.[fF][cC][sS]',recursive=True)
            fnames.sort(key=os.path.getmtime)
            return fnames

        def on_value_change_glob(change):
            with out3:
                out3.clear_output()
                if len(GlobText.fnames):
                    GlobText.globExp = GlobText.value
                    GlobText.patterns = extractFieldsByRegex(GlobText.globExp, GlobText.fnames)
                else:
                    GlobText.value = 'Empty File List!'

        FolderText.observe(on_value_change_path, names='value')
        FolderText.value=pth

        GlobText.observe(on_value_change_glob, names='value')

        out1.append_stdout(pth)

        print('Files:')
        display(out2)

        display(GlobText)
        return(GlobText)

    def _getMappingFromGUI(self, GlobText):
        from ipywidgets import Button, Text, widgets, HBox, Layout
        from IPython.display import clear_output, display
        import itertools
        out1 = widgets.Output()
        out3 = widgets.Output()
        box1 = HBox()
        buttonDone = Button(description='Done',layout=Layout(width='25%', height='80px'),button_style='success',style=dict(
        font_size='48',
        font_weight='bold'))
        def change_traits(change):
            with out1:
                out1.clear_output()
                ordered_list_of_traits=[]
                for i in range(1, len(box1.children), 2):
                    ordered_list_of_traits.append(box1.children[i].value)
                box1.ordered_list_of_traits = ordered_list_of_traits
                print(ordered_list_of_traits)
                print(*GlobText.patterns[0:min(5,len(GlobText.fnames))],sep = '\n')


        def on_change(t):
            with out3:
                out3.clear_output()
                if GlobText.flag:
                    print('Found regular expression with '+str(len(GlobText.patterns[0]))+' :')

                    parts = GlobText.globExp.split('*')
                    options = ['Channel', 'Position', 'frame', 'Zindex', 'IGNORE']
                    layout = widgets.Layout(width='auto', height='40px') #set width and height

                    dddict = {}
                    for i in range(len(parts)-1):
                        # dynamically create key
                        key = parts[i]
                        # calculate value
                        value = [widgets.Label(parts[i]), widgets.Dropdown(options=options,value=options[i], layout=Layout(width='9%'))]
                        dddict[key] = value
                    key = parts[-1]
                    value = [widgets.Label(parts[-1])]
                    dddict[key] = value

                    ddlist = list(itertools.chain.from_iterable(dddict.values()))
                    box1.children = tuple(ddlist)
                    for dd in box1.children:
                        dd.observe(change_traits,names='value')
                    box1.children[1].value = 'frame'
                    box1.children[1].value = 'Channel'
                    display(box1)


        def on_click_done(b):
            b.disabled=True
            self._load_table_FCSs(box1.ordered_list_of_traits, GlobText.fnames, GlobText.patterns, GlobText.path)



        display(out3)
        if GlobText.flag:
            on_change(GlobText)


        box2 = HBox([out1, buttonDone])
        display(box2)

        GlobText.observe(on_change, names='value')
        buttonDone.on_click(on_click_done)

        return box1


    def _load_table_FCSs(self, ordered_list_of_traits, fnames, patterns, pth):
        """
        Helper function to load minimal metadata from tiffs.
        Returns
        -------
        data_table - pd dataframe of metadata image table
        """
        import pandas as pd
        from os.path import join, isdir


        traitdict={}
        for i,trait in enumerate(ordered_list_of_traits):
            if not trait=='IGNORE':
                key=trait
                value=[p[i] for p in patterns]
                traitdict[key]=value

        usecolslist = []
        data_table = pd.DataFrame(columns=usecolslist)


        data_table['filename']=fnames
        data_table['root_pth'] = data_table.filename

        #Default values
        data_table['acq']=pth
        data_table['XY']=[[0,0]]*len(fnames)
        data_table['Z']=0
        data_table['Zindex']=0
        data_table['Channel']=0
        data_table['Position']=0
        data_table['frame']=0
        data_table['PixelSize']=1


        #update whichever values we got from the names
        for trait in traitdict:
            data_table[trait]=traitdict[trait]

        try:
            data_table['frame']=[int(f) for f in data_table['frame']]
        except:
            pass
                    #todo : convert numerical strings to numbers


        data_table.filename = [join(pth, f.split(os.path.sep)[-1]) for f in data_table.filename]
        self.data_table = data_table
