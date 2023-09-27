# Module for managing flow cytometry data (and other tabular data?)
# Based on FlowCytometryTools https://eyurtsev.github.io/FlowCytometryTools/
# Extended to include workflows, sequential gating, easy grouping, etc.
# As a general strategy, we treat samples as immutable. We assign samples to groups by reference.
# We gate/transform/compensate etc without ever changing the original sample.
# Alon Oyler-Yaniv, OY Lab, Talley Lambert, NIC, HMS, 2022


import logging

from os import path

import numpy as np
from FlowCytometryTools import FCMeasurement
from FlowCytometryTools import PolyGate as _PolyG
from FlowCytometryTools import QuadGate as _QuadG
from FlowCytometryTools import ThresholdGate as _ThreshG
from FlowCytometryTools import IntervalGate as _IntG
import copy
import warnings

warnings.simplefilter(action='ignore', category=UserWarning)


md_logger = logging.getLogger(__name__)
md_logger.setLevel(logging.DEBUG)


class Workspace(object):
    """
    Class for experiment workspace.
    Parameters
    ----------
    pth : str full path to saved Workspace file or to directory with fcs files

    Returns
    -------
    Workspace instance

    Class properties
    ----------------
     'type',
     'samples',
     'fluorescent_channel_names',
     'channel_names',
     'groups',
     'datadir',
     'IDs'
     'LinearRange'
     'compmat'

    Class methods
    -------------
     'load',
     'save',
     'addgroup'
    """

    def __init__(self, pth=None):
        loaded = False
        self._transform = "asinh"

        if not path.isdir(pth):
            try:
                r = Workspace.load(pth)
                self.__dict__.update(r.__dict__)
                loaded = True
            except:
                print(pth + " isn't a workspace, trying to load from fcs files")
        if not loaded:
            #try:
            self.type = "Flow"
            self.samples = ordDict()
            self.groups = ordDict()
            self._LinearRange = 150
            if pth is None:
                return

            if path.isdir(pth):
                self.datadir = pth
            else:
                self.datadir, _ = path.split(pth)

            self._rootgate = rootGate()

            self.samples = self._load_samples()
            self.fluorescent_channel_names = list(
                self.samples[0].fluorescent_channel_names
            )
            self.channel_names = list(self.samples[0].channel_names)
            self._compmat = np.eye(len(self.fluorescent_channel_names))
            self.IDs = list(self.samples)
            self.groups["All Samples"] = self._FlowGroup(IDs=self.IDs)
            loaded = True
            print("\nloaded workspace from fcs files")
            #except:
            #    print(pth + " nothing to load here")

    @property
    def LinearRange(self):
        """
        Linear range to use for transforming data to asinh scale
        """
        return self._LinearRange

    @LinearRange.setter
    def LinearRange(self, value):
        for s in self.samples.values():
            s.LinearRange = dict((el, value) for el in self.fluorescent_channel_names)
        self._LinearRange = value

    @property
    def compmat(self):
        """
        Compensation matrix
        """
        return self._compmat

    @compmat.setter
    def compmat(self, value):
        self._compmat = value

    def _load_samples(self):
        """
        Helper function to load fcs files
        """
        import glob
        import os
        from pathlib import Path
        from natsort import natsorted

        fnames = glob.glob(
            self.datadir + "**" + os.path.sep + "*.[fF][cC][sS]", recursive=True
        )

        fnames = natsorted(fnames)

        for file in fnames:
            self.samples[Path(file).stem] = self._FlowMeasurement(
                ID=Path(file).stem, datafile=file
            )
        return self.samples

    def copy(self):
        import copy

        return copy.deepcopy(self)

    def _FlowGroup(self, **kwargs):
        """
        Wrapper for creating a new group
        """
        return _FCG(self, **kwargs)

    def _FlowMeasurement(self, **kwargs):
        """
        Wrapper for creating a new measurement
        """
        return _FCM(self, **kwargs)

    def addgroup(self, **kwargs):
        """
        Function for adding a new group to Workspace

        Parameters:
        ===========
        name : str or None, [None] - name of group, if not given it is taken as the pattern, if no pattern, regexp of the IDs
        IDs : list, subset of Workspace IDs - specific list of IDs to include in the group. Superceded pattern
        pattern : str or None, [None] - pattern of samples to include, if specific IDs aren't given

        """
        g = self._FlowGroup(**kwargs)
        assert g.samples.keys(), "cannot create empty group!"
        self.groups[g.name] = g

    def __call__(self):
        print("Workspace object for  " + str(self.datadir) + ".")
        print("\nData is of type " + self.type + ".")
        print("\nAvailable channels are : " + ", ".join(list(self.channel_names)) + ".")
        print("\nAvailable groups are : " + ", ".join(list(self.groups)) + ".")

    def save(self, filename=None):
        """
        save workspace
        """
        from oyFlow.Utilities import is_path_exists_or_creatable

        if filename == None:
            if self.filename:
                filename = self.filename
        assert filename != None, "Must provide file name if saving for the first time"
        assert is_path_exists_or_creatable(filename), filename + " is not a valid path"
        self.filename = filename
        import dill

        with open(filename, "wb") as dbfile:
            dill.dump(self, dbfile)
            print("saved workspace.")

    @classmethod
    def load(cls, filename):
        """
        load
        """
        from oyFlow.Utilities import is_path_exists_or_creatable

        assert is_path_exists_or_creatable(filename), filename + " is not a valid path"
        import dill

        with open(filename, "rb") as dbfile:
            r = dill.load(dbfile)
            print("loaded workspace from saved file.")
        return r


class _FCG(object):
    """
    Class for sample group.
    Parameters
    ----------
    name : str or None, [None] - name of group, if not given it is taken as the pattern, if no pattern, regexp of the IDs
    IDs : list, subset of Workspace IDs - specific list of IDs to include in the group. Superceded pattern
    pattern : str or None, [None] - pattern of samples to include, if specific IDs aren't given

    Returns
    -------
    Group instance

    Class properties
    ----------------
     'workspace',
     'gates',
     'fluorescent_channel_names',
     'channel_names',
     'IDs',
     'name',
     'samples',
     'gates',
     'LinearRange',
     'compmat'

    Class methods
    -------------
     'append', add sample to group
     'remove', remove sample from group
     'transform',
     'compensate',
     'append_group_gate', - add a gate and all of its parents to each sample. Gates are associated to a sample and can be modified individually
     'remove_group_gate'
    """

    def __init__(self, outer, name=None, IDs=[], pattern=None, **kwargs):
        from oyFlow.Utilities import findregexp

        self.workspace = outer
        self.gates = ordDict({"root": copy.deepcopy(outer._rootgate)})
        self.fluorescent_channel_names = self.workspace.fluorescent_channel_names
        self.channel_names = self.workspace.fluorescent_channel_names
        self.LinearRange = dict((el, 150) for el in self.fluorescent_channel_names)
        self._compmat=np.eye(len(self.workspace.fluorescent_channel_names))
        self._transform = 'None'

        if name is None:
            if pattern is not None:
                name = pattern
            elif len(IDs):
                name = findregexp(IDs)
            else:
                name = "Empty group"
        self.name = name

        if len(IDs) == 0:
            if pattern != None:
                IDs = [id for id in self.workspace.IDs if pattern in id]
            else:
                IDs = []
        self.IDs = IDs

    def __call__(self):
        print("Group object for  experiment " + str(self.workspace.datadir) + ".")
        print("\nGroup name is " + self.name + ".")
        print(
            "\nAvailable channels are : "
            + ", ".join(list(self.fluorescent_channel_names))
            + "."
        )
        print("\nAvailable IDs are : " + ", ".join(list(self.IDs)) + ".")
        print("\nAvailable group gates are : \n")
        self.print_gates()

    @property
    def samples(self):
        import copy
        from oyFlow.Flow import ordDict

        samples = ordDict()
        for id in self.IDs:
            samples[id] = self.workspace.samples[id]
        return samples

    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, value):
        self._workspace = value

    # @property
    # def compmat(self):
    #     return self._compmat

    # @property
    # def _transform(self):
    #     return self._workspace._transform

    def append(self, ID):
        self.IDs.append(ID)

    def remove(self, ID):
        self.IDs.remove(ID)

    def append_group_gate(self, gate):
        from copy import deepcopy

        for s in self.samples.values():
            s.append_gate(gate)
        self.gates = deepcopy(self.samples[0].gates)

    def remove_group_gate(self, gate):
        from copy import deepcopy

        for s in self.samples.values():
            s.remove_gate(gate)
        self.gates = deepcopy(self.samples[0].gates)

    def print_gates(self):
        from anytree import RenderTree
        from anytree.render import ContRoundStyle

        for pre, fill, node in RenderTree(self.gates["root"], style=ContRoundStyle()):
            if node.name in self.gates.keys():
                treestr = "%s%s" % (pre, node.name)
                print(treestr.ljust(8), node.channels)

    def apply_gate(self, gate, apply_parents=True):
        from copy import deepcopy

        tmp = _FCG_nosamps(
            self
        )  # this is a weird hack so that this function returns a group
        tmpsamps = deepcopy(self.samples)
        for name, s in tmpsamps.items():
            tmpsamps[name] = tmpsamps[name].apply_gate(
                gate, apply_parents=apply_parents
            )
        tmp.samples = tmpsamps
        return tmp

    def transform(self, transform=None, channels=None, **kwargs):
        """
        Apply transform to all samples.
        Parameters
        ==========
        transform=[*asinh* with self.LinearRange| 'hlog' | 'tlog' | 'glog' | 'None' | callable]
        direction='forward',
        channels : By default, all fluorescent channels
        return_all=True,
        auto_range=True,
        use_spln=True,
        get_transformer=False,
        ID=None,
        apply_now=True,
        args=(),
        **kwargs,

            returns a group with transformed samples

        """

        from copy import deepcopy

        if transform is None:
            transform = self.workspace._transform
        
        tmp = _FCG_nosamps(
            self
        )  # this is a weird hack so that this function returns a group
        tmpsamps = deepcopy(self.samples)
        for name, s in tmpsamps.items():
            tmpsamps[name] = tmpsamps[name].transform(
                transform=transform, channels=channels, **kwargs
            )
        tmp.samples = tmpsamps
        tmp._transform = transform
        return tmp

    def compensate(self, compmat=None):
        """
        Apply compensation to all sample.
        Parameters
        ==========
        compmat - compensation matrix

        returns a group with compensated samples
        """
        if compmat is None:
            compmat = self.workspace._compmat
        if np.any(compmat == 'None'):
            compmat = np.eye(len(self.fluorescent_channel_names))

        from copy import deepcopy

        tmp = _FCG_nosamps(
            self
        )  # this is a weird hack so that this function returns a group
        tmpsamps = deepcopy(self.samples)
        for name, s in tmpsamps.items():
            tmpsamps[name] = tmpsamps[name].compensate(compmat=compmat)
        tmp.samples = tmpsamps
        tmp._compmat = compmat
        return tmp

    @property
    def mean(self, ChanList=None):
        if ChanList is None:
            ChanList = self.fluorescent_channel_names
        import pandas as pd

        things = self.samples[0][ChanList].mean()
        for id in self.IDs[1:]:
            things = pd.concat(
                [things, self.samples[id][ChanList].mean()], axis=1, ignore_index=True
            )
        things.columns = self.IDs
        return things.T

    @property
    def mode(self, ChanList=None):
        if ChanList is None:
            ChanList = self.fluorescent_channel_names
        import pandas as pd

        things = self.samples[0][ChanList].mode().iloc[[0]].squeeze()
        for id in self.IDs[1:]:
            things = pd.concat(
                [things, self.samples[id][ChanList].mode().iloc[[0]].squeeze()],
                axis=1,
                ignore_index=True,
            )
        things.columns = self.IDs
        return things.T

    @property
    def std(self, ChanList=None):
        if ChanList is None:
            ChanList = self.fluorescent_channel_names
        import pandas as pd

        things = self.samples[0][ChanList].std()
        for id in self.IDs[1:]:
            things = pd.concat(
                [things, self.samples[id][ChanList].std()], axis=1, ignore_index=True
            )
        things.columns = self.IDs
        return things.T

    @property
    def median(self, ChanList=None):
        if ChanList is None:
            ChanList = self.fluorescent_channel_names
        import pandas as pd

        things = self.samples[0][ChanList].median()
        for id in self.IDs[1:]:
            things = pd.concat(
                [things, self.samples[id][ChanList].median()], axis=1, ignore_index=True
            )
        things.columns = self.IDs
        return things.T

    @property
    def count(self, ChanList=None):
        if ChanList is None:
            ChanList = self.fluorescent_channel_names
        import pandas as pd

        things = self.samples[0][ChanList].count()
        for id in self.IDs[1:]:
            things = pd.concat(
                [things, self.samples[id][ChanList].count()], axis=1, ignore_index=True
            )
        things.columns = self.IDs
        return things.T


class _FCG_nosamps(_FCG):
    def __init__(self, other):
        self.__dict__.update(other.__dict__)
        self._samples = ordDict()

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, value):
        self._samples = value


class _FCM(FCMeasurement):
    """
    Class for sample.
    Parameters
    ----------
    workspace :
    LinearRange
    ID :
    datafile=None,
    readdata=False,
    readdata_kwargs={},
    metafile=None,
    readmeta=True,
    readmeta_kwargs={},

    Returns
    -------
    Group instance

    Class properties
    ----------------
     'workspace',
     'gates',
     'fluorescent_channel_names',
     'channel_names',
     'ID',
     'gates',
     'LinearRange',
     'compmat'


    Class methods
    -------------
     'transform',
     'compensate',
     'append_gate', - add a gate and all of its parents to sample.
     'remove_gate'
     'apply_gate' - apply gate and parents
     'print_gates'
    """

    def __init__(self, workspace, LinearRange=150, **kwargs):
        FCMeasurement.__init__(self, **kwargs)

        blankInd = []
        if "FSC-A" in self.channel_names:
            blankInd.append(self.channel_names.index("FSC-A"))
        if "FSC-H" in self.channel_names:
            blankInd.append(self.channel_names.index("FSC-H"))
        if "FSC-W" in self.channel_names:
            blankInd.append(self.channel_names.index("FSC-W"))
        if "SSC-A" in self.channel_names:
            blankInd.append(self.channel_names.index("SSC-A"))
        if "SSC-H" in self.channel_names:
            blankInd.append(self.channel_names.index("SSC-H"))
        if "SSC-W" in self.channel_names:
            blankInd.append(self.channel_names.index("SSC-W"))
        if "Time" in self.channel_names:
            blankInd.append(self.channel_names.index("Time"))

        ChanInd = list(set(np.arange(0, len(self.channel_names))) - set(blankInd))
        ChanList = list(np.array(self.channel_names)[ChanInd])
        self.fluorescent_channel_names = list(ChanList)
        self.LinearRange = dict((el, LinearRange) for el in ChanList)
        self.workspace = workspace
        self.gates = ordDict({"root": copy.deepcopy(workspace._rootgate)})
        self._compmat = np.eye(len(self.fluorescent_channel_names))
        self._transform = 'None'

    def __call__(self):
        print("Sample object for  experiment " + str(self.workspace.datadir) + ".")
        print("\nfor sample ID : " + self.ID + ".")
        print(
            "\nAvailable channels are : "
            + ", ".join(list(self.fluorescent_channel_names))
            + "."
        )
        print("\nAvailable sample gates are : \n")
        self.print_gates()

    @property
    def workspace(self):
        return self._workspace

    @workspace.setter
    def workspace(self, value):
        self._workspace = value

    def copy(self):
        import copy

        return copy.deepcopy(self)

    @property
    def raw_data(self):
        return super().data

    def subsample(self,N=15000):
        tmpsample = self.copy()
        N = min(len(tmpsample.data),N)
        new_data = tmpsample.data.sample(N)
        tmpsample.data = new_data
        return tmpsample

    def compensate(self, compmat=None):
        """
        Apply compensation matrix to sample. If not provided, the self.compmat will be applied
        Parameters
        ==========
        compmat - [self.compmat]
        """
        if compmat is None:
            compmat = self.workspace._compmat
        if np.any(compmat == 'None'):
            compmat = np.eye(len(self.fluorescent_channel_names))
        chan = self.fluorescent_channel_names
        tmpsample = self.copy()
        new_data = tmpsample.data
        new_data[chan] = new_data[chan].dot(compmat)
        tmpsample.data = new_data
        tmpsample._compmat = compmat
        return tmpsample

    def transform(self, transform=None, channels=None, **kwargs):
        """
        Apply transform to sample.
        Parameters
        ==========
        transform=[*asinh* with self.LinearRange| 'hlog' | 'tlog' | 'glog' | callable]
        direction='forward',
        channels : By default, all fluorescent channels
        return_all=True,
        auto_range=True,
        use_spln=True,
        get_transformer=False,
        ID=None,
        apply_now=True,
        args=(),
        **kwargs,


        """
        if transform is None:
            transform = self.workspace._transform
            
        if channels is None:
            channels = self.fluorescent_channel_names
        if transform == 'None':
            tmpsample = self.copy()
        else:
            if transform != "asinh":
                tmpsample = FCMeasurement.transform(
                    self, transform=transform, channels=channels, **kwargs
                )
            elif transform == "asinh":
                tmpsample = self.copy()
                new_data = tmpsample.data
                for chan in channels:
                    new_data[chan] = np.arcsinh(new_data[chan] / self.LinearRange[chan])
                tmpsample.data = new_data
        tmpsample._transform = transform
        return tmpsample

    def apply_gate(self, gate=None, apply_parents=True, verbose=True):
        """
        Apply gate to sample. Applies the gate and all of it's parents(?).
        Parameters
        ==========
        gate : str
        apply_parents : True By default

        """
        tmpsample = self.copy()
        if gate == None or gate == "root":
            if verbose:
                print("No gate supplied, didn't do anything")
            return tmpsample
        assert isinstance(gate, str), "gate must be a string"
        assert gate in self.gates.keys(), "gate must be appended to sample"
        gate = self.gates[gate]
        if apply_parents:
            if isinstance(gate.parent, (_PolyG, _QuadG, _ThreshG, _IntG)):
                tmpsample = tmpsample.apply_gate(gate=gate.parent.name)

        assert tmpsample._transform == gate._transform, 'can\'t apply gate to a sample with a different transform. try to add .transform()'
        assert np.all(tmpsample._compmat == gate._compmat), 'can\'t apply gate to a sample with a different compensation matrix. Try to add .compensate()'
        tmpsample = tmpsample.gate(gate)

        return tmpsample

    def append_gate(self, gate=None, update=False):
        """
        append gate to sample. Parent must be already appended to the sample
        Parameters
        ==========
        gate : str
        """
        import copy

        gate = copy.deepcopy(gate)
        if gate is None:
            print("No gate supplied, didn't do anything")
            return
        assert not isinstance(gate, rootGate), "Can not have multiple roots"

        if gate.name in list(self.gates):
            print("gate already exists in sample, skipping!")
            return

        self.gates[gate.name] = gate
        if gate.parent.name not in list(self.gates):
            self.append_gate(gate.parent)

        gate.children = ()
        gate.parent = self.gates[gate.parent.name]
        return gate

    def remove_gate(self, gate=None):
        """
        remove gate from sample. Removes gate and all of it's children.
        Parameters
        ==========
        gate : str
        """
        if gate is None:
            print("No gate supplied, didn't do anything")
            return
        assert isinstance(gate, str), "gate must be a string"
        assert gate in self.gates.keys(), "gate must be appended to sample"
        assert gate != "root", "Can remove root gate"
        gate = self.gates[gate]
        for c in gate.children:
            if isinstance(c, (_PolyG, _QuadG, _ThreshG, _IntG)):
                if c.name in self.gates.keys():
                    self.remove_gate(gate=c.name)
        if gate.name in self.gates.keys():
            self.gates[gate.name].parent = None
            self.gates[gate.name].children = ()
            del self.gates[gate.name]

    def print_gates(self):
        """
        print tree of available gates
        """
        from anytree import RenderTree
        from anytree.render import ContRoundStyle

        for pre, fill, node in RenderTree(self.gates["root"], style=ContRoundStyle()):
            if node.name in self.gates.keys():
                treestr = "%s%s" % (pre, node.name)
                print(treestr.ljust(8), node.channels)


class ordDict(dict):
    """
    ordered dictionary class, modifies dictionary by adding retrieval by index
    """

    def __getitem__(self, itm):
        if isinstance(itm, int):
            return list(self.values())[itm]
        else:
            return super().__getitem__(itm)


## Adding parentage to gates so that a gating tree can be constructed.

from anytree import NodeMixin


class rootGate(NodeMixin):
    def __init__(self):
        self.parent = None
        self.name = "root"
        self.channels = ""


class PolyGate(_PolyG, NodeMixin):
    def __init__(self, vert, channels, parent=None, workspace=None, **kwargs):
        _PolyG.__init__(self, vert, channels, **kwargs)
        self.parent = parent
        self._compmat = workspace._compmat
        self._transform = workspace._transform

    def selector(self, ax, alpha=0.2, color="tab:blue"):
        from matplotlib.widgets import PolygonSelector
        from oyFlow.Gating import inverse_transform, transform

        def _onselectpoly(verts):
            self.vert = verts

        span = PolygonSelector(
            ax,
            _onselectpoly,
            useblit=True,
            props=dict(alpha=alpha, color=color),
        )

        span._selection_completed = True
        verts_tmp = self.vert

        span.verts = verts_tmp

        return span


class QuadGate(_QuadG, NodeMixin):
    def __init__(self, vert, channels, parent=None, workspace=None, **kwargs):
        _QuadG.__init__(self, vert, channels, **kwargs)
        self.parent = parent
        self._compmat = workspace._compmat
        self._transform = workspace._transform


class ThresholdGate(_ThreshG, NodeMixin):
    def __init__(self, vert, channels, parent=None, workspace=None, **kwargs):
        _ThreshG.__init__(self, vert, channels, **kwargs)
        self.parent = parent
        self._compmat = workspace._compmat
        self._transform = workspace._transform


class IntervalGate(_IntG, NodeMixin):
    def __init__(self, vert, channels, parent=None, workspace=None, **kwargs):
        _IntG.__init__(self, vert, channels, **kwargs)
        self.parent = parent
        self._compmat = workspace._compmat
        self._transform = workspace._transform

    def selector(self, ax, alpha=0.2, color="tab:blue"):
        from matplotlib.widgets import SpanSelector
        from oyFlow.Gating import inverse_transform, transform

        def _onselectspan(xmin, xmax):
            self.vert = (xmin, xmax)

        span = SpanSelector(
            ax,
            _onselectspan,
            "horizontal",
            useblit=True,
            props=dict(alpha=alpha, facecolor=color, edgecolor=color),
            interactive=True,
            drag_from_anywhere=False,
            ignore_event_outside=True,
        )
        span._selection_completed = True

        span.extents = self.vert

        return span
