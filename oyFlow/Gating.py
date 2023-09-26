from oyFlow.Flow import Workspace
from PyQt5.QtWidgets import QTreeWidgetItem
import numpy as np


def _get_ticks(W):
    """
    Helper function for calculating position of tickmarks
    W is a workspace that has a W._transform attribute
    """
    a1 = (
        np.expand_dims(np.arange(1, 10), 0) * 10.0 ** np.expand_dims(np.arange(1, 6), 1)
    ).flatten()
    a2 = -np.flip(
        (
            np.expand_dims(np.arange(1, 10), 0)
            * 10.0 ** np.expand_dims(np.arange(1, 2), 1)
        ).flatten()
    )
    a = np.insert(a1, 0, 0)
    a = np.hstack((a2, a))
    a = np.append(a, 10**6)
    emptyStr = [""] * 8
    import itertools

    XTickLabel = list(
        itertools.chain.from_iterable(
            [
                emptyStr,
                [""],
                ["0"],
                [""],
                emptyStr,
                [""],
                emptyStr,
                ["1000"],
                emptyStr,
                ["10000"],
                emptyStr,
                ["100000"],
                emptyStr,
                ["1000000"],
            ]
        )
    )
    XTicks = transform(a, W)
    return XTicks, XTickLabel


def inverse_transform(x, W):
    """
    Helper function for simple inverse transform (back to linear scale)
    W is a workspace that has a W._transform attribute
    """
    if W._transform == None:
        return x
    else:
        if W._transform == "asinh":
            return W._LinearRange * np.sinh(x)
        else:
            from FlowCytometryTools.core.transforms import Transformation

            transformer = Transformation(W._transform, "inverse")
            return transformer(x)


def transform(x, W):
    """
    Helper function for simple transform (to transformed space)
    W is a workspace that has a W._transform attribute
    """
    if W._transform == None:
        return x
    else:
        if W._transform == "asinh":
            return np.arcsinh(x / W._LinearRange)
        else:
            from FlowCytometryTools.core.transforms import Transformation

            transformer = Transformation(W._transform, "forward")
            return transformer(x)


def _prompt_for_gate_name(gate_name, W, widget):
    """
    helper widget for naming gates
    """
    from magicgui.widgets import request_values

    name_wid = request_values(
        name={
            "annotation": str,
            "value": gate_name,
            "label": "Enter unique gate name:",
        },
    )
    if name_wid:
        gate_name = name_wid["name"]
        if gate_name in [
            ch.name
            for ch in W.groups[widget.group.value]
            .samples[widget.sample.value]
            .gates[widget.population.value]
            .children
        ]:
            gate_name = _prompt_for_gate_name(gate_name, W, widget)
        return gate_name
    else:
        return None


def add_tree_to_qtreewidget(qtree_widget, anytree_node, parent_item=None):
    """
    Helper function for populating QTreeWidget with anytree tree
    """
    # Create a QTreeWidgetItem for the current anytree node
    item = QTreeWidgetItem(parent_item if parent_item else qtree_widget)
    item.setText(0, anytree_node.name)  # Set the text for the item

    # Recursively add child nodes
    for child in anytree_node.children:
        add_tree_to_qtreewidget(qtree_widget, child, item)

    # Expand the item to show its children
    item.setExpanded(True)

    return item


def set_item_by_name(qtree_widget, item_name="root", parent_item=None):
    """
    Helper function for setting a chosen item by it's name only in a QTreeWidget
    """
    for item_index in (
        range(parent_item.childCount())
        if parent_item
        else range(qtree_widget.topLevelItemCount())
    ):
        item = (
            parent_item.child(item_index)
            if parent_item
            else qtree_widget.topLevelItem(item_index)
        )
        if (
            item.text(0) == item_name
        ):  # Assuming you want to match against the text in the first column
            qtree_widget.setCurrentItem(item)
            return True  # Exit the function once the item is found and selected
        # Recursively search in child items
        set_item_by_name(qtree_widget, item_name, item)


def gater(W):
    """
    gater app
    W is a workspace
    """
    from typing import List
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from magicgui import magicgui
    from magicgui.widgets import Checkbox, Container, Label
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
    )
    from natsort import natsorted
    from scipy.stats import gaussian_kde
    import qtawesome as qta
    from matplotlib.widgets import SpanSelector, PolygonSelector, RectangleSelector
    from oyFlow.Flow import IntervalGate, PolyGate

    from itertools import cycle
    import matplotlib.style as mplstyle

    mplstyle.use("fast")

    cmaps = cycle(["green", "blue", "red", "cyan", "magenta", "yellow"])

    # make axis where things will be drawn
    matplotlib.use("Agg")
    mpl_fig, ax = plt.subplots()
    fc = FigureCanvasQTAgg(mpl_fig)
    ax.set_position([0.18, 0.18, 0.8, 0.8])

    # callbacks for initial gate selections
    def _onselectspan(xmin, xmax):
        """
        Callback for when a span is selected using the SpanSelector widget
        """
        gate_name = (
            "I_"
            + widget.channel_1.value
            + "_"
            + widget.sample.value
            + "_"
            + widget.group.value
            + "_"
            + widget.population.value
        )
        gate_name = _prompt_for_gate_name(gate_name, W, widget)
        if not gate_name:
            return 0

        # make gate
        gate = IntervalGate(
            (xmin, xmax),
            channels=widget.channel_1.value,
            region="in",
            name=gate_name,
            parent=W.groups[widget.group.value]
            .samples[widget.sample.value]
            .gates[widget.population.value],
            workspace=W,
        )

        # append to sample
        W.groups[widget.group.value].samples[widget.sample.value].append_gate(
            gate=gate, update=True
        )

        gate.parent = None
        widget._active_gate.set_active(False)
        widget._active_gate.set_visible(False)
        widget._active_gate.clear()
        _on_sample_changed()

    def _onspanbutton():
        """
        Callback for when span gate button is pushed, created a SpanSelector
        """
        span = SpanSelector(
            ax,
            _onselectspan,
            "horizontal",
            useblit=False,
            props=dict(alpha=0.5, facecolor="tab:blue"),
            interactive=True,
            drag_from_anywhere=False,
            ignore_event_outside=True,
        )
        widget._active_gate = span

    def _onselectpoly(verts):
        """
        Callback for when a span is selected using the SpanSelector widget
        """
        gate_name = (
            "P_"
            + widget.channel_1.value
            + "_"
            + widget.channel_2.value
            + "_"
            + widget.sample.value
            + "_"
            + widget.group.value
            + "_"
            + widget.population.value
        )
        gate_name = _prompt_for_gate_name(gate_name, W, widget)
        if not gate_name:
            return 0

        # make gate
        gate = PolyGate(
            verts,
            channels=(widget.channel_1.value, widget.channel_2.value),
            region="in",
            name=gate_name,
            parent=W.groups[widget.group.value]
            .samples[widget.sample.value]
            .gates[widget.population.value],
            workspace=W,
        )

        # append to sample
        W.groups[widget.group.value].samples[widget.sample.value].append_gate(
            gate=gate, update=True
        )

        gate.parent = None
        widget._active_gate.set_active(False)
        widget._active_gate.set_visible(False)
        widget._active_gate.clear()
        _on_sample_changed()

    def _onpolybutton():
        """
        Callback for when span gate button is pushed, created a SpanSelector
        """
        span = PolygonSelector(
            ax,
            _onselectpoly,
            useblit=True,
            props=dict(alpha=0.5, color="tab:blue"),
        )
        widget._active_gate = span

    def _onselectrect(eclick, erelease):
        """
        Callback for when a span is selected using the SpanSelector widget
        """

        verts = [
            (eclick.xdata, eclick.ydata),
            (eclick.xdata, erelease.ydata),
            (erelease.xdata, erelease.ydata),
            (erelease.xdata, eclick.ydata),
        ]
        gate_name = (
            "P_"
            + widget.channel_1.value
            + "_"
            + widget.channel_2.value
            + "_"
            + widget.sample.value
            + "_"
            + widget.group.value
            + "_"
            + widget.population.value
        )
        gate_name = _prompt_for_gate_name(gate_name, W, widget)
        if not gate_name:
            return 0

        # make gate
        gate = PolyGate(
            verts,
            channels=(widget.channel_1.value, widget.channel_2.value),
            region="in",
            name=gate_name,
            parent=W.groups[widget.group.value]
            .samples[widget.sample.value]
            .gates[widget.population.value],
            workspace=W,
        )

        # append to sample
        W.groups[widget.group.value].samples[widget.sample.value].append_gate(
            gate=gate, update=True
        )

        gate.parent = None
        widget._active_gate.set_active(False)
        widget._active_gate.set_visible(False)
        widget._active_gate.clear()
        _on_sample_changed()

    def _onrectbutton():
        """
        Callback for when span gate button is pushed, created a SpanSelector
        """
        span = RectangleSelector(
            ax,
            _onselectrect,
            useblit=True,
            interactive=True,
            props=dict(alpha=0.5, color="tab:blue"),
        )
        widget._active_gate = span

    def _onremovesamplegate():
        """
        Callback for when clicking the remove sample gate button
        """
        if widget.pop_tree_widget.currentItem():
            current_item_name = widget.pop_tree_widget.currentItem().text(0)
            assert current_item_name != "root", "Can't delete root!"
            W.groups[widget.group.value].samples[widget.sample.value].remove_gate(
                current_item_name
            )
            _on_sample_changed()

    def _onaddgroupgate():
        """
        Callback for when clicking the add group gate button
        """
        if widget.pop_tree_widget.currentItem():
            current_item_name = widget.pop_tree_widget.currentItem().text(0)
            g = W.groups[widget.group.value]
            s = W.groups[widget.group.value].samples[widget.sample.value]
            g.append_group_gate(s.gates[current_item_name])
            _on_sample_changed()

    def _onremovegroupgate():
        """
        Callback for when clicking the remove group gate button
        """
        if widget.pop_tree_widget.currentItem():
            current_item_name = widget.pop_tree_widget.currentItem().text(0)
            assert current_item_name in list(
                W.groups[widget.group.value].gates
            ), "Not a group gate! You can delete it from the sample."
            assert current_item_name != "root", "Can't delete root!"
            W.groups[widget.group.value].remove_group_gate(current_item_name)
            _on_sample_changed()

    ### Define main widget
    @magicgui(
        auto_call=False,
        group={"choices": natsorted([a for a in W.groups.keys()])},
        sample={
            "choices": natsorted([a for a in W.groups["All Samples"].samples.keys()])
        },
        channel_1={"choices": list(W.channel_names)},
        channel_2={"choices": list(W.channel_names)},
        population={"choices": list(W.groups[0].samples[0].gates)},
    )
    def widget(
        group: List[str],
        sample: List[str],
        channel_1: List[str],
        channel_2: List[str],
        population: str,
    ):
        """Main widget for gating"""
        # get curent sample
        sample_d = (
            W.groups[group]
            .samples[sample]
            .compensate(W._compmat)
            .transform(W._transform)
            .apply_gate(population, verbose=False)
            .data
        )

        # clear whatever gates might be showing:
        [g.set_active(False) for g in widget._gates]
        [g.set_visible(False) for g in widget._gates]
        [g.clear() for g in widget._gates]
        widget._gates = []

        # clear current figure and get channel choices
        ax.cla()
        ch1_choice = channel_1
        ch2_choice = channel_2
        # if 2 different channels, make a scatter plot
        if ch1_choice != ch2_choice:
            x = sample_d[ch1_choice]
            y = sample_d[ch2_choice]

            


            # Evaluate density at each data point
            if len(x) > 20:
                xy = np.vstack([x, y])
                kde = gaussian_kde(xy[:, : np.min((250, len(x)))])
                density = kde(xy)
                normalized_density = (density - density.min()) / (
                    density.max() - density.min()
                )
            else:
                normalized_density = "k"

            # Create scatter plot with colored points based on local density
            ax.scatter(x, y, c=normalized_density, cmap="viridis", s=5)
            ax.set_xlabel(ch1_choice)
            ax.set_ylabel(ch2_choice)

            ax.set_xlim([-500, 250000])
            ax.set_ylim([-500, 250000])
            # fix axis labels
            if W._transform:
                if ch1_choice in W.fluorescent_channel_names:
                    ax.set_xlim(transform(np.array([-500, 250000]), W))
                    XTicks, XTickLabel = _get_ticks(W)
                    ax.set_xticks(XTicks)
                    ax.set_xticklabels(XTickLabel)
                if ch2_choice in W.fluorescent_channel_names:
                    ax.set_ylim(transform(np.array([-500, 250000]), W))
                    XTicks, XTickLabel = _get_ticks(W)
                    ax.set_yticks(XTicks)
                    ax.set_yticklabels(XTickLabel)

            # put buttons in place. polygon and quad for scatter
            [act.setVisible("scat" in key) for key, act in widget.added_actions.items()]

            # put relevant gates in place
            child_gates = (
                W.groups[widget.group.value]
                .samples[widget.sample.value]
                .gates[widget.population.value]
                .children
            )
            # put gates in place if they have just one channel (span or threshold) and it is the current shown channel.
            widget._gates = [
                ch.selector(ax, color=next(cmaps))
                for ch in child_gates
                if len(ch.channels) == 2
                and ch.channels[0] == ch1_choice
                and ch.channels[1] == ch2_choice
            ]

        # if ch1 and ch2 are the same make a histogram
        elif ch1_choice == ch2_choice:
            x = sample_d[ch1_choice]
            ax.hist(x, 100)
            ax.set_xlabel(ch1_choice)

            ax.set_xlim([-500, 250000])
            if W._transform:
                if ch1_choice in W.fluorescent_channel_names:
                    ax.set_xlim(transform(np.array([-500, 250000]), W))
                    XTicks, XTickLabel = _get_ticks(W)
                    ax.set_xticks(XTicks)
                    ax.set_xticklabels(XTickLabel)

            # put buttons in place. span and thresh for hist
            [act.setVisible("hist" in key) for key, act in widget.added_actions.items()]

            # put relevant gates in place
            child_gates = (
                W.groups[widget.group.value]
                .samples[widget.sample.value]
                .gates[widget.population.value]
                .children
            )
            # put gates in place if they have just one channel (span or threshold) and it is the current shown channel.
            widget._gates = [
                ch.selector(ax, color=next(cmaps))
                for ch in child_gates
                if len(ch.channels) == 1 and ch.channels[0] == ch1_choice
            ]

        fc.draw()

    ### make toolbar and special buttons for selecting gates
    toolbar = NavigationToolbar(fc, widget.parent)
    widget.added_actions = {}

    # span gate button
    act_range = toolbar.addAction(
        qta.icon("mdi.arrow-expand-horizontal"), "add span gate"
    )
    act_range.triggered.connect(_onspanbutton)
    widget.added_actions["hist_act_range"] = act_range

    # polygon gate button
    act_poly = toolbar.addAction(qta.icon("mdi.shape-polygon-plus"), "add polygon gate")
    act_poly.triggered.connect(_onpolybutton)  ##TODO
    widget.added_actions["scat_act_poly"] = act_poly

    # rect gate button
    act_quad = toolbar.addAction(
        qta.icon("mdi.shape-rectangle-plus"), "add rectangle gate"
    )
    act_quad.triggered.connect(_onrectbutton)  ##TODO
    widget.added_actions["scat_act_quad"] = act_quad

    ### Make population tree widget showing the different gates in a nice hirarchy
    from PyQt5.QtWidgets import QVBoxLayout, QTreeWidget, QToolBar, QAction

    widget.pop_tree_widget = QTreeWidget()
    widget.pop_tree_widget.setHeaderLabel("population")
    # make the native population choice invisible
    widget.population.visible = False

    # make toolbar for population tree widget!
    pop_toolbar = QToolBar(widget.pop_tree_widget)
    delete_action = QAction("Delete Gate", pop_toolbar)
    delete_action.triggered.connect(_onremovesamplegate)
    delete_action.setIcon(qta.icon("mdi.delete-forever-outline"))
    pop_toolbar.addAction(delete_action)
    add_to_group_action = QAction("Add Gate to Group", pop_toolbar)
    add_to_group_action.triggered.connect(_onaddgroupgate)
    add_to_group_action.setIcon(qta.icon("mdi6.table-plus"))
    pop_toolbar.addAction(add_to_group_action)
    remove_from_group_action = QAction("Remove Gate from Group", pop_toolbar)
    remove_from_group_action.triggered.connect(_onremovegroupgate)
    remove_from_group_action.setIcon(qta.icon("mdi6.table-minus"))
    pop_toolbar.addAction(remove_from_group_action)

    # callback for population change on tree widget
    def _on_item_select_change():
        selected_items = widget.pop_tree_widget.selectedItems()
        if selected_items:
            widget.population.value = selected_items[0].text(0)
        else:
            set_item_by_name(widget.pop_tree_widget)

    # connect callback
    widget.pop_tree_widget.itemSelectionChanged.connect(_on_item_select_change)

    # Callbacks for channel changes:
    def _on_ch_changed():
        widget()

    widget.channel_1.changed.connect(_on_ch_changed)
    widget.channel_2.changed.connect(_on_ch_changed)

    # Callbacks for sample change:
    @widget.sample.changed.connect
    def _on_sample_changed():
        widget.population.choices = list(
            W.groups[widget.group.value].samples[widget.sample.value].gates.keys()
        )
        popname = widget.population.value
        widget.pop_tree_widget.clear()
        add_tree_to_qtreewidget(
            widget.pop_tree_widget,
            W.groups[widget.group.value].samples[widget.sample.value].gates[0],
        )
        set_item_by_name(widget.pop_tree_widget, item_name=popname)
        widget()

    # Callbacks for group change:
    @widget.group.changed.connect
    def _on_group_changed():
        widget.sample.choices = list(W.groups[widget.group.value].samples.keys())
        widget.sample.value = widget.sample.choices[0]

    # Callbacks for population change:
    @widget.population.changed.connect
    def _on_pop_changed():
        widget()

    ### layout
    container = Container(layout="horizontal")
    layout = container.native.layout()

    l = QVBoxLayout()
    l.addWidget(toolbar)
    l.addWidget(fc)
    widget.root_native_widget.layout().addWidget(widget.pop_tree_widget)
    widget.root_native_widget.layout().addWidget(pop_toolbar)
    layout.addWidget(widget.native)  # adding native, because we're in Qt
    layout.addLayout(l)

    container.show()
    matplotlib.use("Qt5Agg")
    widget.call_button.visible = False
    widget.max_width = 250

    ##inits
    widget._gates = (
        []
    )  # _gates are the active gates currently plotted. initialized to an empty list.

    widget()

    add_tree_to_qtreewidget(
        widget.pop_tree_widget,
        W.groups[widget.group.value].samples[widget.sample.value].gates[0],
    )
    set_item_by_name(widget.pop_tree_widget)

    # exposes, should be removed eventually
    container.pop_tree_widget = widget.pop_tree_widget
    container.widget = widget
    container.ax = ax

    return container
