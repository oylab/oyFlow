from oyFlow.Flow import Workspace
import numpy as np


def _get_ticks(LinearRange):
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
    XTicks = np.arcsinh(a / LinearRange)
    return XTicks, XTickLabel


def inverse_transform(x, W):
    if W._transform == None:
        return x
    else:
        if W._transform == "asinh":
            return W._LinearRange * np.sinh(x)
        else:
            from FlowCytometryTools.core.transforms import Transformation

            transformer = Transformation(W._transform, "inverse")
            return transformer(x)


def _prompt_for_gate_name(gate_name, W, widget):
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
    from PyQt5.QtWidgets import QTreeWidgetItem
    # Create a QTreeWidgetItem for the current anytree node
    item = QTreeWidgetItem(parent_item if parent_item else qtree_widget)
    item.setText(0, anytree_node.name)  # Set the text for the item

    # Recursively add child nodes
    for child in anytree_node.children:
        add_tree_to_qtreewidget(qtree_widget, child, item)

    # Expand the item to show its children
    item.setExpanded(True)

    return item

def set_item_by_name(qtree_widget, item_name='root', parent_item=None):
    for item_index in range(parent_item.childCount()) if parent_item else range(qtree_widget.topLevelItemCount()):
        item = parent_item.child(item_index) if parent_item else qtree_widget.topLevelItem(item_index)
        if item.text(0) == item_name:  # Assuming you want to match against the text in the first column
            qtree_widget.setCurrentItem(item)
            return True  # Exit the function once the item is found and selected
        # Recursively search in child items
        set_item_by_name(qtree_widget, item_name, item)


def gater(W):
    """
    gater app
    """

    from typing import List
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from magicgui import magicgui
    from magicgui.widgets import Checkbox, Container, PushButton, Label
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg,
        NavigationToolbar2QT as NavigationToolbar,
    )
    from natsort import natsorted
    from scipy.stats import gaussian_kde
    import qtawesome as qta
    from matplotlib.widgets import SpanSelector
    from oyFlow.Flow import IntervalGate

    cmaps = ["cyan", "magenta", "yellow", "red", "green", "blue"]

    # make axis where things will be drawn
    matplotlib.use("Agg")
    mpl_fig, ax = plt.subplots()
    fc = FigureCanvasQTAgg(mpl_fig)
    ax.set_position([0.18, 0.18, 0.8, 0.8])


    # callbacks for selections
    def _onselectspan(xmin, xmax):
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

        if widget.channel_1.value in W.fluorescent_channel_names:
            (xmin, xmax) = inverse_transform((xmin, xmax), W)

        gate = IntervalGate(
            (xmin, xmax),
            channels=widget.channel_1.value,
            region="in",
            name=gate_name,
            parent=W.groups[widget.group.value]
            .samples[widget.sample.value]
            .gates[widget.population.value],
        )

        W.groups[widget.group.value].samples[widget.sample.value].append_gate(
            gate=gate, update=False
        )

        #_update_gates_label()
        _on_sample_changed()

    def _onspanbutton():
        span = SpanSelector(
            ax,
            _onselectspan,
            "horizontal",
            useblit=True,
            props=dict(alpha=0.5, facecolor="tab:blue"),
            interactive=True,
            drag_from_anywhere=True,
        )
        widget._gates = span

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
        sample_d = (
            W.groups[group]
            .samples[sample]
            .apply_gate(population, verbose=False)
            .tc_data
        )

        ax.cla()
        ch1_choice = channel_1
        ch2_choice = channel_2
        if ch1_choice != ch2_choice:
            x = sample_d[ch1_choice]
            y = sample_d[ch2_choice]

            # Evaluate density at each data point
            xy = np.vstack([x, y])
            kde = gaussian_kde(xy[:, : np.min((250, len(x)))])
            density = kde(xy)
            normalized_density = (density - density.min()) / (
                density.max() - density.min()
            )

            # Create scatter plot with colored points based on local density
            ax.scatter(x, y, c=normalized_density, cmap="viridis", s=10)
            ax.set_xlabel(ch1_choice)
            ax.set_ylabel(ch2_choice)

            # fix axis labels
            if ch1_choice in W.fluorescent_channel_names:
                XTicks, XTickLabel = _get_ticks(W.LinearRange)
                ax.set_xticks(XTicks)
                ax.set_xticklabels(XTickLabel)
            if ch2_choice in W.fluorescent_channel_names:
                XTicks, XTickLabel = _get_ticks(W.LinearRange)
                ax.set_yticks(XTicks)
                ax.set_yticklabels(XTickLabel)

            # put buttons in place. polygon and quad for scatter
            [act.setVisible("scat" in key) for key, act in widget.added_actions.items()]

        # if ch1 and ch2 are the same make a histogram
        elif ch1_choice == ch2_choice:
            x = sample_d[ch1_choice]
            # h, c = np.histogram(x, 100)
            # ax.plot((c[:-1] + c[1:]) / 2, h)
            ax.hist(x, 100)
            ax.set_xlabel(ch1_choice)
            if ch1_choice in W.fluorescent_channel_names:
                XTicks, XTickLabel = _get_ticks(W.LinearRange)
                ax.set_xticks(XTicks)
                ax.set_xticklabels(XTickLabel)

            # put buttons in place. span and thresh for hist
            [act.setVisible("hist" in key) for key, act in widget.added_actions.items()]

        fc.draw()


    # widget._gates=[]

    @widget.channel_1.changed.connect
    def _on_ch1_changed():
        widget()

    @widget.channel_2.changed.connect
    def _on_ch2_changed():
        widget()

    @widget.sample.changed.connect
    def _on_sample_changed():
        widget.population.choices = list(
            W.groups[widget.group.value].samples[widget.sample.value].gates.keys()
        )
        popname = widget.population.value
        widget.pop_tree_widget.clear()
        add_tree_to_qtreewidget(widget.pop_tree_widget, W.groups[widget.group.value].samples[widget.sample.value].gates[0])
        set_item_by_name(widget.pop_tree_widget,item_name=popname)
        widget()

    @widget.group.changed.connect
    def _on_group_changed():
        widget.population.choices = list(
            W.groups[widget.group.value].samples[widget.sample.value].gates.keys()
        )
        widget.sample.choices = list(W.groups[widget.group.value].samples.keys())
        widget.pop_tree_widget.clear()
        add_tree_to_qtreewidget(widget.pop_tree_widget, W.groups[widget.group.value].samples[widget.sample.value].gates[0])
        set_item_by_name(widget.pop_tree_widget)
        widget()

    @widget.population.changed.connect
    def _on_pop_changed():
        # widget._gates.set_visible(False)
        # widget._gates.set_active(False)
        widget()
        pass

    # make updating label showing available gates
    gate_labels = Label(value="")
    widget.append(gate_labels)


    ### make toolbar and special buttons
    toolbar = NavigationToolbar(fc, widget.parent)
    widget.added_actions = {}

    act_range = toolbar.addAction(
        qta.icon("mdi.arrow-expand-horizontal"), "add span gate"
    )
    act_range.triggered.connect(_onspanbutton)
    widget.added_actions["hist_act_range"] = act_range

    act_thresh = toolbar.addAction(
        qta.icon("ph.arrows-out-line-horizontal"), "add threshold gate"
    )
    act_thresh.triggered.connect(lambda: print("action triggered"))
    widget.added_actions["hist_act_thresh"] = act_thresh

    act_poly = toolbar.addAction(qta.icon("mdi.shape-polygon-plus"), "add polygon gate")
    act_poly.triggered.connect(lambda: print("action triggered"))
    widget.added_actions["scat_act_poly"] = act_poly

    act_quad = toolbar.addAction(qta.icon("ph.arrows-out-light"), "add quadrants gate")
    act_quad.triggered.connect(lambda: print("action triggered"))
    widget.added_actions["scat_act_quad"] = act_quad

    ### Make population tree widget
    from PyQt5.QtWidgets import QVBoxLayout, QTreeWidget, QTreeWidgetItem
    widget.pop_tree_widget = QTreeWidget()

    def _on_item_select_change():
        selected_items = widget.pop_tree_widget.selectedItems()
        if selected_items:
            widget.population.value = selected_items[0].text(0)
        else:
            set_item_by_name(widget.pop_tree_widget)
        
    
    widget.pop_tree_widget.itemSelectionChanged.connect(_on_item_select_change)
    add_tree_to_qtreewidget(widget.pop_tree_widget, W.groups[widget.group.value].samples[widget.sample.value].gates[0])
    set_item_by_name(widget.pop_tree_widget)
    


    ### layout
    container = Container(layout="horizontal")
    layout = container.native.layout()


    l = QVBoxLayout()
    l.addWidget(toolbar)
    l.addWidget(fc)
    widget.root_native_widget.layout().addWidget(widget.pop_tree_widget)

    layout.addWidget(widget.native)  # adding native, because we're in Qt
    layout.addLayout(l)

    container.show()
    matplotlib.use("Qt5Agg")
    widget.call_button.visible = False
    widget.max_width = 250
    widget()
    container.pop_tree_widget = widget.pop_tree_widget
    return container
