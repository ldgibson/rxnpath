import os.path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from adjustText import adjust_text


class ReactionDiagram(nx.DiGraph):
    def __init__(self, colors=['black', 'blue', 'green',
                               'red', 'purple', 'orange'], **kwargs):
        nx.DiGraph.__init__(self, **kwargs)
        self.subgraphs = {}
        self.step_size = None
        self.colors = colors.copy()
        return

    def add_pathway(self, energies, labels, name, color=None, positions=[]):
        """
        Adds a reaction pathway to reaction diagram.

        Generates a subgraph labeled as `name`. Can be accessed with
        self.subgraphs[`name`]['graph']. Other attributes include the
        labels, energies, positions, and color for each node in the
        subgraph.

        Parameters
        ----------
        energies : array-like
        labels : list
            Labels of each state that correspond to values in `energies`.
        name : str
            Name of the pathway.
        color : str, optional
            Color of the pathway. If `None`, then default colors are used.
        positions : dict or list, optional
            Dictates the locations of each state. Default places each state
            equally spaced. Indexing starts at 0. If using a dict, then the
            keys are labels and values are the position of that label.

        Examples
        --------
        >>> rxn = ReactionDiagram()
        >>> rxn.add_pathway(energies=[0, 4, 1, 2],
                            labels=['A', 'B', 'C', 'D'],
                            name='path1',
                            color='black',
                            positions=[0, 1, 2, 4])  # skips position 3
        >>> rxn.add_pathway(energies=[0, 3, 1],
                            labels=['E', 'F', 'G'],
                            name='path2',
                            color='blue',
                            positions={'E': 0, 'F': 2, 'G': 4})  # skips 1 + 3
        """
        if len(labels) == len(energies):
            pass
        else:
            raise Exception("Length of labels must match length of energies.")

        if positions:
            if len(labels) == len(positions):
                pass
            else:
                raise Exception("Length of positions must "
                                "match length of labels.")

            if isinstance(positions, list):
                positions = dict((label, i) for i, label
                                 in zip(positions, labels))
            elif isinstance(positions, dict):
                pass
            else:
                raise Exception("Positions must be a list or dictionary.")
        else:
            positions = dict((label, i) for i, label in enumerate(labels))

        assert set(labels) == set(positions.keys()), \
            "Position keys and labels must match."

        if color:
            pass
        else:
            color = self.colors.pop(0)
        self.add_nodes_from(labels, color=color)
        self.add_path(labels, color=color)
        self.subgraphs[name] = dict(graph=self.subgraph(labels),
                                    color=color,
                                    labels=labels.copy(),
                                    energies=energies,
                                    positions=positions)
        return

    def add_state(self, label, energy, position,
                  edges, subgraph_name, color=None):
        """
        Adds a new state in the reaction diagram.

        Parameters
        ----------
        label : str
        energy : float
        position : int or float
        edges : list of tuple of str
            List of tuples that have names of nodes to be connected.
            Edges are directed, so first node name is the source and
            the second node name is the target.
            E.g., edge=[('C', 'new_node'), ('new_node', 'E')] would connect
            'new_node' to both 'C' and 'E' in the order: C -> new_node -> E.
        subgraph_name : str
            Name of the subgraph to which the state belongs. Can be a
            subgraph that already exists, or a new subgraph.
        color : str, optional
            Name of the color for the new state.
        """
        for (u, v) in edges:
            assert label in [u, v],\
                "'{}' not included in the edge declaration.".format(label)
            if u == label:
                # Check if label already exists in graph.
                assert not self.has_node(u), "'{}' already exists".format(u)
            else:
                assert self.has_node(u), "'{}' does not exist".format(u)

            if v == label:
                # Check if label already exists in graph.
                assert not self.has_node(v), "'{}' already exists".format(v)
            else:
                assert self.has_node(v), "'{}' does not exist".format(v)

        self.add_node(label, color=color)
        if subgraph_name in self.subgraphs.keys():
            new_subgraph = self.subgraphs[subgraph_name]['labels'] + [label]
            self.subgraphs[subgraph_name]['graph'] = \
                self.subgraph(new_subgraph)
            self.subgraphs[subgraph_name]['labels'].append(label)
            self.subgraphs[subgraph_name]['positions'][label] = position
            self.subgraphs[subgraph_name]['energies'] = \
                np.append(self.subgraphs[subgraph_name]['energies'], energy)
        else:
            self.subgraphs[subgraph_name] = dict(graph=self.subgraph(label),
                                                 color=color,
                                                 labels=[label],
                                                 energies=np.array([energy]),
                                                 positions={label: position})
        self.add_edges_from(edges, color=color)
        return

    def prepare_diagram(self, step_size):
        self.step_size = step_size

        for subG in self.subgraphs.keys():
            for i, label in enumerate(self.subgraphs[subG]['labels']):
                G = self.subgraphs[subG]['graph']
                data = G.node[label]
                data['energy'] = self.subgraphs[subG]['energies'][i]
                if self.subgraphs[subG]['positions'].values:
                    pos = self.subgraphs[subG]['positions'][label]
                else:
                    pos = i
                data['position'] = pos
                line_coords = np.array([[pos, pos + self.step_size],
                                        [data['energy'], data['energy']]],
                                       dtype=float)
                label_coords = np.array([[np.sum(line_coords[0]) / 2],
                                         [line_coords[1, 0]]])
                data['line'] = np.insert(line_coords, [1],
                                         label_coords, axis=1)
                # data['linewidth'] = linewidth
                # data['color'] = self.subgraphs[subG]['color']
                data['label_coords'] = label_coords
        return

    def plot_diagram(self, figsize=(10, 8), fontsize=12, ylabel="Energy",
                     xlabel="Reaction Progress", margins=(0.1, 0.1),
                     step_size=0.75, fname=None, prefix=None,
                     saveparams=dict(),
                     state_line_attr=dict(linewidth=3, linestyle='-'),
                     edge_line_attr=dict(linewidth=1, linestyle='--'),
                     show_energies=False, show_positions=False,
                     adjust=dict()):
        """
        Renders the reaction diagram.

        Parameters
        ----------
        figsize : tuple, optional
        fontsize : int, optional
            Font size of the state labels
        xlabel, ylabel : str, optional
        margins : float or tuple of float, optional
            Sets margins on diagram. Float specifies margin for both
            x- and y-axes. A tuple of float specifies margin for x- and
            y-axis, respectively.
        step_size : float, optional
            Controls each state's step width. Ranges from (0, 1).
            A value of `1` will set each consecutive state to start
            from the end of the previous state's x-position and will
            lead to all connecting edges being vertical lines for
            consecutive states.
        fname : str, optional
            File name of diagram to be saved. Specifying a file name
            will cause the diagram to be saved.
        prefix : str, optional
            Prefix to be prepended to `fname`. Can be used to specify
            a path different from working directory.
        saveparams : dict, optional
            kwargs of `matplotlib.pyplot.savefig()`
        state_line_attr : dict, optional
            Adjusts line style for states.
            Kwargs of `matplotlib.axes.Axes.plot()`
        edge_line_attr : dict, optional
            Adjusts line style for edges.
            Kwargs of `matplotlib.axes.Axes.plot()`
        show_energies : bool, optional
            Displays energies in labels if `True`.
        show_positions : bool, optional
            Displays numerical position locations on x-axis if `True`.
        adjust : dict, optional
            Kwargs of `adjustText.adjust_text()`
        """
        self.prepare_diagram(step_size)
        fig, ax = plt.subplots(figsize=figsize)
        text = []
        sign = 1.0
        counter = 0
        for n, data in self.nodes(data=True):
            ax.plot(*data['line'], color=data['color'], **state_line_attr)
            if show_energies:
                label = n + " ({:.2f})".format(data['label_coords'][1, 0])
            else:
                label = n

            # Shift label position in y-direction by an offset if at
            # position 0. This makes the autoalignment of text easier
            # for adjustText.
            label_x, label_y = data['label_coords']
            if data['position'] == 0:
                label_y += sign * counter
                sign *= -1.0
                counter += 0.1
            text.append(ax.text(label_x, label_y, label,
                                fontsize=fontsize))

        for n1, n2, data in self.edges(data=True):
            if data['color']:
                color = data['color']
            else:
                color = self.node[n2]['color']
            begin = self.node[n1]['line'][:, -1].reshape(2, 1)
            end = self.node[n2]['line'][:, 0].reshape(2, 1)
            edge_line = np.append(begin, end, axis=1)
            ax.plot(*edge_line, color=color, **edge_line_attr)

        if isinstance(margins, int) or isinstance(margins, float):
            ax.margins(margins)
        elif isinstance(margins, list) or isinstance(margins, tuple):
            ax.margins(*margins)
        else:
            raise TypeError("'margins' should be a number, "
                            "or a list or tuple of two numbers.")

        if show_positions:
            max_distance = 0
            for n, data in self.nodes(data=True):
                if max_distance < data['line'][0, -1] + 1:
                    max_distance = int(data['line'][0, -1] + 1)
                else:
                    pass
            pos = [i + self.step_size / 2 for i in range(max_distance)]
            ax.set_xticks(pos)
            ax.set_xticklabels([str(i) for i in range(max_distance)])
            ax.tick_params(axis='x', labelsize=16)
        else:
            ax.tick_params(axis='x',
                           which='both',
                           bottom=False,
                           top=False,
                           labelbottom=False)

        if ylabel:
            ax.tick_params(axis='y', labelsize=16)
            ax.set_ylabel(ylabel, fontsize=20)
        else:
            pass

        if xlabel:
            ax.set_xlabel(xlabel, fontsize=20)
        else:
            pass

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        adjust_text(text, autoalign='y',
                    only_move={'text': 'y', 'points': 'y'}, **adjust)
        if fname:
            if prefix:
                fname = os.path.join(prefix, fname)
            else:
                pass
            plt.savefig(fname, transparent=True, **saveparams)
        return
