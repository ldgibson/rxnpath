import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from adjustText import adjust_text


class ReactionDiagram(nx.DiGraph):
    def __init__(self, colors=['black', 'blue', 'green',
                               'red', 'purple', 'orange'], **kwargs):
        nx.DiGraph.__init__(self, **kwargs)
        self.subgraphs = {}
        self.step_length = None
        self.colors = colors.copy()
        return

    def add_pathway(self, energies, labels, name, color=None, positions=[]):
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
        for (u, v) in edges:
            assert label in [u, v],\
                "'{}' not included in the edge declaration.".format(label)
            if u == label:
                # New label doesn't exist yet.
                pass
            else:
                assert self.has_node(u), "'{}' does not exist".format(u)

            if v == label:
                # New label doesn't exist yet.
                pass
            else:
                assert self.has_node(v), "'{}' does not exist".format(v)

        self.add_node(label, color=color)
        if subgraph_name in self.subgraphs.keys():
            new_subgraph = self.subgraphs[subgraph_name]['labels'] + [label]
            self.subgraphs[subgraph_name]['graph'] = self.subgraph(new_subgraph)
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
        self.add_edges_from(edges)
        return


    def prepare_diagram(self, step_width=3, between=1):
        self.step_length = step_width / (step_width + between)

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
                line_coords = np.array([[pos, pos + self.step_length],
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
                     fname=None, prefix=None, saveparams=dict(),
                     state_line_attr=dict(linewidth=3, linestyle='-'),
                     edge_line_attr=dict(linewidth=1, linestyle='--'),
                     show_energies=False,
                     show_positions=False, adjust=dict()):
        fig, ax = plt.subplots(figsize=figsize)
        text = []
        for n, data in self.nodes(data=True):
            ax.plot(*data['line'], color=data['color'], **state_line_attr)
            if show_energies:
                label = n + " ({:.2f})".format(data['label_coords'][1, 0])
            else:
                label = n
            text.append(ax.text(*data['label_coords'], label,
                                fontsize=fontsize))

        for n1, n2 in self.edges:
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
            pos = [i + self.step_length / 2 for i in range(max_distance)]
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
