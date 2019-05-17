import matplotlib.pyplot as plt
import numpy as np

from adjustText import adjust_text


class ReactionDiagram:
    """Simple class for generating reaction diagrams."""
    def __init__(self, figsize=(8, 6), label_font_size=12,
                 step_width=2.5, between=1):
        """Inits a `ReactionDiagram` object.

        Parameters
        ----------
        figsize : tuple of float or int, optional
            Defines the size of the figure.
        label_font_size : int, optional
            Font size for all labels in reaction diagram.
        step_width : float or int, optional
            Width of line for each state in reaction pathway.
        between : float or int, optionl
            Distance between reaction states in pathway.
        """
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.label_font_size = label_font_size
        self.step_width = step_width
        self.between = between
        self.x_points = np.array([])
        self.label_x_points = np.array([])
        self.text = []
        return

    def _generate_x_points(self, energies):
        _p = 0
        coords = []
        label_coords = []
        for _ in energies:
            coords.append(_p)
            _p += self.step_width / 2
            coords.append(_p)
            label_coords.append(_p)
            _p += self.step_width / 2
            coords.append(_p)
            _p += self.between

        self.x_points = np.array(coords)
        self.label_x_points = np.array(label_coords)
        return

    def _generate_pathway_points(self, energies):
        if self.x_points.size:
            pass
        else:
            self._generate_x_points(energies)
        pathway = []
        counter = 0
        for E in energies:
            for _ in range(3):
                pathway.append([self.x_points[counter], E])
                counter += 1
        return np.array(pathway)

    def _generate_label_points(self, pathway):
        return pathway[1::3, :]

    def add_pathway(self, energies, labels, fontsize=None, **kwargs):
        pathway = self._generate_pathway_points(energies)
        label_points = self._generate_label_points(pathway)
        self.ax.plot(pathway[:, 0], pathway[:, 1], **kwargs)
        for label, (x, y) in zip(labels, label_points):
            # offset exists so labels do not share exact same coordinates
            # in any of the dimensions so `adjust_text` can effectively
            # push them apart
            offset = np.random.rand(1)[0]
            x += offset / 100
            y += offset / 2
            if fontsize:
                pass
            else:
                fontsize = self.label_font_size
            self.text.append(self.ax.text(x, y, label, fontsize=fontsize))
        return

    def set_xlabel(self, xlabel, **kwargs):
        self.ax.set_xlabel(xlabel, **kwargs)
        return

    def set_ylabel(self, ylabel, **kwargs):
        self.ax.set_ylabel(ylabel, **kwargs)
        return

    def set_title(self, title, **kwargs):
        self.ax.set_title(title, **kwargs)

    def set_figsize(self, figsize=(8, 6)):
        self.fig.set_size_inches(figsize)
        return

    def draw_diagram(self, move_text=True, margin=0.1, **kwargs):
        self.ax.margins(margin)
        self.ax.tick_params(axis='x',
                            which='both',
                            bottom=False,
                            top=False,
                            labelbottom=False)
        if move_text:
            if 'force_text' in kwargs.keys():
                force_text = kwargs['force_text']
                del kwargs['force_text']
            else:
                force_text = (0, 0.2)

            if 'force_points' in kwargs.keys():
                force_points = kwargs['force_points']
                del kwargs['force_points']
            else:
                force_points = (0, 0.2)
            count = adjust_text(self.text, ax=self.ax, autoalign='y',
                                force_points=force_points,
                                force_text=force_text,
                                only_move={'text': 'y', 'points': 'y'},
                                **kwargs)
        else:
            count = 0
        return count
