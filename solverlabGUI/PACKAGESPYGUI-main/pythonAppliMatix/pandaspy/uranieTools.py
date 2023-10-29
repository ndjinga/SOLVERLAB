#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import pandas as pd
from pandas.compat import range, lrange, lmap, map, zip, string_types
import pandas.compat as compat
import pandas.core.common as com

def _get_standard_kind(kind):
    return {'density': 'kde'}.get(kind, kind)

def _get_standard_colors(num_colors=None, colormap=None, color_type='default',
                         color=None):
    import matplotlib.pyplot as plt

    if color is None and colormap is not None:
        if isinstance(colormap, compat.string_types):
            import matplotlib.cm as cm
            cmap = colormap
            colormap = cm.get_cmap(colormap)
            if colormap is None:
                raise ValueError("Colormap {0} is not recognized".format(cmap))
        colors = lmap(colormap, np.linspace(0, 1, num=num_colors))
    elif color is not None:
        if colormap is not None:
            warnings.warn("'color' and 'colormap' cannot be used "
                          "simultaneously. Using 'color'")
        colors = color
    else:
        if color_type == 'default':
            # need to call list() on the result to copy so we don't
            # modify the global rcParams below
            colors = list(plt.rcParams.get('axes.color_cycle',
                                           list('bgrcmyk')))
            if isinstance(colors, compat.string_types):
                colors = list(colors)
        elif color_type == 'random':
            import random
            def random_color(column):
                random.seed(column)
                return [random.random() for _ in range(3)]

            colors = lmap(random_color, lrange(num_colors))
        else:
            raise ValueError("color_type must be either 'default' or 'random'")

    if isinstance(colors, compat.string_types):
        import matplotlib.colors
        conv = matplotlib.colors.ColorConverter()
        def _maybe_valid_colors(colors):
            try:
                [conv.to_rgba(c) for c in colors]
                return True
            except ValueError:
                return False

        # check whether the string can be convertable to single color
        maybe_single_color = _maybe_valid_colors([colors])
        # check whether each character can be convertable to colors
        maybe_color_cycle = _maybe_valid_colors(list(colors))
        if maybe_single_color and maybe_color_cycle and len(colors) > 1:
            msg = ("'{0}' can be parsed as both single color and "
                   "color cycle. Specify each color using a list "
                   "like ['{0}'] or {1}")
            raise ValueError(msg.format(colors, list(colors)))
        elif maybe_single_color:
            colors = [colors]
        else:
            # ``colors`` is regarded as color cycle.
            # mpl will raise error any of them is invalid
            pass

    if len(colors) != num_colors:
        multiple = num_colors//len(colors) - 1
        mod = num_colors % len(colors)

        colors += multiple * colors
        colors += colors[:mod]

    return colors



def parallel_coordinates(frame, class_column, cols=None, ax=None, color=None,
                         use_columns=False, xticks=None, colormap=None,
                         axvlines=True, axvlines_kwds={'linewidth':1,'color':'black'}, **kwds):
    """Parallel coordinates plotting.

    Parameters
    ----------
    frame: DataFrame
    class_column: str
        Column name containing class names
    cols: list, optional
        A list of column names to use
    ax: matplotlib.axis, optional
        matplotlib axis object
    color: list or tuple, optional
        Colors to use for the different classes
    use_columns: bool, optional
        If true, columns will be used as xticks
    xticks: list or tuple, optional
        A list of values to use for xticks
    colormap: str or matplotlib colormap, default None
        Colormap to use for line colors.
    axvlines: bool, optional
        If true, vertical lines will be added at each xtick
    axvlines_kwds: keywords, optional
        Options to be passed to axvline method for vertical lines
    kwds: keywords
        Options to pass to matplotlib plotting method

    Returns
    -------
    ax: matplotlib axis object

    Examples
    --------
    >>> from pandas import read_csv
    >>> from pandas.tools.plotting import parallel_coordinates
    >>> from matplotlib import pyplot as plt
    >>> df = read_csv('https://raw.github.com/pydata/pandas/master/pandas/tests/data/iris.csv')
    >>> parallel_coordinates(df, 'Name', color=('#556270', '#4ECDC4', '#C7F464'))
    >>> plt.show()
    """
    import matplotlib.pyplot as plt

    print("uranieTools.parallel_coordinates")
    n = len(frame)
    classes = frame[class_column].drop_duplicates()
    class_col = frame[class_column]

    if cols is None:
        df = frame.drop(class_column, axis=1)
    else:
        df = frame[cols]

    #used_legends = set([])

    ncols = len(df.columns)

    # determine values to use for xticks
    if use_columns is True:
        if not np.all(np.isreal(list(df.columns))):
            raise ValueError('Columns must be numeric to be used as xticks')
        x = df.columns
    elif xticks is not None:
        if not np.all(np.isreal(xticks)):
            raise ValueError('xticks specified must be numeric')
        elif len(xticks) != ncols:
            raise ValueError('Length of xticks must match number of columns')
        x = xticks
    else:
        x = lrange(ncols)

    if ax is None:
        ax = plt.gca()

    color_values = _get_standard_colors(num_colors=len(classes),
                                        colormap=colormap, color_type='random',
                                        color=color)

    colors = dict(list(zip(classes, color_values)))

    print("uranieTools.parallel_coordinates color done %s" % len(classes))
    for i in range(n):
        y = df.iloc[i].values
        kls = class_col.iat[i]
        #label = com.pprint_thing(kls)
        #used_legends.add(label)
        ax.plot(x, y, color=colors[kls], **kwds)
        """
        if label not in used_legends:
            used_legends.add(label)
            ax.plot(x, y, color=colors[kls], label=label, **kwds)
        else:
            ax.plot(x, y, color=colors[kls], **kwds)
        """

    if axvlines:
        for i in x:
            ax.axvline(i, **axvlines_kwds)

    ax.set_xticks(x)
    ax.set_xticklabels(df.columns)
    ax.set_xlim(x[0], x[-1])
    #ax.legend(loc='upper right')
    ax.grid()
    print("uranieTools.parallel_coordinates return")
    return ax



