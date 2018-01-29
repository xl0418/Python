import plotly.plotly as py
from plotly.grid_objs import Grid, Column

import time

column_1 = Column([0.9, 1.1], 'x')
column_2 = Column([1.0, 1.0], 'y')
column_3 = Column([0.8, 1.2], 'x2')
column_4 = Column([1.2, 0.8], 'y2')
column_5 = Column([0.7, 1.3], 'x3')
column_6 = Column([0.7, 1.3], 'y3')
column_7 = Column([0.6, 1.4], 'x4')
column_8 = Column([1.5, 0.5], 'y4')
column_9 = Column([0.4, 1.6], 'x5')
column_10 = Column([1.2, 0.8], 'y5')

grid = Grid([column_1, column_2, column_3, column_4, column_5,
             column_6, column_7, column_8, column_9, column_10])
py.grid_ops.upload(grid, 'points_changing_size_grid'+str(time.time()), auto_open=False)

# create figure
figure = {
    'data': [
        {
            'xsrc': grid.get_column_reference('x'),
            'ysrc': grid.get_column_reference('y'),
            'mode': 'markers',
            'marker': {'color': '#48186a', 'size': 10}
        }
    ],
    'layout': {'title': 'Growing Circles',
               'xaxis': {'range': [0, 2], 'autorange': False},
               'yaxis': {'range': [0, 2], 'autorange': False},
               'updatemenus': [{
                   'buttons': [
                       {'args': [None],
                        'label': 'Play',
                        'method': 'animate'}
               ],
               'pad': {'r': 10, 't': 87},
               'showactive': False,
               'type': 'buttons'
                }]},
    'frames': [
        {
            'data': [
                {
                    'xsrc': grid.get_column_reference('x2'),
                    'ysrc': grid.get_column_reference('y2'),
                    'mode': 'markers',
                    'marker': {'color': '#3b528b', 'size': 25}
                }
            ]
        },
        {
            'data': [
                {
                    'xsrc': grid.get_column_reference('x3'),
                    'ysrc': grid.get_column_reference('y3'),
                    'mode': 'markers',
                    'marker': {'color': '#26828e', 'size': 50}
                }
            ]
        },
        {
            'data': [
                {
                    'xsrc': grid.get_column_reference('x4'),
                    'ysrc': grid.get_column_reference('y4'),
                    'mode': 'markers',
                    'marker': {'color': '#5ec962', 'size': 80}
                }
            ]
        },
        {
            'data': [
                {
                    'xsrc': grid.get_column_reference('x5'),
                    'ysrc': grid.get_column_reference('y5'),
                    'mode': 'markers',
                    'marker': {'color': '#d8e219', 'size': 100}
                }
            ]
        }
    ]
}
py.icreate_animations(figure, 'points_changing_size'+str(time.time()))