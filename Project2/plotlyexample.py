from plotly.offline import init_notebook_mode, iplot
from IPython.display import display, HTML

init_notebook_mode(connected=True)

figure = {'data': [{'x': [0, 1], 'y': [0, 1]}],
          'layout': {'xaxis': {'range': [0, 5], 'autorange': False},
                     'yaxis': {'range': [0, 5], 'autorange': False},
                     'title': 'Start Title'},
          'frames': [{'data': [{'x': [1, 2], 'y': [1, 2]}]},
                     {'data': [{'x': [1, 4], 'y': [1, 4]}]},
                     {'data': [{'x': [3, 4], 'y': [3, 4]}],
                      'layout': {'title': 'End Title'}}]}

iplot(figure)