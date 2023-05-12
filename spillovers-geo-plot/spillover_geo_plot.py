#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import plotly.express as px

if __name__ == '__main__':
    df = pd.read_csv('spillovers_geo.csv')
    legend_title = '# pdm introductions'
    df.rename({'introductions': legend_title}, axis=1, inplace=True)

    fig = px.choropleth(df,
                        locations='state',
                        locationmode="USA-states",
                        scope="usa",
                        color=legend_title,
                        color_continuous_scale=[(0, 'white'), (1, 'blue')],
                        )
    fig.update_layout(
        coloraxis_colorbar=dict(
            title='',
        ),
        font=dict(size=16),
        height=500,
        width=700,
        margin=dict(l=0, r=0, t=0, b=0),
    )
    # fig.show()
    fig.write_image('spillovers_geo.pdf')
