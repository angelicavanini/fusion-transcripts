#figure 3F funnel chart 

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Data for Fusion Transcripts
sample_counts_1 = [717, 672, 51]

# Data for Normal Transcripts
stages = ['A', 'B', 'C']
sample_counts_2 = [717, 680, 212]

# Create a subplot figure with 1 row and 2 columns
fig = make_subplots(rows=1, cols=2, subplot_titles=("Fusion Transcripts", "Normal Transcripts"))

# Funnel chart for FT
funnel_1 = go.Funnel(
    y=stages,
    x=sample_counts_1,
    textinfo="value+percent initial",
    texttemplate="%{value} (%{percentInitial:.1%})",
    marker=dict(color=["#154360", "#1F618D", "#7FB3D5"])
)

# Funnel chart for NT
funnel_2 = go.Funnel(
    y=stages,
    x=sample_counts_2,
    textinfo="value+percent initial",
    texttemplate="%{value} (%{percentInitial:.1%})",
    marker=dict(color=["#154360", "#1F618D", "#7FB3D5"])
)

# Add funnel charts to the subplot figure
fig.add_trace(funnel_1, row=1, col=1)
fig.add_trace(funnel_2, row=1, col=2)

# Update layout for better visualization
fig.update_layout(
    title="Fusion Transcripts Ribosomal Profiling vs Control",
    showlegend=False
)

# Set x-axis titles
fig.update_xaxes(title_text="Number of Samples", row=1, col=1)
fig.update_xaxes(title_text="Number of Samples", row=1, col=2)

# Show the figure
fig.show()
