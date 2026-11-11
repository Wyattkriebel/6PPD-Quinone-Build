"""Render MD coordinate-analysis outputs; never infer a simulation was performed."""
import argparse
import json
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--fixture', action='store_true', help='Prominently label synthetic data')
    args = parser.parse_args()
    data = json.loads(args.analysis.read_text())
    fig = make_subplots(rows=2, cols=2, subplot_titles=('Backbone RMSD', 'Ligand displacement and pose RMSD', 'Backbone RMSF', 'Atom-pair contact occupancy'))
    fig.add_trace(go.Scatter(x=data['time_ps'], y=data['backbone_rmsd'], name='Backbone RMSD'), row=1, col=1)
    for key in ('ligand_pose_rmsd','ligand_centroid_displacement'):
        fig.add_trace(go.Scatter(x=data['time_ps'], y=data[key], name=key), row=1, col=2)
    fig.add_trace(go.Bar(x=list(range(len(data['backbone_rmsf']))), y=data['backbone_rmsf'], name='Backbone RMSF'), row=2, col=1)
    fig.add_trace(go.Heatmap(z=data['contact_occupancy_by_atom_pair'], zmin=0, zmax=1, showscale=False, name='Contact occupancy'), row=2, col=2)
    fig.update_xaxes(title_text='Time (ps)', row=1)
    fig.update_yaxes(title_text='Distance (angstrom)', row=1)
    fig.update_xaxes(title_text='Protein atom index', row=2, col=1)
    fig.update_yaxes(title_text='RMSF (angstrom)', row=2, col=1)
    fig.update_xaxes(title_text='Ligand atom index', row=2, col=2)
    fig.update_yaxes(title_text='Protein atom index', row=2, col=2)
    label = 'SYNTHETIC FIXTURE — NO MD EXECUTED' if args.fixture else 'Coordinate analysis — verify source provenance'
    fig.update_layout(title=label, height=850, legend=dict(orientation='h', y=-.15))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        handle.write(fig.to_html(full_html=True, include_plotlyjs=True))

if __name__ == '__main__': main()
