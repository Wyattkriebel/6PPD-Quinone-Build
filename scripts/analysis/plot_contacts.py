"""Create an interactive distance plot from observed geometry JSON."""
import argparse
import json
import hashlib
from pathlib import Path
import plotly.graph_objects as go
p=argparse.ArgumentParser();p.add_argument('metrics');p.add_argument('--output',required=True);p.add_argument('--evidence-label',required=True);a=p.parse_args()
data=json.loads(Path(a.metrics).read_text())
c=data['contacts']
fig=go.Figure(go.Scatter(x=['/'.join(map(str,x['protein_atom'])) for x in c],y=['/'.join(map(str,x['ligand_atom'])) for x in c],mode='markers',marker={'size':12,'color':[x['distance_A'] for x in c],'colorscale':'Viridis','showscale':True,'colorbar':{'title':'Distance Å'}},text=[f"distance={x['distance_A']:.3f} Å; overlap={x['overlap_A']:.3f} Å; clash={x['clash']}" for x in c],hovertemplate='%{text}<extra></extra>'))
fig.update_layout(title=f'{a.evidence_label}: geometric contact pairs',xaxis_title='Protein atom identity',yaxis_title='Ligand atom identity',height=700)
provenance={'source_metrics':str(Path(a.metrics).resolve()), 'source_metrics_sha256':hashlib.sha256(Path(a.metrics).read_bytes()).hexdigest(), 'source_structure':data.get('source_structure'), 'parameters':data.get('parameters'), 'evidence_label':a.evidence_label}
sidecar=Path(a.output+'.provenance.json')
if Path(a.output).exists() or sidecar.exists():raise FileExistsError('Figure or provenance already exists')
with Path(a.output).open('x') as f:f.write(fig.to_html(include_plotlyjs=True,full_html=True))
with sidecar.open('x') as f:json.dump(provenance,f,indent=2)
