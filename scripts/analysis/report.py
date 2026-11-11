"""Render an auditable candidate CSV as standalone HTML; no invented metric values."""
import argparse
import csv
import html
from pathlib import Path

FIELDS=['candidate_id','evidence_class','parent_id','strategy','ligand_id','sequence','source_structure','checkpoint_sha256','seed','config_sha256','provenance_path','metrics_path','clash_pairs','backbone_rmsd_A','plddt','docking_score_kcal_mol']

def write_report(table, output):
    with Path(table).open() as handle:
        reader=csv.DictReader(handle); columns=reader.fieldnames or FIELDS; rows=list(reader)
    heading=''.join('<th>'+html.escape(x)+'</th>' for x in columns)
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(row.get(x) or 'not measured'))+'</td>' for x in columns)+'</tr>' for row in rows)
    content='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Binder candidate evidence</title>
<style>body{font:16px system-ui;margin:2rem;color:#172b3a}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:.6rem;text-align:left}th{background:#e4f0f1}section{overflow:auto}p{max-width:80ch}</style>
<h1>Binder candidate evidence</h1><p>Observed calculations must be interpreted separately from biological hypotheses. Blank fields mean not measured. Synthetic fixtures are software tests and are never candidate evidence. Docking scores are model outputs, not experimentally measured affinities. No arbitrary combined score is used.</p>'''
    content+=f'<p>Source table: {html.escape(str(Path(table).resolve()))}. Rows: {len(rows)}.</p>'
    if not rows:
        content+='<p><strong>No scientific candidates generated. Phase I reproduction remains a prerequisite.</strong></p>'
    content+='<section><table><thead><tr>'+heading+'</tr></thead><tbody>'+body+'</tbody></table></section></html>'
    with Path(output).open('x') as handle: handle.write(content)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--table',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    write_report(a.table,a.output)
