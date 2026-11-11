"""Import genuine LigandMPNN inference and passed interface checks into a separate table."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta_records(path):
    records=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            fields=line[1:].split(', ')
            records.append({'name':fields[0], 'metadata':dict(x.split('=',1) for x in fields[1:] if '=' in x),'sequence':''})
        elif line.strip():
            if not records:raise ValueError('Sequence before FASTA header')
            records[-1]['sequence']+=line.strip()
    return records


def collect(root, run):
    root,run=Path(root).resolve(),Path(run).resolve()
    if not run.is_relative_to(root):raise ValueError('Run must be inside repository')
    rel=lambda p:str(Path(p).resolve().relative_to(root))
    provenance_path=run/'sequence_design/adapter_provenance.json'
    provenance=json.loads(provenance_path.read_text())
    if provenance.get('status')!='success':raise ValueError('Inference did not succeed')
    fastas=list((run/'sequence_design/seqs').glob('*.fa'))
    if len(fastas)!=1:raise ValueError('Expected one reference FASTA file')
    records=fasta_records(fastas[0])
    if len(records)<2:raise ValueError('No designed sequences')
    native=records[0]['sequence']
    rows=[];mutations=[];inputs={rel(provenance_path):sha256(provenance_path),rel(fastas[0]):sha256(fastas[0])}
    for record in records[1:]:
        metadata=record['metadata'];design_id=int(metadata['id'])
        validation_path=run/f'validate_sequence_{design_id}/validation.json'
        validation=json.loads(validation_path.read_text())
        if validation.get('status')!='passed':raise ValueError('Interface validation did not pass')
        if validation['designed_sequence']!=record['sequence']:raise ValueError('FASTA and validation sequence disagree')
        if len(native)!=len(record['sequence']):raise ValueError('Sequence lengths differ')
        for path,digest in validation['sources'].items():
            original=Path(path)
            # Portable bundles may retain the original absolute paths in immutable provenance.
            if not original.exists():
                anchors=('data/reference/','results/reproduction/')
                matches=[root/(anchor+path.split(anchor,1)[1]) for anchor in anchors if anchor in path]
                if len(matches)!=1:raise ValueError(f'Cannot resolve source {path}')
                original=matches[0]
            if sha256(original)!=digest:raise ValueError(f'Source checksum changed: {path}')
            inputs[rel(original)]=digest
        inputs[rel(validation_path)]=sha256(validation_path)
        backbone=run/f'sequence_design/backbones/{record["name"]}_{design_id}.pdb'
        digest=hashlib.sha256((sha256(fastas[0])+str(design_id)+provenance['checkpoint_sha256']).encode()).hexdigest()[:12]
        identity=sum(a==b for a,b in zip(native,record['sequence']))/len(native)
        row={'candidate_id':f'reference-lmpnn-{digest}','evidence_class':'reference_sequence_inference','design_id':design_id,'strategy':'LigandMPNN reference functionality','sequence':record['sequence'],'sequence_length':len(record['sequence']),'native_sequence_identity':identity,'mutation_count':sum(a!=b for a,b in zip(native,record['sequence'])),'overall_sequence_model_confidence':float(metadata['overall_confidence']),'ligand_context_sequence_model_confidence':float(metadata['ligand_confidence']),'sequence_model_confidence_is_plddt':False,'source_structure':rel(backbone),'structure_status':'copied input backbone with designed residue names; not predicted or sidechain-packed','context_residue_names':';'.join(validation['context_residue_names']),'context_heavy_atoms':validation['context_heavy_atoms'],'checkpoint_sha256':provenance['checkpoint_sha256'],'seed':int(metadata['seed']),'sampling_temperature':float(metadata['T']),'provenance_path':rel(provenance_path),'metrics_path':rel(validation_path),'source_fasta':rel(fastas[0]),'interface_status':validation['status'],'plddt':None,'backbone_prediction_rmsd_A':None,'binding_affinity':None}
        rows.append(row)
        for index,(a,b) in enumerate(zip(native,record['sequence']),1):
            mutations.append({'candidate_id':row['candidate_id'],'sequence_position':index,'reference_aa':a,'designed_aa':b,'changed':a!=b})
    return rows,mutations,inputs,native


def write_outputs(root,run,output):
    rows,mutations,inputs,native=collect(root,run)
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    for name,records in [('reference_sequence_designs.csv',rows),('sequence_mutations.csv',mutations)]:
        with (output/name).open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    pairs=[{'candidate_a':a['candidate_id'],'candidate_b':b['candidate_id'],'ungapped_identity':sum(x==y for x,y in zip(a['sequence'],b['sequence']))/len(a['sequence'])} for i,a in enumerate(rows) for b in rows[i+1:]]
    record={'evidence_class':'reference_sequence_inference','inputs_sha256':inputs,'candidate_count':len(rows),'pairwise_sequence_identity':pairs,'limitations':['Copied backbones are not structure predictions.','Sequence-model confidence is not pLDDT or binding confidence.','Reference DNA/zinc context does not reproduce organic small-molecule binder generation.','No 6PPD-Q candidates or experimental affinities.']}
    (output/'provenance.json').write_text(json.dumps(record,indent=2))
    body='<h1>Actual LigandMPNN reference sequence inference</h1><p>DNA/zinc context. Residue labels were redesigned on copied input coordinates; no new fold, sidechain packing or binding validation was calculated. Colors indicate sequence mutations only.</p><pre>Reference  '+html.escape(native)+'\n'
    for row in rows:
        sequence=''.join(f'<span style="background:#f4b183">{b}</span>' if a!=b else b for a,b in zip(native,row['sequence']))
        body+=f'Design {row["design_id"]}   '+sequence+'\n'
    body+='</pre><p>Source files, model checkpoint, inference seed, interface validation, and checksums are linked by the adjacent CSV and provenance.json. Sequence-model confidence is not structural confidence.</p>'
    (output/'sequence_comparison.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Reference sequence inference</title><style>body{font:16px system-ui;margin:2rem}pre{font:14px monospace;line-height:2;overflow:auto}</style>'+body+'</html>')
    return record

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',default='.');p.add_argument('--run',required=True);p.add_argument('--output',default='results/reference_tables');a=p.parse_args()
    print(json.dumps(write_outputs(a.root,a.run,a.output),indent=2))
