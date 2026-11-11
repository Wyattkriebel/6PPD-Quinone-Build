"""Launch from repository root: streamlit run dashboard/app.py -- --root ."""
import argparse
import json
from pathlib import Path
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

p=argparse.ArgumentParser();p.add_argument('--root',default='.');args,_=p.parse_known_args()
root=Path(args.root).resolve()
st.set_page_config(page_title='6PPD-Q candidate evidence',layout='wide')
st.title('6PPD-Q candidate evidence')
st.caption('Geometry, model confidence, approximate energies, and experimental affinity are distinct evidence types.')
path=root/'results/candidate_tables/candidates.csv'
reference_path=root/'results/reference_tables/reference_sequence_designs.csv'
is_reference=False
if reference_path.exists():
    dataset=st.sidebar.selectbox('Dataset',['6PPD-Q candidates','Reference inference'])
    is_reference=dataset=='Reference inference'
    if is_reference:
        path=reference_path
        st.warning('Reference DNA/zinc-context sequence inference. These structures are copied input backbones with designed residue names, not predicted folds or sidechain-packed complexes. Sequence-model confidence is not pLDDT or binding confidence. These are not 6PPD-Q candidates.')
if not path.exists():
    st.info('No candidate table yet. Run the validated pipeline first.');st.stop()
table=pd.read_csv(path)
if table.empty:
    st.info('No scientific candidates generated. Phase I reproduction has not produced validated candidates.');st.dataframe(table);st.stop()
if 'evidence_class' in table:
    options=sorted(table.evidence_class.dropna().unique())
    choices=st.sidebar.multiselect('Evidence classes',options,default=[x for x in options if x!='synthetic_fixture'])
    table=table[table.evidence_class.isin(choices)]
st.dataframe(table,width="stretch")
if table.empty: st.stop()
selected=st.selectbox('Reference sequence' if is_reference else 'Candidate',table.candidate_id)
row=table.loc[table.candidate_id==selected].iloc[0]
st.json(row.dropna().to_dict())

def safe_path(value):
    if not isinstance(value,str) or not value: return None
    path=Path(value)
    if not path.is_absolute(): path=root/path
    path=path.resolve()
    if not path.is_relative_to(root) or not path.is_file(): return None
    return path

structure=safe_path(row.get('source_structure'))
if structure:
    text=structure.read_text()
    st.download_button('Download sequence-labeled copied backbone' if is_reference else 'Download original structure',text,file_name=structure.name)
    try:
        import py3Dmol
        viewer=py3Dmol.view(width=900,height=540)
        viewer.addModel(text,'pdb');viewer.setStyle({'hetflag':False},{'cartoon':{'color':'spectrum'}})
        viewer.setStyle({'hetflag':True},{'stick':{'colorscheme':'greenCarbon'}})
        viewer.addStyle({'within':{'distance':5,'sel':{'hetflag':True}},'hetflag':False},{'stick':{}})
        viewer.zoomTo();components.html(viewer._make_html(),height=560)
        st.caption('Viewer JavaScript requires browser network access to 3Dmol. Structure downloads work offline.')
    except ImportError:
        st.warning('Install py3Dmol to enable the interactive viewer.')
for field in ('metrics_path','provenance_path'):
    path=safe_path(row.get(field))
    if path:
        st.download_button(f'Download {field}',path.read_bytes(),file_name=path.name)
        if path.suffix=='.json':
            with st.expander(field): st.json(json.loads(path.read_text()))
metric_names=('native_sequence_identity','mutation_count','overall_sequence_model_confidence','ligand_context_sequence_model_confidence') if is_reference else ('clash_pairs','backbone_rmsd_A','plddt','docking_score_kcal_mol')
metrics=[x for x in metric_names if x in table and table[x].notna().any()]
if metrics:
    metric=st.selectbox('Compare one evidence type',metrics)
    st.bar_chart(table.set_index('candidate_id')[metric])

if is_reference:
    st.caption("Sequence-model values describe sampled residue preferences; they do not validate folding or affinity.")
    mutation_path=root/"results/reference_tables/sequence_mutations.csv"
    if mutation_path.exists():
        mutations=pd.read_csv(mutation_path)
        with st.expander("Residue-by-residue changes from the reference sequence"):
            st.dataframe(mutations[mutations.candidate_id==selected],width="stretch")
    fasta_path=safe_path(row.get("source_fasta"))
    if fasta_path:
        st.download_button("Download original inference FASTA",fasta_path.read_bytes(),file_name=fasta_path.name)
