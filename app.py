"""Interactive local research dashboard. Run: streamlit run app.py"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(page_title='GeneScope AI',layout='wide')
st.title('GeneScope AI | Transcriptome Explorer')
root=Path(st.sidebar.text_input('Artifacts directory','artifacts'))
@st.cache_data
def load(root):
    p=Path(root)
    return pd.read_csv(p/'nodes.csv'),pd.read_csv(p/'edges.csv'),pd.read_csv(p/'expression.csv',index_col=0),json.loads((p/'metrics.json').read_text())
try: n,e,x,m=load(str(root))
except FileNotFoundError:
    st.error('Run pipeline.py first.'); st.stop()
a,b,c,d=st.columns(4)
a.metric('Genes',len(n)); b.metric('Edges',len(e)); c.metric('Communities',m['network']['communities']); d.metric('Isolates',m['network']['isolates'])
st.caption('Exploratory co-expression, not physical locations or causal interactions. Anomaly scores are empirical ranks.')
explore, performance, inference_tab = st.tabs(['Gene universe','Model performance','Pretrained inference'])
with explore:
    threshold=st.sidebar.slider('Ensemble anomaly threshold',0,100,95)
    time=st.sidebar.select_slider('Time point',options=list(x.columns))
    selected=st.sidebar.selectbox('Inspect gene',n.gene.tolist())
    mode=st.sidebar.radio('Color by',['Community','Expression','Anomaly score'])
    show=st.sidebar.checkbox('Show selected gene connections',True)
    n=n.copy(); n['expression']=x.loc[n.gene,time].to_numpy(); n['group']=n.community.astype(str); n['size']=2+np.log1p(n.degree); n['flag']=np.where(n.anomaly_score>=threshold,'Flagged','Other')
    color={'Community':'group','Expression':'expression','Anomaly score':'anomaly_score'}[mode]
    fig=px.scatter_3d(n,x='x',y='y',z='z',color=color,size='size',symbol='flag',hover_name='gene',hover_data=['degree','anomaly_score','expression'],template='plotly_dark',color_continuous_scale='RdBu_r',size_max=8)
    if show:
        sub=e[(e.source==selected)|(e.target==selected)]; pos=n.set_index('gene'); xx=[]; yy=[]; zz=[]
        for row in sub.itertuples():
            for arr,col in [(xx,'x'),(yy,'y'),(zz,'z')]: arr.extend([pos.loc[row.source,col],pos.loc[row.target,col],None])
        fig.add_trace(go.Scatter3d(x=xx,y=yy,z=zz,mode='lines',line=dict(color='white',width=2),name='Selected edges'))
    fig.update_layout(height=720,scene=dict(xaxis_title='Dimension 1',yaxis_title='Dimension 2',zaxis_title='Dimension 3'))
    st.plotly_chart(fig,use_container_width=True)
    st.subheader(f'Gene inspector: {selected}'); st.line_chart(x.loc[selected].rename('Log ratio')); st.dataframe(n[n.gene==selected]); st.dataframe(e[(e.source==selected)|(e.target==selected)])
    st.subheader('Anomaly ranking'); st.dataframe(n.sort_values('anomaly_score',ascending=False).head(100)); st.download_button('Export node results',n.to_csv(index=False),'genes.csv')
with performance:
    st.header('Experimental model performance')
    st.caption('Synthetic paired holdout; not biological ground truth. Threshold 95 is fixed, not optimized.')
    table=pd.DataFrame(m['synthetic_holdout']).T.drop(columns='confusion_matrix',errors='ignore').astype(float)
    winner=table.roc_auc.idxmax(); a,b,c=st.columns(3); a.metric('Best detector',winner); b.metric('ROC-AUC',f"{table.loc[winner,'roc_auc']:.4f}"); c.metric('F1 at 95',f"{table.loc[winner,'f1_at_95']:.4f}")
    st.dataframe(table.style.format('{:.4f}'),use_container_width=True)
    st.plotly_chart(px.bar(table.reset_index().melt(id_vars='index'),x='index',y='value',color='variable',barmode='group',range_y=[0,1],labels={'index':'Detector','value':'Score','variable':'Metric'}),use_container_width=True)
    if 'per_perturbation' in m: st.subheader('ROC-AUC by perturbation'); st.dataframe(pd.DataFrame(m['per_perturbation']).T)
    for name in ['roc_curves','pr_curves','network_stability','pca_variance','anomaly_distribution','communities','edge_correlations']:
        path=root/'figures'/f'{name}.png'
        if path.exists(): st.image(str(path)); st.download_button(f'Download {name}',path.read_bytes(),path.name,key=name)
    st.download_button('Download metrics JSON',json.dumps(m,indent=2),'metrics.json')
    with st.expander('Network and experimental details'): st.json(m)
with inference_tab:
    st.header('Score compatible expression profiles')
    st.caption('No retraining. Use the same timepoints, organism and measurement protocol. Only trusted local Joblib models are loaded.')
    upload=st.file_uploader('Expression CSV (first column: gene ID)',type=['csv'])
    if upload is not None and st.button('Run inference'):
        try:
            from inference import load_bundle,predict
            result=predict(pd.read_csv(upload,index_col=0),load_bundle(root/'models'))
            st.dataframe(result); st.download_button('Download predictions',result.to_csv(),'predictions.csv')
        except Exception as exc: st.error(str(exc))
