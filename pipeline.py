"""Reproducible exploratory transcriptomics pipeline."""
import argparse, json, hashlib, platform
import joblib
from importlib.metadata import version
from inference import raw_scores
from sklearn.metrics import precision_score, recall_score, confusion_matrix
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import rankdata, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score, roc_auc_score, adjusted_rand_score, f1_score
from sklearn.manifold import trustworthiness


def features(x):
    d=np.diff(x,axis=1); z=x-x.mean(1,keepdims=True)
    power=abs(np.fft.rfft(z,axis=1))**2
    return np.column_stack([x.mean(1),x.std(1),np.ptp(x,axis=1),x.min(1),x.max(1),d.std(1),abs(d).max(1),np.mean(abs(d),axis=1),np.sum(z[:,:-1]*z[:,1:],axis=1)/(np.sum(z*z,axis=1)+1e-12),power[:,1:].max(1)/(power[:,1:].sum(1)+1e-12)])


def graph(x, ids, threshold=.8, k=10):
    c=np.nan_to_num(np.corrcoef(x)); np.fill_diagonal(c,0)
    g=nx.Graph(); g.add_nodes_from(ids)
    for i,row in enumerate(c):
        for j in np.argsort(-abs(row),kind='stable')[:k]:
            if abs(row[j])>=threshold:
                g.add_edge(ids[i],ids[j],weight=float(abs(row[j])),correlation=float(row[j]))
    return g


def communities(g,seed):
    if not g.number_of_edges(): return []
    # Isolates have no evidence of community membership.
    h=g.subgraph([v for v in g if g.degree(v)>0]); return list(nx.community.louvain_communities(h,weight='weight',seed=seed))


def fit_detectors(x,seed):
    fs=StandardScaler().fit(features(x)); p=StandardScaler().fit(x)
    a=fs.transform(features(x)); b=p.transform(x)
    iso=IsolationForest(n_estimators=150,random_state=seed,n_jobs=1).fit(a)
    lof=LocalOutlierFactor(n_neighbors=25,novelty=True).fit(a)
    pc=PCA(n_components=8,random_state=seed).fit(b)
    return {'feature_scaler':fs,'expression_scaler':p,'isolation_forest':iso,'lof':lof,'pca_anomaly':pc}


def percentile(raw,ref):
    return np.column_stack([100*np.searchsorted(np.sort(ref[:,i]),raw[:,i],side='right')/len(ref) for i in range(raw.shape[1])])


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',type=Path,required=True); ap.add_argument('--out',type=Path,default=Path('artifacts')); ap.add_argument('--seed',type=int,default=42); ap.add_argument('--threshold',type=float,default=.8); ap.add_argument('--k',type=int,default=10); ap.add_argument('--umap',action='store_true'); args=ap.parse_args()
    if not 0<args.threshold<=1 or args.k<1: ap.error('Invalid graph parameters')
    out=args.out; out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.input,index_col=0); df.columns=df.columns.astype(int); df=df.sort_index(axis=1); x=df.to_numpy(float); ids=df.index.astype(str).tolist()
    if len(set(ids))!=len(ids) or not np.isfinite(x).all() or df.columns.has_duplicates or x.shape[1]!=23 or len(x)<60: raise ValueError('Invalid data: duplicate IDs, nonfinite values or too few timepoints')
    rng=np.random.default_rng(args.seed)
    g=graph(x,ids,args.threshold,args.k); groups=communities(g,args.seed); labels={v:i for i,s in enumerate(groups) for v in s}; lab=np.array([labels.get(v,-1) for v in ids]); pr=nx.pagerank(g,weight='weight');
    edges=pd.DataFrame([(a,b,d['correlation'],d['weight']) for a,b,d in g.edges(data=True)],columns=['source','target','correlation','weight']); edges.to_csv(out/'edges.csv',index=False)
    pca=PCA(n_components=3).fit(x); coords=pca.transform(x); method='PCA'
    if args.umap:
        import umap
        coords=umap.UMAP(n_components=3,random_state=args.seed,n_jobs=1).fit_transform(x); method='UMAP'
    train,test=train_test_split(np.arange(len(x)),test_size=.25,random_state=args.seed)
    train,cal=train_test_split(train,test_size=.20,random_state=args.seed)
    models=fit_detectors(x[train],args.seed); scorer=lambda y: raw_scores(models,y)
    ref=scorer(x[cal]); raw=scorer(x); scores=percentile(raw,ref); ensemble=scores.mean(1)
    # Paired synthetic corruption; original held-out profiles are proxy negatives, not verified normal genes.
    clean=x[test].copy(); corrupt=clean.copy(); n=len(clean); scale=np.maximum(clean.std(1),.15); kind=np.arange(n)%3
    for i in range(n):
        if kind[i]==0: corrupt[i,rng.integers(x.shape[1])]+=6*scale[i]
        elif kind[i]==1: corrupt[i]+=np.linspace(0,5*scale[i],x.shape[1])
        else: corrupt[i]=rng.permutation(clean[i])
    y=np.r_[np.zeros(n),np.ones(n)]; eval_scores=percentile(scorer(np.vstack([clean,corrupt])),ref); eval_scores=np.column_stack([eval_scores,eval_scores.mean(1)])
    names=['isolation_forest','lof','pca_error','ensemble']; synthetic={name:{'pr_auc_average_precision':float(average_precision_score(y,s)), 'roc_auc':float(roc_auc_score(y,s)), 'f1_at_95':float(f1_score(y,s>=95))} for name,s in zip(names,eval_scores.T)}
    for name,s in zip(names,eval_scores.T):
        synthetic[name].update(precision=float(precision_score(y,s>=95,zero_division=0)),recall=float(recall_score(y,s>=95,zero_division=0)),confusion_matrix=confusion_matrix(y,s>=95,labels=[0,1]).tolist())
    ev=pd.DataFrame(eval_scores,columns=names); ev['label']=y.astype(int); ev['gene']=np.tile(np.asarray(ids)[test],2); ev['perturbation']=np.r_[np.repeat('original',n),np.array(['spike','drift','shuffle'])[kind]]; ev.to_csv(out/'evaluation.csv',index=False)
    splits=pd.DataFrame({'gene':ids,'split':'train'}); splits.loc[cal,'split']='calibration'; splits.loc[test,'split']='test'; splits.to_csv(out/'splits.csv',index=False)

    base_edges={frozenset(e) for e in g.edges}; stability=[]
    for start in [0,7,14,20]:
        keep=np.ones(x.shape[1],bool); keep[start:start+3]=False; h=graph(x[:,keep],ids,args.threshold,args.k); he={frozenset(e) for e in h.edges}; hl={v:i for i,s in enumerate(communities(h,args.seed)) for v in s}; other=np.array([hl.get(v,-1) for v in ids]); mask=(lab>=0)&(other>=0); stability.append({'removed_columns':df.columns[start:start+3].tolist(),'edge_jaccard':len(base_edges&he)/max(1,len(base_edges|he)), 'community_ari_shared_nonisolates':float(adjusted_rand_score(lab[mask],other[mask])) if mask.sum()>1 else None,'shared_nonisolates':int(mask.sum())})
    nodes=pd.DataFrame({'gene':ids,'x':coords[:,0],'y':coords[:,1],'z':coords[:,2],'community':lab,'degree':[g.degree(v) for v in ids],'pagerank':[pr[v] for v in ids],'anomaly_score':ensemble,'mean':x.mean(1),'std':x.std(1)});
    for i,name in enumerate(names[:3]): nodes[name]=scores[:,i]
    nodes.to_csv(out/'nodes.csv',index=False); df.to_csv(out/'expression.csv');
    connected=[v for v in g if g.degree(v)>0]; h=g.subgraph(connected)
    metrics={'seed':args.seed,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),'genes':len(x),'timepoints':x.shape[1],'missing':0,'embedding':method,'pca_variance_ratio':pca.explained_variance_ratio_.tolist(),'trustworthiness_10':float(trustworthiness(x,coords,n_neighbors=10)),'network':{'edges':g.number_of_edges(),'isolates':len(x)-len(connected),'communities':len(groups),'density':nx.density(g),'modularity':float(nx.community.modularity(h,groups,weight='weight')) if groups else None,'threshold':args.threshold,'top_k_union':args.k,'block_deletion_stability':stability},'synthetic_holdout':synthetic,'anomaly_score_rank_correlations':np.asarray(spearmanr(scores).statistic).tolist(),'flagged_at_95':int((ensemble>=95).sum())}
    metrics['split']={'train':len(train),'calibration':len(cal),'test':len(test)}
    metrics['protocol']='v2 independent calibration, synthetic paired holdout'
    metrics['per_perturbation']={}
    for j,kname in enumerate(['spike','drift','shuffle']):
        mask=np.r_[kind==j,kind==j]; metrics['per_perturbation'][kname]={name:float(roc_auc_score(y[mask],s[mask])) for name,s in zip(names,eval_scores.T)}
    md=out/'models'; md.mkdir(exist_ok=True); models.update(reference=ref,timepoints=df.columns.tolist(),pca_embedding=pca)
    joblib.dump(models,md/'detectors.joblib',compress=3)
    versions={p:version(p) for p in ['numpy','pandas','scipy','scikit-learn','networkx','joblib','matplotlib']}
    (md/'manifest.json').write_text(json.dumps({'sklearn':version('scikit-learn'),'versions':versions,'python':platform.python_version(),'sha256':hashlib.sha256((md/'detectors.joblib').read_bytes()).hexdigest(),'input_sha256':metrics['input_sha256'],'seed':args.seed,'inference_embedding':'PCA'},indent=2))
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2,allow_nan=False)); print(json.dumps(metrics,indent=2))
    from figures import generate
    generate(out)
if __name__=='__main__': main()
