import argparse,json,hashlib
from pathlib import Path
import joblib,numpy as np,pandas as pd

def raw_scores(m,x):
    from pipeline import features
    f=m['feature_scaler'].transform(features(x)); q=m['expression_scaler'].transform(x); p=m['pca_anomaly']
    return np.column_stack([-m['isolation_forest'].score_samples(f),-m['lof'].score_samples(f),np.mean((q-p.inverse_transform(p.transform(q)))**2,axis=1)])

def load_bundle(root):
    root=Path(root); meta=json.loads((root/'manifest.json').read_text()); path=root/'detectors.joblib'
    import sklearn
    if sklearn.__version__!=meta['sklearn']: raise ValueError('Install the scikit-learn version in manifest.json')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']: raise ValueError('Checksum mismatch')
    return joblib.load(path) # Trusted local artifacts only; pickle can execute code.

def predict(df,m):
    from pipeline import percentile
    df=df.copy(); df.columns=df.columns.astype(int)
    if df.columns.has_duplicates or df.index.has_duplicates or set(df.columns)!=set(m['timepoints']): raise ValueError('Expected unique IDs and exact training timepoints')
    x=df[m['timepoints']].to_numpy(float)
    if not np.isfinite(x).all(): raise ValueError('Nonfinite input')
    scores=percentile(raw_scores(m,x),m['reference']); out=pd.DataFrame(scores,index=df.index,columns=['isolation_forest','lof','pca_error']); out['anomaly_score']=scores.mean(1); out['flagged']=out.anomaly_score>=95
    out[['x','y','z']]=m['pca_embedding'].transform(x); return out
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--models',default='artifacts/models'); p.add_argument('--output',default='predictions.csv'); a=p.parse_args(); predict(pd.read_csv(a.input,index_col=0),load_bundle(a.models)).to_csv(a.output)
