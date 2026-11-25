import numpy as np
from pipeline import features, graph, percentile, communities

def test_features():
    x=np.arange(46,dtype=float).reshape(2,23)
    assert features(x).shape==(2,10)
    assert np.isfinite(features(x)).all()

def test_graph_signed():
    x=np.array([np.arange(23),-np.arange(23),np.ones(23)])
    with np.errstate(invalid='ignore',divide='ignore'):
        g=graph(x,['a','b','c'],.8,2)
    assert g.number_of_edges()==1
    assert g['a']['b']['correlation']<-.99
    assert g.degree('c')==0
    assert len(communities(g,42))==1

def test_percentile():
    assert percentile(np.array([[2.]]),np.array([[1.],[2.],[3.]])).item()==100*2/3

def test_saved_inference(tmp_path):
    import joblib,pandas as pd
    from sklearn.decomposition import PCA
    from pipeline import fit_detectors
    from inference import raw_scores,predict
    x=np.random.default_rng(42).normal(size=(100,23)); m=fit_detectors(x[:60],42)
    m.update(reference=raw_scores(m,x[60:80]),timepoints=list(range(23)),pca_embedding=PCA(3).fit(x[:60]))
    df=pd.DataFrame(x[80:]); before=predict(df,m); path=tmp_path/'model.joblib'; joblib.dump(m,path)
    pd.testing.assert_frame_equal(before,predict(df,joblib.load(path)))
