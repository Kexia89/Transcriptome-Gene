import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd,numpy as np
from sklearn.metrics import roc_curve,precision_recall_curve

def generate(root):
    root=Path(root); out=root/'figures'; out.mkdir(exist_ok=True); m=json.loads((root/'metrics.json').read_text()); e=pd.read_csv(root/'evaluation.csv'); n=pd.read_csv(root/'nodes.csv'); edges=pd.read_csv(root/'edges.csv')
    plt.style.use('seaborn-v0_8-whitegrid'); names={'isolation_forest':'Isolation Forest','lof':'LOF','pca_error':'PCA reconstruction','ensemble':'Ensemble'}
    def save(fig,name):
        fig.tight_layout(); fig.savefig(out/(name+'.png'),dpi=200,bbox_inches='tight'); fig.savefig(out/(name+'.svg'),bbox_inches='tight'); plt.close(fig)
    table=pd.DataFrame(m['synthetic_holdout']).T.drop(columns='confusion_matrix'); table.to_csv(root/'model_metrics.csv')
    fig,ax=plt.subplots(figsize=(10,5)); table.astype(float).rename(index=names,columns={'roc_auc':'ROC-AUC','pr_auc_average_precision':'AP','f1_at_95':'F1'}).plot.bar(ax=ax,rot=0); ax.set_ylim(0,1.12); ax.set_title('Synthetic holdout | detector performance'); ax.legend(ncol=3); save(fig,'model_comparison')
    for mode in ['roc','pr']:
        fig,ax=plt.subplots(figsize=(7,5))
        for col,label in names.items():
            if mode=='roc': x,y,_=roc_curve(e.label,e[col]); metric=m['synthetic_holdout'][col]['roc_auc']
            else: y,x,_=precision_recall_curve(e.label,e[col]); metric=m['synthetic_holdout'][col]['pr_auc_average_precision']
            ax.plot(x,y,label=f'{label} ({metric:.3f})')
        if mode=='roc': ax.plot([0,1],[0,1],'--',color='gray'); ax.set(xlabel='False positive rate',ylabel='True positive rate',title='ROC | synthetic holdout')
        else: ax.axhline(e.label.mean(),ls='--',color='gray'); ax.set(xlabel='Recall',ylabel='Precision',title='Precision–recall | synthetic holdout')
        ax.legend(); save(fig,mode+'_curves')
    fig,ax=plt.subplots(figsize=(8,5)); ax.hist(n.anomaly_score,bins=40); ax.axvline(95,color='red',ls='--'); ax.set(xlabel='Mean empirical percentile',ylabel='Genes',title='Ensemble score distribution'); save(fig,'anomaly_distribution')
    fig,ax=plt.subplots(figsize=(8,5)); s=pd.DataFrame(m['network']['block_deletion_stability']); s.index=s.removed_columns.map(lambda x:f'{x[0]}–{x[-1]}'); s[['edge_jaccard','community_ari_shared_nonisolates']].plot.bar(ax=ax,rot=0); ax.set_ylim(0,1); ax.set_title('Network sensitivity to time-block deletion'); save(fig,'network_stability')
    fig,ax=plt.subplots(figsize=(7,5)); ax.bar(['PC1','PC2','PC3'],np.array(m['pca_variance_ratio'])*100); ax.set(ylabel='Explained variance (%)',title='PCA dimensionality reduction'); save(fig,'pca_variance')
    fig=plt.figure(figsize=(9,7)); ax=fig.add_subplot(111,projection='3d'); ax.scatter(n.x,n.y,n.z,c=n.community,cmap='turbo',s=5,alpha=.55); f=n[n.anomaly_score>=95]; ax.scatter(f.x,f.y,f.z,c='red',s=15,label='Ensemble ≥95'); ax.set_title(f"{m['embedding']} expression space (not physical positions)"); ax.legend(); save(fig,'gene_universe')
    fig,ax=plt.subplots(figsize=(7,5)); ax.hist(edges.correlation,bins=40); ax.set(xlabel='Pearson correlation',ylabel='Edges',title='Signed co-expression associations'); save(fig,'edge_correlations')
    fig,ax=plt.subplots(figsize=(8,5)); n[n.community>=0].community.value_counts().head(15).sort_values().plot.barh(ax=ax); ax.set(xlabel='Genes',ylabel='Community',title='Largest Louvain communities'); save(fig,'communities')
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--artifacts',default='artifacts'); generate(p.parse_args().artifacts)
