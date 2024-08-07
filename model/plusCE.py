import warnings
import argparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy.stats import spearmanr
from Bio import SeqIO
import subprocess
import gzip
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os

warnings.filterwarnings("ignore")


STOPS = set(["TAA", "TAG", "TGA"])
STARTS = set(["ATG", "TTG", "GTG", "CTG"])


parser = argparse.ArgumentParser()
parser.add_argument("--data", help="raw data path",required=True)
parser.add_argument("--annot", help="annotation file name",required=True)
parser.add_argument("--feature", help="prepared features path",required=True)
parser.add_argument("--rna", help="min rna-seq rpkm", default=5)
parser.add_argument("--ribo", help="min ribo-seq rpkm", default=0.1)
parser.add_argument("--querydb", help="query database path",required=True)

args = parser.parse_args()


features_in = os.path.join(args.feature, "input.fa.sparseFeature.txt.gz")

with gzip.open(features_in, "rb") as f:
    data = f.read()
    data = data.decode("utf-8")
    data = data.split("\n")
    data = data[:-1]
    data = [i.split("\t") for i in data]
    data = [(int(i[0]), int(i[1]), float(i[2])) for i in data]

row = [i[0] for i in data]
col = [i[1] for i in data]
value = [i[2] for i in data]

X = sp.csr_matrix((value, (row, col)), shape=(57415, 5493))


annot_in = os.path.join(args.data, args.annot)
min_rna_rkpm = float(args.rna)
min_ribo_rkpm = float(args.ribo)
df = pd.read_csv(annot_in, sep=" ", index_col=0)
df = df[(df["rpkm_rnaseq"] > min_rna_rkpm) & (df["rpkm_riboseq"] > min_ribo_rkpm)]

rownames = os.path.join(args.feature, "input.fa.sparseFeature.rowname")
df_feat = pd.read_csv(rownames, header=None, sep="\t")

query_enst = [i.split(".")[0] for i in set(df_feat[1]).intersection(set(df.index))]

train_idx = []
train_idx2enst = {}
train_enst2idx = {}

for idx, enst in enumerate(df_feat[1]):
    enst = enst.split(".")[0]
    if enst not in query_enst:
        continue
    train_idx.append(idx)
    train_idx2enst[idx] = enst
    train_enst2idx[enst] = idx

train_idx = sorted(train_idx)
train_enst = [train_idx2enst[i] for i in train_idx]

input_fasta = os.path.join(args.feature, "input.filter.fa")
records = SeqIO.parse(input_fasta, "fasta")
enst2utr = {}
for record in records:
    enst = record.id.split(".")[0]
    if enst not in query_enst:
        continue
    enst2utr[enst] = str(record.seq).upper()

with open("tmp.fa", "w") as w:
    for enst in train_enst:
        w.write(">%s\n%s\n" % (enst, enst2utr[enst]))


db = args.querydb
cmd = f"blastn -query tmp.fa -db {db}/ac1 -outfmt 6 -out align_class1.out & blastn -query tmp.fa -db {db}/ac5 -outfmt 6 -out align_class5.out"
subprocess.call(cmd, shell=True)


newX = X[train_idx, :]

df = pd.read_csv("align_class1.out", sep="\t", header=None)
lowphast = list(df[0])

df = pd.read_csv("align_class5.out", sep="\t", header=None)
highphast = list(df[0])

FeaturesX = np.zeros((len(train_idx), 5495))
FeaturesX[:, :-2] = newX.todense()

for row, enst in enumerate(train_enst):
    utr = enst2utr[enst]
    ustart_cnt = 0
    for start in STARTS:
        ustart_cnt += utr.count(start)
    ustop_cnt = 0
    for stop in STOPS:
        ustop_cnt += utr.count(stop)

    FeaturesX[row, 9] = ustart_cnt
    FeaturesX[row, 10] = ustop_cnt

    FeaturesX[row, -2] = lowphast.count(enst)
    FeaturesX[row, -1] = highphast.count(enst)

enst2te = {}
for enst, te_log in zip(df.index, df["te"]):
    enst = enst.split(".")[0]
    if enst not in train_enst2idx:
        continue
    row = train_enst2idx[enst]
    enst2te[enst] = te_log

Y = np.array([enst2te[i] for i in train_enst])
Y = np.log(Y)


kf = KFold(n_splits=10, shuffle=True, random_state=42)

spearman_rf = []

for idx, train_index, test_index in enumerate(kf.split(FeaturesX)):

    X_train, X_test = FeaturesX[train_index], FeaturesX[test_index]
    y_train, y_test = Y[train_index], Y[test_index]

    rf_model = RandomForestRegressor(
        n_estimators=200,
        max_features="sqrt",
        min_samples_leaf=5,
        bootstrap=False,
        random_state=42,
        n_jobs=40,
    )

    rf_model.fit(X_train, y_train)

    y_pred = rf_model.predict(X_test)

    spearman_rf.append(spearmanr(y_test, y_pred)[0])

    print(f"Fold {idx}, spearman: {spearman_rf[-1]}")

print(
    f"mean spearman: {np.mean(spearman_rf)}, std: {np.std(spearman_rf)}"
)
