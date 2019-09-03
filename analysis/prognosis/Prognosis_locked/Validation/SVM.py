import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.model_selection import ShuffleSplit, GridSearchCV,train_test_split
from sklearn.model_selection import cross_validate

from sksurv.datasets import load_veterans_lung_cancer
from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastSurvivalSVM,FastKernelSurvivalSVM
from sksurv.kernels import clinical_kernel
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['Status'], y['Survival_in_days'], prediction)
    return result[0]
    
  df_final = pd.read_table("full_data_validation.tsv",sep=" ")

eln = [112,113,114]
comp =list(range(88,112)) 
#comp_overlap = list(range(167,197))
age = [83]

all_gen = list(range(0,57))
tmp = df_final.iloc[:,all_gen][df_final.iloc[:,all_gen] >0].count()
gen = [df_final.columns.get_loc(c) for c in tmp[tmp>df_final.shape[0]*0.02].keys() if c in df_final]

all_cyto = list(range(57,80))
tmp = df_final.iloc[:,all_cyto][df_final.iloc[:,all_cyto] >0].count()
cyto = [df_final.columns.get_loc(c) for c in tmp[tmp>df_final.shape[0]*0.02].keys() if c in df_final]

clin=list(range(84,87))
demo=[82,83]
demo_without_age = [82]

eln_comp = eln + comp
eln_gen = eln + gen
eln_cyto = eln + cyto
eln_clin = eln + clin
eln_demo = eln + demo

# USEFUL FOR ELN COMPARISON
# with comp
eln_comp_gen = eln_comp + gen
eln_comp_cyto = eln_comp + cyto
eln_comp_clin = eln_comp + clin
eln_comp_demo = eln_comp + demo




eln_comp_gen_cyto = eln_comp_gen + cyto
eln_comp_gen_clin = eln_comp_gen + clin
eln_comp_gen_demo = eln_comp_gen + demo

eln_comp_cyto_clin = eln_comp_cyto + clin
eln_comp_cyto_demo = eln_comp_cyto + demo


eln_comp_clin_demo = eln_comp_clin + demo


eln_comp_gen_cyto_clin_demo = eln_comp_gen_cyto + clin + demo
eln_comp_gen_cyto_clin_demo_without_age = eln_comp_gen_cyto + clin + demo_without_age
# without comp


eln_gen_cyto = eln_gen + cyto
eln_gen_clin = eln_gen + clin
eln_gen_demo = eln_gen + demo


eln_cyto_clin = eln_cyto + clin
eln_cyto_demo = eln_cyto + demo

eln_clin_demo = eln_clin + demo
eln_clin_demo_without_age = eln_clin + demo_without_age


eln_gen_cyto_clin_demo = eln_gen_cyto + clin + demo

# USEFUL FOR COMP

comp_gen = comp + gen
comp_cyto = comp + cyto
comp_clin = comp + clin
comp_demo = comp + demo
comp_gen_cyto = comp_gen + cyto
comp_clin_demo = comp_clin + demo
comp_gen_cyto_clin_demo = comp_gen_cyto + clin + demo

#USEFUL FOR GEN
gen_cyto = gen + cyto
gen_clin = gen + clin
gen_demo = gen + demo
gen_clin_demo = gen_clin + demo
gen_cyto_clin_demo = gen_cyto + clin + demo

#USEFUL FOR CYTO 
cyto_clin = cyto + clin
cyto_demo = cyto + demo
gen_demo_without_age = gen + demo_without_age
cyto_clin_demo = cyto_clin + demo
cyto_gen_demo = gen_cyto + demo


clin_demo  = clin + demo



dict_features_type_final_comp = dict(zip(("cyto_gen_demo","gen_cyto_clin_demo","demo","clin","gen","cyto","comp","eln","gen","cyto","gen_cyto","eln_gen_cyto","comp_gen_cyto","eln_comp","eln_clin_demo","comp_clin_demo","eln_comp_gen_cyto_clin_demo"),
                                         (cyto_gen_demo,gen_cyto_clin_demo,demo,clin,gen,cyto,comp,eln,gen,cyto,gen_cyto,eln_gen_cyto,comp_gen_cyto,eln_comp,eln_clin_demo,comp_clin_demo,eln_comp_gen_cyto_clin_demo)))
estimator = FastSurvivalSVM(max_iter=1000, tol=1e-6, random_state=17)
param_grid = {'alpha': 10. ** np.array([-6,-5,-4,-3,-2,-1,0]),'optimizer':["rbtree"]}
cv = ShuffleSplit(n_splits=5,random_state=17)
gcv = GridSearchCV(estimator, param_grid, scoring=score_survival_model,
                   n_jobs=10, iid=False, refit=True,
                   cv=cv)
df=pd.DataFrame(columns=dict_features_type_final_comp.keys())
for key,item in dict_features_type_final_comp.items():
    x = df_final.iloc[:,item]
    y = np.array(list(zip(df_final.OS_Status, df_final.OS)),dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
    ci=[]
    for i in range(25):
        X_train, X_test, y_train, y_test = train_test_split(pd.DataFrame(x), y, test_size=0.2, random_state=i)
        gcv = gcv.fit(X_train,y_train)
        print(gcv.best_params_)
        ci.append(concordance_index_censored(y_test['Status'], y_test['Survival_in_days'], gcv.predict(X_test))[0])
        print(ci)
    df[key] = ci
    
df.to_csv("SVMs.csv")