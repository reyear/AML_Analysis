import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.model_selection import ShuffleSplit, GridSearchCV,train_test_split

from sksurv.datasets import load_veterans_lung_cancer
from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastSurvivalSVM,FastKernelSurvivalSVM
from sksurv.kernels import clinical_kernel
def score_survival_model(model, X, y):
    prediction = model.predict(X)
    result = concordance_index_censored(y['Status'], y['Survival_in_days'], prediction)
    return result[0]
### new columns indexes for the dataset with final components for analysis for features containing the component
eln_clin_demo_comp = [0]+list(range(153,180))
eln_cyto_gen_comp = list(range(171))
eln_cyto_comp = [0] + list(range(84,171))
eln_gen_comp = list(range(84)) + list(range(153,171))
clin_demo_comp = list(range(153,180))
cyto_gen_comp = list(range(1,171))
cyto_comp = list(range(84,171))
gen_comp = list(range(1,84))+list(range(153,171))
clin_demo_cyto_gen_comp = list(range(1,180))
comp = list(range(153,171))
df_final = pd.read_table("df_prognosis_features_ready_final_component.tsv")
estimator = FastSurvivalSVM(max_iter=1000, tol=1e-6, random_state=17)
param_grid = {'alpha': 10. ** np.array([-6,-5.5,-5,-4.5,-2.5,-1,0]),'optimizer':["avltree"]}
cv = ShuffleSplit(n_splits=5,random_state=17)
gcv = GridSearchCV(estimator, param_grid, scoring=score_survival_model,
                   n_jobs=4, iid=False, refit=True,
                   cv=cv)
dict_features_type_final_comp = dict(zip(("eln_clin_demo_comp_final_comp","eln_cyto_gen_comp_final_comp","eln_cyto_comp_final_comp","eln_gen_comp_final_comp",
                                          "clin_demo_comp_final_comp","cyto_gen_comp_final_comp","cyto_comp_final_comp","gen_comp_final_comp",
                                          "clin_demo_cyto_gen_comp_final_comp","comp_final_comp"),
                                         (eln_clin_demo_comp,eln_cyto_gen_comp,eln_cyto_comp,eln_gen_comp,clin_demo_comp,cyto_gen_comp,
                                          cyto_comp,gen_comp,clin_demo_cyto_gen_comp,comp)))
df=pd.DataFrame(columns=dict_features_type_final_comp.keys())

for key,item in dict_features_type_final_comp.items():
    x = df_final.iloc[:,item]
    y = np.array(list(zip(df_final.os_status, df_final.os)),dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
    ci=[]
    for i in range(25):
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=i)
        gcv = gcv.fit(X_train,y_train)
        print(gcv.best_params_)
        ci.append(concordance_index_censored(y_test['Status'], y_test['Survival_in_days'], gcv.predict(X_test))[0])
        print(ci)
    df[key] = ci
    
df.to_csv("SVM_final_comp.csv")