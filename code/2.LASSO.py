'''---------------------------LASSO---------------------------------'''
from sklearn.linear_model import Lasso,LassoCV
import warnings
warnings.filterwarnings(action='ignore')
import pandas as pd
import numpy as np
import pandas
def cancer_LASSO(n,g):
    '''evo'''
    DATA_path = './RNA-seq/TCGA_%s_exp_classmean_p0.05_evo_quartile.txt'%(n)
    DATA = pd.read_table(filepath_or_buffer = DATA_path,low_memory=False)
    x = DATA.iloc[:,:-1]
    y=DATA.iloc[:,-1]
    X_train=x
    y_train=y
    alpha_range = np.logspace(-4, -0.8, 500, base=10)
    lasso_cv = LassoCV(alphas=alpha_range,max_iter=10000,random_state=123)
    lasso_cv.fit(X_train, y_train)
    lasso_best_alpha = lasso_cv.alpha_
    lasso001=Lasso(alpha=lasso_best_alpha).fit(X_train,y_train)
    print(n+'The number of evo genes with non-zero correlation coefficients is：',np.sum(lasso001.coef_ != 0))
    mask = lasso001.coef_ != 0
    new_reg_data = x.iloc[:, mask]
    dic = {'feature': x.columns, 'coefficient': lasso001.coef_}
    df = pd.DataFrame(dic)
    df1 = df[df['coefficient'] != 0]
    df1.to_csv('./Input features/%s_feature_mean.csv'%(n),index=0,header=True,encoding="utf-8")
    feature=pd.DataFrame()
    dict1 = dict(zip(df1['feature'],df1['coefficient']))
    for i in new_reg_data.columns:
        feature[i]=new_reg_data[i]*dict1[i]
    merge=pd.concat([feature ,y],axis=1)
    merge.to_csv('./Input features/%s_exp_best_feature_mean.csv'%(n),header=True, encoding="utf-8")
    '''nonevo'''
    DATA_path = './RNA-seq/TCGA_%s_exp_classmean_p0.05_nonevo_quartile.txt'%(n)
    DATA = pd.read_table(filepath_or_buffer = DATA_path,low_memory=False)
    result_path='./Input features/'
    x = DATA.iloc[:,:-1]
    y=DATA.iloc[:,-1]
    lasso_cv = LassoCV(alphas=alpha_range,max_iter=10000,random_state=123)
    for j in range(1,6):
        ran=x.sample(n=int(g),axis=1,random_state=j)
        X_train=ran
        y_train=y
        lasso_cv.fit(X_train, y_train)
        lasso_best_alpha = lasso_cv.alpha_
        print(lasso_best_alpha)
        lasso001=Lasso(alpha=lasso_best_alpha).fit(X_train,y_train)
        print(n+'The number of non-evo genes with non-zero correlation coefficients is：',np.sum(lasso001.coef_ != 0))
        mask = lasso001.coef_ != 0
        new_reg_data = ran.iloc[:, mask]
        dic = {'feature': ran.columns, 'coefficient': lasso001.coef_}
        df = pd.DataFrame(dic)
        df1 = df[df['coefficient'] != 0]
        df1.to_csv(result_path+'%s_non%d_feature_mean.csv'%(n,j),index=0,header=True,encoding="utf-8")
        feature=pd.DataFrame()
        dict1 = dict(zip(df1['feature'],df1['coefficient']))
        for i in new_reg_data.columns:
            feature[i]=new_reg_data[i]*dict1[i]
        merge=pd.concat([feature ,y],axis=1)
        merge.to_csv(result_path+'%s_exp_non%d_best_feature_mean.csv'%(n,j),header=True, encoding="utf-8")
cancer={"HNSC":471,"LIHC":281,"LUAD":545,"LUSC":480}
for n,g in cancer.items():
    cancer_LASSO(n,g)

d=input("")
















































'''-------------------------------画图---------------------------------'''
import matplotlib.pyplot as plt
Lambdas = np.logspace(-5, 2, 200)  # 10的-5到10的2次方
# 构造空列表，用于存储模型的偏回归系数
lasso_cofficients = []
for Lambda in Lambdas:
    lasso = Lasso(alpha=Lambda, normalize=True, max_iter=10000)
    lasso.fit(X_train, y_train)
    lasso_cofficients.append(lasso.coef_)
# 绘制Lambda与回归系数的关系
plt.plot(Lambdas, lasso_cofficients)
# 对x轴作对数变换
plt.xscale('log')
# 设置折线图x轴和y轴标签
plt.xlabel('Lambda')
plt.ylabel('Cofficients')
# 显示图形
plt.show()




import numpy as np # 快速操作结构数组的工具
import pandas as pd
import warnings
warnings.filterwarnings(action='ignore')
from sklearn.linear_model import Lasso,LassoCV

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LassoCV,lasso_path
from sklearn.datasets import load_diabetes

# 加载糖尿病数据集
X, y = load_diabetes(return_X_y=True)

# 创建LassoCV对象，指定正则化参数(alpha)的候选值
alphas_lasso, coefs_lasso, _ = lasso_path(X, y, eps=0.001,n_alphas=1000,random_state=123456)

# 训练模型
lasso.fit(X, y)

# 输出最佳正则化参数和对应的系数
print("Best alpha:", lasso.alpha_)
print("Coefficients:", lasso.coef_)

# 绘制系数的路径图
alphas = alphas_lasso
coefs = coefs_lasso
plt.figure(figsize=(10, 6))
plt.plot(np.log10(alphas), coefs.T)
plt.axvline(np.log10(lasso.alpha_), linestyle='--', color='k',
            label='alpha: CV estimate')
plt.xlabel('log(alpha)')
plt.ylabel('coefficients')
plt.title('Lasso coefficients as a function of the regularization')
plt.legend()
plt.axis('tight')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import lasso_path,LassoCV
import pandas as pd
import warnings
warnings.filterwarnings(action='ignore')
from itertools import cycle

DATA_path = r'F:\postgraduate\predict model\random forest\TCGA_CRC_exp_class3_suvival0.05_evo_quartile.txt'
DATA = pd.read_table(filepath_or_buffer = DATA_path,low_memory=False)
x = DATA.iloc[:,:-1]
#x=np.asarray(x)
y=DATA.iloc[:,-1]
from sklearn.model_selection import train_test_split
#X_train,X_test,y_train,y_test = train_test_split(x,y,random_state=1369)
X_train=x
y_train=y
lasso_cv = LassoCV(n_alphas=10000,cv=5, max_iter=10000,random_state=123456)
lasso_cv.fit(X_train, y_train)
lasso_best_alpha = lasso_cv.alpha_
lasso_alpha=lasso_cv.alphas_
lasso001=Lasso(alpha=lasso_best_alpha).fit(X_train,y_train)
mask = lasso001.coef_ != 0  #返回一个相关系数是否为零的布尔数组
# Compute paths
#alphas_lasso, coefs_lasso, _ = lasso_path(x, y, alphas=np.logspace(-3,1,50),random_state=123456)
# neg_log_alphas_enet = -np.log10(alphas_enet)
# for coef_l in coefs_lasso:
#     l1 = plt.plot(neg_log_alphas_lasso, coef_l)

alphas,coefs,_ = lasso_cv.path(x, y, alphas=np.logspace(-3,1,100),random_state=123456)
coefs=coefs.T
coefs = coefs[:,mask == True]
neg_log_alphas_lasso = np.log10(alphas)
plt.figure(1)
coefs=coefs.T
colors = cycle(["black", "Coral", "Teal", "blue", "Crimson","Olive","Indigo"])
for coef_l, c in zip(coefs, colors):
    l1 = plt.plot(neg_log_alphas_lasso, coef_l, c=c,linewidth=1.5)
#plt.plot(neg_log_alphas_lasso, coefs,linewidth='2')
plt.xlabel("Log(alpha)")
plt.ylabel("coefficients")
plt.title("Lasso Paths")
plt.axis("tight")
plt.axvline(np.log10(lasso_best_alpha), linestyle=(0,(3,1,1,1)), color='crimson',label='Best alpha',linewidth='3')




from sklearn.linear_model import LassoCV
import matplotlib.pyplot as plt
import numpy as np

# 创建一个示例数据集
X = np.random.randn(100, 10)
y = np.random.randn(100)

# 初始化LassoCV模型，设置alpha值的范围和步长
model = LassoCV(alphas=np.logspace(-4, 0, 50), cv=5)

# 拟合数据
model.fit(X, y)

# 可视化alpha值与平均误差的关系
plt.figure(figsize=(10, 6))
plt.plot(model.alphas_, model.mse_path_.mean(axis=-1))
plt.xlabel('alpha')
plt.ylabel('Average MSE')
plt.title('LassoCV Alpha Selection')
plt.legend()
plt.show()
c=model.mse_path_