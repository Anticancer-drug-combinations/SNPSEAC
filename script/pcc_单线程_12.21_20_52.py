import numpy as np
import pandas as pd
import re,math
import concurrent.futures

data1 = pd.read_excel(r'pcbi.1006752.s002.xlsx', sheet_name="Sheet1",index_col='Drug').drop('SMILES',axis=1)
data2 = pd.read_excel(r'583drugs39cell Zscore行标准化 synergy.xlsx')
drug_names = list(data1.index)

# 合并drugA和drugB为一列作为索引
# 顺序
temp1 = [] # 分子指纹合并值 【【】，【】，【】，】
index_1 = [] # 对应的名字列表 【【】，【】，【】，】
# 逆序
temp2 = []
index_2 = []
for index,row in data2.iterrows():
    ct = list((np.array(data1.loc[row["DrugA"],:])))
    ct.extend(list((np.array(data1.loc[row["DrugB"],:]))))
    temp1.append(ct)
    index_1.append((row["DrugA"],row["DrugB"]))
    ct = list((np.array(data1.loc[row["DrugB"],:])))
    ct.extend(list((np.array(data1.loc[row["DrugA"],:]))))
    temp2.append(ct)
    index_2.append((row["DrugB"],row["DrugA"]))
    
# 两组数据的tanimoto计算
def tanimoto(p,q):
    tep1=0
    tep2=0
    lenthp=int(len(p))
    for i in range(lenthp):
        a=p[i]
        b=q[i]
        if (a==1)|(b==1): # 并
            tep1=tep1+1
        if (a==1)&(b==1): # 交
            tep2=tep2+1
    c=round((tep2 / tep1),4)      #取值4位数
    return c

def compute_tanimoto(i, j, i_n, j_n, temp1, temp2, index_1, index_2, idx):
    tmp=[]
    a = np.array(i) #取出一分子指纹 转换列表为向量
    b = np.array(j)
    for k,m,k_n,m_n in zip(temp1,temp2,index_1,index_2): # 遍历每一个药物对的分子指纹
        c = np.array(k)
        d = np.array(m)
        tp1 = tanimoto(a,c)
        tp2 = tanimoto(a,d)
        tp3 = tanimoto(b,c)
        tp4 = tanimoto(b,d)
        tps = [tp1, tp2, tp3, tp4]
        index_list = [(i_n,k_n,tp1),(i_n,m_n,tp2),(j_n,k_n,tp3),(j_n,m_n,tp4)]
        index_location = tps.index(max(tps))
        tmp.append(max(tps))
        Tanimoto_index_dist[f"{i_n},{k_n}"] = index_list[index_location]
    return idx, tmp

temp=[]
Tanimoto_index_dist = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    # 提交任务时使用 enumerate ，并将索引值作为参数传递给任务函数，最后使用这个索引值将结果插入到 temp 列表的正确位置
    futures = [executor.submit(compute_tanimoto, i, j, i_n, j_n, temp1, temp2, index_1, index_2, idx)
               for idx, (i,j,i_n,j_n) in enumerate(zip(temp1,temp2,index_1,index_2))]
    for future in concurrent.futures.as_completed(futures):
        temp.insert(future.result()[0], future.result()[1])
Tanimoto = pd.DataFrame(temp, columns=index_1, index=index_1)

# 获取栏目数组
ttpp = []
for i in list(Tanimoto.columns):
    ttpp.append(str(i)[1:-1]) # 去除首尾括号的药物对名字
data2 = pd.DataFrame(data2.iloc[:,2:].T.values, columns=list(ttpp)) # 药物对协同得分
columns = data2.columns # 修改后的药物对名字

# 删去对应行后的字典
dict_lines = {} # {键1：值1} 谷本系数字典
ct = 0
for i in Tanimoto.index:
    y = list(Tanimoto.iloc[:,ct])
    y.pop(ct)
    dict_lines[columns[ct]] = y
    ct += 1

# 高相似多少 - 优化后不使用 dict_columns 了
coum = 200
dict_columns_split = {}
for col in columns:
    y = data2[col]
    temp = [c for c in columns if c != col]
    x = data2[temp]
    df = pd.DataFrame(dict_lines[col], index=temp).sort_values(by=0, ascending=True)
    max_cols = df.index[-coum:]
    max_values = data2[max_cols]
    for key, y_line in y.items():
        dict_columns_split[f"{col},{key}"] = {
            "y_line": y_line,
            "max_x_line": max_values.loc[key].values,
        }

# 高相似谷本系数字典 (删去对应行后的字典)
dict_high_lines = {} 
count = 0
cct = 0
for i in Tanimoto.index: 
    y = list(Tanimoto.iloc[:,cct])
    y.pop(count) 
    count += 1
    dict_high_lines[str(i)[1:-1]] = y
    tp = []
    for j in columns:
        if str(i)[1:-1] != j:
            tp.append(j)
    df = pd.DataFrame(dict_high_lines[str(i)[1:-1]], index=tp)
    df.sort_values(by=0 , inplace=True, ascending=True) 
    max_l = df.index[-coum:]  
    tp = [] 
    for j in max_l:
        tp.append(Tanimoto.iloc[list(columns).index(j), :])
    max_5y = list(pd.DataFrame(tp).T.iloc[[list(columns).index(str(i)[1:-1])]].values)[0]
    data = {
        "y": y,
        "max_relevance": max_l,
        "high_5": max_5y
    }
    dict_high_lines[str(i)[1:-1]] = data
    cct += 1
dict_high_lines_split = {}
for i in dict_high_lines: # i 是目标药物对名字
    count = 0
    for j in dict_high_lines[i]['max_relevance']: # 遍历高相似的名字
        dict_high_lines_split[f'{i},{j}']=dict_high_lines[i]['high_5'][count]
        count += 1
TTT = [] # 值列表
for i in dict_high_lines_split:
    TTT.append(dict_high_lines_split[i])


# 定义函数来计算 MSE 和 RMSE
def compute_mse_rmse(y_pred, y_true):
    y_pred = np.array(y_pred)
    y_true = np.array(y_true)
    mse = np.mean(np.square(y_pred - y_true))
    rmse = math.sqrt(mse)
    return mse, rmse

# 使用 list comprehension 和 zip 函数来计算 PCC
def compute_pcc(y_pred, y_true):
    y_pred_mean = np.mean(y_pred)
    y_true_mean = np.mean(y_true)
    numerator = sum((y_pred - y_pred_mean) * (y_true - y_true_mean))
    denominator = math.sqrt(sum(np.square(y_pred - y_pred_mean)) * sum(np.square(y_true - y_true_mean)))
    return numerator / denominator

# 定义函数来计算预测值
def compute_prediction(TTT, alp, data2, dict_columns_split, linename, coum):
    pre = []
    for elem in range(39):
        pp = []
        for j in TTT[(200 * a - coum):(200 * a)]: 
            pp.append(np.exp(-pow(1 - j, 2) / (2 * pow(alp, 2))))
        pre.append(np.sum(pp / np.sum(pp) * dict_columns_split[f"{linename},{elem}"]['max_x_line'][-coum:]))
    return pre

Name = []
a = 1
for linename in columns:
    PCC = 0
    MSE = 1000
    result = []
    Name_data = []
    for coum in range(2, 201, 1):  
        #求PCC
        alpha = []
        result_pcc = []
        lin = np.linspace(1,3,401)
        pcc_split = []
        for alp in lin:
            # 计算预测值
            y_pred = compute_prediction(TTT, alp, data2, dict_columns_split, linename, coum)
            # 计算 PCC
            pcc_split.append(compute_pcc(y_pred, data2[linename]))
        result_pcc.append(np.max(pcc_split))
        alpha.append(lin[pcc_split.index(np.max(pcc_split))])
        pcc=result_pcc[0]
        #求MSE
        alp = alpha[0]
        # 计算预测值
        y_pred = compute_prediction(TTT, alp, data2, dict_columns_split, linename, coum)
        # 计算真实值
        y_true = []
        for elem in range(39):
            y_true.append(dict_columns_split[f"{linename},{elem}"]['y_line'])
        # 计算 MSE 和 RMSE
        mse, rmse = compute_mse_rmse(y_pred, y_true)
        if pcc > PCC:
            Name_coum, Name_mse, Name_pcc, Name_alp, PCC, Name_rmse = coum, mse, pcc, alp, pcc, rmse
        print("name:{},coum:{:0d},alp={},mse={},pcc={},rmse={}".format(linename, coum, alp, mse, pcc,rmse))
        if coum == 200:
            Name_data.append(linename)
            Name_data.append(Name_coum)
            Name_data.append(Name_alp)
            Name_data.append(Name_mse)
            Name_data.append(Name_pcc)
            Name_data.append(Name_rmse)
        a = a+1
        Name.append(Name_data)
pd.DataFrame(Name).to_csv('Tanimotomax synergy anlieduqu zscore result all.csv', header=None, index=False)