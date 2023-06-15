import numpy as np
import pandas as pd
import re,math
import multiprocessing
import multiprocessing.pool
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

temp = [None] * len(temp1)
Tanimoto_index_dist = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    # 提交任务时使用 enumerate ，并将索引值作为参数传递给任务函数
    futures = [executor.submit(compute_tanimoto, i, j, i_n, j_n, temp1, temp2, index_1, index_2, idx)
               for idx, (i,j,i_n,j_n) in enumerate(zip(temp1,temp2,index_1,index_2))]
    for future in concurrent.futures.as_completed(futures):
        idx, result = future.result()
        temp[idx] = result
Tanimoto = pd.DataFrame(temp, columns=index_1, index=index_1)

# 获取columns数组 - 处理后的
ttpp = []
for i in list(Tanimoto.columns):
    ttpp.append(str(i)[1:-1]) # 去除首尾括号的药物对名字
data2 = pd.DataFrame(data2.iloc[:,2:].T.values, columns=list(ttpp)) # 药物对协同得分
columns = data2.columns # 修改后的药物对名字

# 删去对应行后的字典
dict_lines = {} # {键1：值1} 谷本系数字典
count = 0
ct = 0
for i in Tanimoto.index:
    y = list(Tanimoto.iloc[:,ct]) # 每一列谷本系数 [:,1] 第一列的每一行 实际就是第一列
    y.pop(count)
    count += 1
    dict_lines[columns[ct]] = y # 谷本系数药物对的名字做键1，y是值1
    ct += 1
    
coum=200
# 分离x y数据构造字典
dict_columns_split = {} # 高相似字典
# 计算x数据
x = data2.drop(columns=columns)  # 其他所有的列
for i in columns:
    # b 要预测协同得分那一列
    y = data2.loc[:, i]
    # A 权重
    x_i = x.copy()
    x_i['y'] = y

    # 高相似的    取高相似的药物对名字max1和对应的值max5
    tp = []
    for j in columns:
        if i != j:
            tp.append(j)
    df = pd.DataFrame(dict_lines[i], index = tp)
    df.sort_values(by=0, inplace=True, ascending=True)
    max_l = df.index[-coum:]
    max_5x = pd.DataFrame([data2.loc[:, j] for j in max_l]).T

    # 存储 对应数据
    count = 0 # 第几行 int
    for j in list(y.index): # 35个索引 0~34 数字（字符串-str） 【0，1，2，3.。。。】
        dict_columns_split[f'{i},{j}']={} # 空字典 键-空字典
        dict_columns_split[f'{i},{j}']['y_line']=y[count]
        dict_columns_split[f'{i},{j}']['max_x_line']=max_5x.iloc[count,:].values
        count += 1

TTT = [] # 值列表
cct = 0
for i in Tanimoto.index: 
    y = list(Tanimoto.iloc[:,cct])
    y.pop(cct) 
    tp_columns = list(columns.copy())
    tp_columns.pop(cct)
    df = pd.DataFrame(y, index=tp_columns)
    df.sort_values(by=0 , inplace=True, ascending=True) 
    max_l = df.index[-coum:]  
    # 计算max_5y
    max_5y = []
    for j in max_l:
        max_5y.append(Tanimoto.iloc[list(columns).index(j), cct])
    # 将max_5y加入TTT
    TTT += max_5y
    # 修改cct变量的值
    cct += 1

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
            pre = []
            for elem in range(39):
                pp = []
                for j in TTT[(200 * a - coum):(200 * a)]: 
                    pp.append(np.exp(-pow(1 - j, 2) / (2 * pow(alp, 2))))#改函数
                pre.append(np.sum(pp / np.sum(pp) * dict_columns_split[f"{linename},{elem}"]['max_x_line'][-coum:]))
            pcc_split.append(np.corrcoef(pre, data2[linename])[0, 1])
        result_pcc.append(np.max(pcc_split))
        alpha.append(lin[pcc_split.index(np.max(pcc_split))])
        pcc=result_pcc[0]
        #求MSE
        r1 = []
        ct = 0
        alp = alpha[0]
        for elem in range(39):
            pp = []
            for j in TTT[(200 * a - coum):(200 * a)]: # dict_high_lines_split 里一组的数量是80
                pp.append(np.exp(-pow(1 - j, 2) / (2 * pow(alp, 2))))  # 改函数
            r1.append(np.sum((pp/np.sum(pp)) * dict_columns_split[f"{linename},{elem}"]['max_x_line'][-coum:]) -
                            dict_columns_split[f"{linename},{elem}"]['y_line'])
        mse=np.mean(np.square(r1))
        rmse=math.sqrt(mse)
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