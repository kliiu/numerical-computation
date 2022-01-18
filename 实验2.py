'''
文件名：西华大学数值计算C++实验报告——实验二
当前版本：1.0
完成作者：李子涵
完成日期：2021.11.10
注：所有系数列表均为升幂排列
'''
import numpy as np
from numpy.core.fromnumeric import transpose
import math
import seaborn
def solve(coefficient):
    """
    解三元（二元）一次方程组
    1x + 2y + 3z = 4
    5x + 7y + 6z = 7
    9x + 3y + 6z = 9
    示例：
    coefficient=[[1, 2, 3, 4], [5, 7, 6, 7], [9, 3, 6, 9]]
    """
    #coe=[[5,-0.5,1.875,9.32],[-0.5,1.875,-0.6875,2.2925],[1.875,-0.6875,1.3828,2.7138]]#系数列表
    #coe=[2,-1,1,10,3,2,-1,16,1,6,-1,28]
    
    #分别表示每一个式子的系数
    flag=len(coefficient)#二元or三元标志
    first=coefficient[0]
    second=coefficient[1]
    if(flag==3):
        last=coefficient[2]

    def  elimination(first_coe,second_coe):
        """
        消去一项
        输入参数：第一个式子和第二个式子的系数
        返回参数：消去第一个未知数后的一个式子的系数
        """
        new=[]
        #若第二个式子无此未知数则直接用第二个式子的系数
        if(second_coe[0]==0):
            second_coe.pop(0)
            new=second_coe
        #若有就先消去未知数
        else:
            cnt=0
            for i in second_coe:
                new.append(first_coe[cnt]-i*first_coe[0]/second_coe[0])
                cnt+=1
            new.pop(0)#去掉0项
        return new

    #对第一二个式子进行运算，消去x
    new_first=elimination(first,second)#用1，2式消去第一个未知数得到的一个式子
    
    if(flag==3):
        #对第一三个式子进行运算，消去x
        new_second=elimination(first,last)#消元x多项式后的第二个式子

        #对新得到的两个式子消去y
        new_last=elimination(new_first,new_second)#消去y后的式子系数

    #三元代入计算z,y,x
    if(flag==3):
        z=new_last[-1]/new_last[-2]
        y=(new_first[-1]-new_first[-2]*z)/new_first[0]#由前式计算y
        x=(first[-1]-first[-2]*z-first[-3]*y)/first[0]#由前式计算x
        solution=[x,y,z]
    #二元代入计算y,x
    if(flag==2):
        y=new_first[-1]/new_first[-2]
        x=(first[-1]-first[-2]*y)/first[0]
        solution=[x,y]
    return solution

def Predict(ans_list,x):
    """
    由拟合函数预测第x年的人口数
    输入参数：
    ans_list:拟合函数的系数
    """
    ans=0
    for i in ans_list:
        ans=ans+i*math.pow(x,ans_list.index(i))
    return ans
def Quadratic(x,y,n):
    """
    二次多项式拟合
    """
    Y=y.T#将Y转为列向量
    transposed_A=np.matrix([np.ones(n),x,x*x])#构造矩阵A(转置)
    matrix_A=transposed_A.T#转置得到A
    matrix_ATA=transposed_A@matrix_A#矩阵乘法AAT，法方程的左边系数
    matrix_ATY=transposed_A@Y#法方程右边系数
    coe=matrix_ATA.tolist()#将左边系数矩阵的未知数部分转为矩阵
    
    #将右边系数添加至对应的左边系数行中，构成方程的系数矩阵
    cnt=0
    for i in coe:
        i.append(matrix_ATY[cnt,0])
        cnt+=1
    so=solve(coe)#解三元一次方程
    print("二次多项式拟合结果为:\n"+str(so[0]),end="")
    if(so[1]>0):
        print("+",end="")
    print(str(so[1])+"x",end="")
    if(so[2]>0):
        print("+",end="")
    print(str(so[2])+"x²")
    pre=Predict(so,1965.)
    print("1965年的人口数约为："+str(int(pre)))
    pre=Predict(so,2002.)
    print("2002年的人口数约为："+str(int(pre)))
def Linear(x,y,n):
    """
    线性拟合
    """
    Y=y.T#将Y转为列向量
    transposed_A=np.matrix([np.ones(n),x])#构造矩阵A(转置)
    matrix_A=transposed_A.T#转置得到A
    matrix_ATA=transposed_A@matrix_A#矩阵乘法AAT，法方程的左边系数
    matrix_ATY=transposed_A@Y#法方程右边系数
    coe=matrix_ATA.tolist()#将左边系数矩阵的未知数部分转为矩阵
    
    #将右边系数添加至对应的左边系数行中，构成方程的系数矩阵
    cnt=0
    for i in coe:
        i.append(matrix_ATY[cnt,0])
        cnt+=1
    so=solve(coe)#解二元一次方程
    print("线性拟合结果为:\n"+str(so[0]),end="")
    if(so[1]>0):
        print(" + ",end="")
    print(str(so[1])+"x")
    pre=Predict(so,1965)
    print("1965年的人口数约为："+str(int(pre)))
    pre=Predict(so,2012)
    print("2012年的人口数约为："+str(int(pre)))
def Exponent(x,y,n):
    """
    指数拟合
    """
    #将y值转为ln(y)
    for i in range(n):
        y[0,i]=math.log(y[0,i])
        a=1
    Y=y.T#将Y转为列向量
    transposed_A=np.matrix([np.ones(n),x])#构造矩阵A(转置)
    matrix_A=transposed_A.T#转置得到A
    matrix_ATA=transposed_A@matrix_A#矩阵乘法AAT，法方程的左边系数
    matrix_ATY=transposed_A@Y#法方程右边系数
    coe=matrix_ATA.tolist()#将左边系数矩阵的未知数部分转为矩阵
    
    #将右边系数添加至对应的左边系数行中，构成方程的系数矩阵
    cnt=0
    for i in coe:
        i.append(matrix_ATY[cnt,0])
        cnt+=1
    so=solve(coe)#解二元一次方程

    print("指数增长拟合结果为:\ne^"+str(so[0]),end="")
    if(so[1]>0):
        print(" + ",end="")
    print(str(so[1])+"x")
    pre=Predict(so,1965)
    print("1965年的人口数约为："+str(int(math.exp(pre))))
    pre=Predict(so,2012)
    print("2012年的人口数约为："+str(int(math.exp(pre))))

if __name__=='__main__':
    pairs=[1920.,105711.,1930.,123203.,1940.,131669.,1950.,150697.,1960.,179323.,1970.,203212.]#依次列出各点的x,y值
    #pairs=[-1,0.22,-0.5,0.8,0,2.,0.25,2.5,0.75,3.8]
    #pairs=[-3,4,-2,2,-1,3,0,0,1,-1,2,-2,3,-5]
    #pairs=[-1.15,0.22,-0.5,0.8,0.1,2,0.25,2.5,0.75,3.8]
    if(len(pairs)%2!=0):
        print("输入有误")
        exit()
    n=int(len(pairs)/2)
    x=np.array([pairs[i] for i in range(len(pairs)) if i%2==0])
    y=np.matrix([pairs[i] for i in range(len(pairs)) if i%2==0])
    Quadratic(x,y,n)
    Linear(x,y,n)
    Exponent(x,y,n)


