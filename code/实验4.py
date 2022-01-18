'''
文件名：西华大学数值计算C++实验报告——实验四
内容：求解线性方程组的Dolittle分解法和Crout分解法
当前版本：1.0
完成作者：李子涵
学号：3120190971401
完成日期：2021.11.20
'''
import numpy as np
from numpy.core.fromnumeric import transpose
import math

from numpy.matrixlib.defmatrix import matrix
def Doolittle(coef,b):
    """
    用Doolittle求解方程式(系数矩阵为方阵)
    """
    def add_1(i,j):
        """
        传入第i行第j列，求解由r=0到i-1求和的值：l_kr*u_rk
        """
        sum=0
        for r in range(i):
            sum+=matrix_L[i,r]*matrix_U[r,j]
        return sum
    def add_2(i):
        """
        求和
        """
        sum=0
        for j in range(i):
            sum+=matrix_L[i,j]*y[j]
        return sum
    def add_3(i):
        """
        求解求和从j=i+1到shape :u_ij*x_j
        """
        sum=0
        for j in range(i,shape):
            sum+=matrix_U[i,j]*x[j]
        return sum
    
    matrix_a=coef#系数矩阵
    matrix_U=np.zeros(matrix_a.shape)#U矩阵
    matrix_L=np.eye(matrix_a.shape[0])#L矩阵
    if(matrix_a.shape[0]!=matrix_a.shape[1]):
        print("非方阵")
        return#若不为方阵打印信息并退出
    shape=matrix_a.shape[0]#shape为方阵的行（列）数
    
    for i in range(shape):
        #求解U中第i行
        for j in range(i,shape):
            matrix_U[i,j]=matrix_a[i,j]-add_1(i,j)
        #求解L中第i列
        for k in range(i+1,shape):
            matrix_L[k,i]=(matrix_a[k,i]-add_1(i,j))/matrix_U[i,i]
    
    #求解方程组LY=b
    y=[0]*shape#创建y初始值为0
    for i in range(shape):
        y[i]=b[i]-add_2(i)#第二个求和

    #求解方程组UX=Yx
    x=[0]*shape#创建x初始值为0
    for i in range(shape-1,-1,-1):#range最后一个参数表示倒着取值
        x[i]=(y[i]-add_3(i))/matrix_U[i,i]
    print("Doolittle求解方程答案为：")
    print(x)

def Crout(coef,b):
    """
    用Crout法求解线性方程组（系数矩阵为方阵）
    """
    def add_1(i,k):
        """
        计算r=0到k-1对L_ik*U_rk求和
        """
        sum=0
        for r in range(k):
            sum+=matrix_L[i,r]*matrix_U[r,k]
        return sum
    def add_2(j,k):
        """
        计算r=0到k-1对L_kr*U_rj求和
        """
        sum=0
        for r in range(k):
            sum+=matrix_L[k,r]*matrix_U[r,j]
        return sum
    def add_3(i):
        """
        求和
        """
        sum=0
        for j in range(i+1,shape):
            sum+=matrix_U[i,j]*x[j]
        return sum
    def add_4(i):
        """
        求解求和从j=i+1到shape :u_ij*x_j
        """
        sum=0
        for j in range(i):
            sum+=matrix_L[i,j]*y[j]
        return sum
    
    matrix_a=coef#系数矩阵
    matrix_U=np.eye(matrix_a.shape[0])#U矩阵(单位上三角矩阵)
    matrix_L=np.zeros(matrix_a.shape)#L矩阵(下三角)
    if(matrix_a.shape[0]!=matrix_a.shape[1]):
        print("非方阵")
        return#若不为方阵打印信息并退出
    shape=matrix_a.shape[0]#shape为方阵的行（列）数
    
    for k in range(shape):#计算L与U
        for i in range(k,shape):#先算L中第K列
            matrix_L[i,k]=matrix_a[i,k]-add_1(i,k)
        for j in range(k+1,shape):#再算U中第K行
            matrix_U[k,j]=(matrix_a[k,j]-add_2(j,k))/matrix_L[k,k]
    
    #求解方程组LY=b
    y=[0]*shape#创建y初始值为0
    for i in range(shape):
        y[i]=(b[i]-add_4(i))/matrix_L[i,i]
    #求解方程组UX=Yx
    x=[0]*shape#创建x初始值为0
    for i in range(shape-1,-1,-1):#range最后一个参数表示倒着取值
        x[i]=y[i]-add_3(i)
    print("Crout求解方程答案为：")
    print(x)
if __name__=='__main__':
    #给出系数矩阵
    coef=np.matrix([(2,1,1),(1,3,2),(1,2,2)])
    y=[4,6,5]
    Doolittle(coef,y)
    Crout(coef,y)
