'''
文件名：西华大学数值计算C++实验报告——实验五
内容：求解线性方程组的Jaccobi迭代和Gauss-Seidel迭代法
当前版本：1.0
完成作者：李子涵
学号：3120190971401
完成日期：2021.11.20
'''
import numpy as np
import math
from numpy.core.fromnumeric import transpose
from numpy.matrixlib.defmatrix import matrix


def Jacobi(matrix_iter,X_for,X_next,epsilon):
    """
    Jacobi方法求解方程组
    输入：系数矩阵,X_0,X_1向量,ε误差下界
    """
    shape=matrix_iter.shape
    while(max(abs(X_next-X_for)) > epsilon): #直到结果与上一代逐项相减，最大项不大于epsilon  
        X_for=X_next
        X_next=matrix_iter@X_for#更新下一代
        X_next=np.insert(X_next, shape[0], values=1, axis=0)#加一行1
    print("Jacobi求解结果为：") 
    print(X_next[:shape[0]])
def Gauss_Seidel(coef,b,matrix_iter,X_next,epsilon):
    """
    Gauss_Seidel方法求解方程组
    输入：系数矩阵,X_0,X_1向量,ε误差下界
    """
    shape=matrix_iter.shape
    X_for=X_next*100000+epsilon#初始化为一个很大的数
    while(max(abs(coef@X_next[:-1]-matrix_Y)) >=epsilon): #直到结果与上一代逐项相减，最大项不大于epsilon  
        X_for=X_next
        for i in range(shape[0]):#一共有i次计算
            X_next[i] = matrix_iter[i]@X_for
            X_for[i]=X_next[i]
    print("Gauss_Seidel求解结果为：") 
    print(X_next[:shape[0]])

if __name__=='__main__':
    #给出系数矩阵
    coef=np.matrix([(2,-1,-1),(1,5,-1),(1,1,10)])
    shape=coef.shape#系数矩阵形状
    b=[-5,8,11]#给定y值
    epsilon=10e-4#给定epsilon
    X_0=np.matrix([1.,1.,1.]).T#给定初始值
    #初始值多加一行1，使常数项始终为1
    X_0=np.insert(X_0, shape[0], values=1., axis=0)
    #开始求解
    matrix_iter=np.zeros(shape)
    y=np.zeros((shape[0],1))#存放y值向量
    for i in range(shape[0]):
        y[i]=b[i]/coef[i,i]#常数项除以a_ii
        matrix_iter[i]=coef[i,:]/coef[i,i]*(-1)#将系数矩阵每一行除以a_ii
        matrix_iter[i,i]=0#a_ii变为0
    matrix_iter=np.concatenate((matrix_iter,y),axis=1)#将系数矩阵与常数项合并
    X_1=matrix_iter@X_0#下一代=系数矩阵与初代向量作矩阵乘法（3，4）（4，1）
    #将得到的（3，1）矩阵插入1，变为（4，1）继续做矩阵乘法
    X_1=np.insert(X_1, shape[0], values=1, axis=0)

    Jacobi(matrix_iter,X_0,X_1,epsilon)
    #将b转为列向量
    matrix_Y=np.zeros(y.shape)
    for i in range(matrix_Y.shape[0]):
        matrix_Y[i]=b[i]
    Gauss_Seidel(coef,matrix_Y.T,matrix_iter,X_0,epsilon)
