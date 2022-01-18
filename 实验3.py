'''
文件名：西华大学数值计算C++实验报告——实验三
当前版本：1.0
完成作者：李子涵
完成日期：2021.11.11
'''
import math
MAXREPT=1000000
def f(x):
    """
    函数表达式
    """
    y=0
    for i in coef.keys():
        y+=pow(x,i)*coef[i]
    return y
def f_d(x):
    """
    求导数
    """
    y=0
    for i in coef.keys():
        if(i!=0):
            y+=pow(x,i-1)*coef[i]*i
    return y


def Newton_iter(x_k0,epsilon):
    def g(x):
        y=x-f(x)/f_d(x)
        return y
    for i in range(1,MAXREPT):#最大迭代数
        x_k1=g(x_k0)
        if(abs(x_k1-x_k0)<epsilon):
            print("牛顿迭代法近似解为:",end="")
            print(x_k1)
            return
        x_k0=x_k1
    print("牛顿迭代法在x0附近无根")
def String(x_k0,x_k1,epsilon):
    """
    弦截法求根
    """
    #初始值
    x_k=x_k0#x_k
    x_k_min1=x_k1#x_k-1
    for i in range(1,MAXREPT):#最大迭代数
        x_k_add1=x_k-f(x_k)*(x_k-x_k_min1)/(f(x_k)-f(x_k_min1))
        if(abs(x_k_add1-x_k)<epsilon or abs(f(x_k_add1)<epsilon)):#当精度达到时，输出解
            print("弦截法近似解为:",end="")
            print(x_k_add1)
            return
        x_k_min1=x_k#否则继续循环，更新当前xk和xk-1的值
        x_k=x_k_add1
    print("弦截法在x0附近无根")

if __name__=='__main__':
    #多项式格式为系数在x前，次数在x后,且已合并S所有同类项
    coef=dict(zip([3,2,1,0],[1,-7.7,19.2,-15.3]))
    Newton_iter(1,pow(10,-9))
    String(1.5,4.0,pow(10,-9))
