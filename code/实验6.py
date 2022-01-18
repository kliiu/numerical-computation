'''
文件名：西华大学数值计算C++实验报告——实验六
内容：复化梯形及复化辛普森
当前版本：1.0
完成作者：李子涵
学号：3120190971401
完成日期：2021.11.20
'''
import math
#梯形公式
def tixing(f,a,b,n,h):
    x=a
    s=f(x)+f(b)
    for k in range(0,n-1):
        x=x+h
        s=s+2*f(x)
    result=(h/2)*s
    return result

#辛普森公式
def simprson(f,a,b,m,h):
    """ f被积函数，a,b分别为积分的上下线，m为子区间的个数，s是梯形总面积"""
    x=a
    h = (b-a)/m
    s=f(x)+f(b)
    x=x+h/2
    for i in range (0,int(m)):
        s=s+4*f(x)
        x=x+h
    x=a+h
    for i in range (1,int(m)):
        s=s+2*f(x)
        x=x+h
    result = (h/6)*s
    return result

if __name__ == '__main__':
    a=dict({0.6:5.7,0.8:4.6,1:3.5,1.2:3.7,1.4:4.9,1.6:5.2,1.8:5.5})
    h=0.2
    x=list(a.keys())
    n=int((x[-1]-x[0])/h)#区间数
    f=lambda x :a[round(x,2)]
    print("复化梯形：")
    print(tixing(f, x[0], x[-1], n, h))
    print("复化Simpson:")
    print(simprson(f, x[0], x[-1], n/2, h))
