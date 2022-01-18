'''
文件名：西华大学数值计算C++实验报告——实验七
内容：欧拉公式与龙格公式
当前版本：1.0
完成作者：李子涵
完成日期：2021.12.22
'''

import math
def dfx(x,y):
    return x*x+y

def Euler(x0,y0,h,N):
    """
	x0,x1 旧值
	y0y1 新值
	h 步长
	N 步数
	yp yc 平均化的中间值
	"""
    n=1
    x1=0
    while(1):  
        x1=x0+h
        x1 = x0 + h
        yp = y0 + h * dfx(x0, y0)
        yc = y0 + h * dfx(x1, yp)
        y1 = (yp + yc) / 2
        print("x"+str(n)+"   "+"y"+str(n))
        print(str(round(x1,2))+"   "+str(y1))
        if n == N :
            break
        else:
            x0 = x1
            y0 = y1
            n+=1

def Longe_kutta(x0,y0,h,N):
    """
	x0,x1 旧值
	y0y1 新值
	h 步长
	N 步数
	k1 k2 k3 k4 平均化的中间值
	"""

    n=1
    x1=0
    while(1):  
        x1=x0+h
        K1 =dfx(x0,y0);
        K2 = dfx(x0 + h / 2, y0 + h / 2*K1);
        K3 = dfx(x0 + h / 2, y0 + h / 2 * K2);
        K4 = dfx(x1, y0 + h * K3);
        y1 = y0 + h * (K1 + 2*K2 + 2*K3 + K4) / 6;
        print("x"+str(n)+"   "+"y"+str(n))
        print(str(round(x1,2))+"   "+str(y1))
        if n == N :
            break
        else:
            x0 = x1
            y0 = y1
            n+=1
if __name__ == '__main__':

    x0=0
    y0=1
    h=0.2
    N=5
    print("改进的欧拉公式:")
    Euler(x0,y0,h,N)
    print("龙格公式:")
    Longe_kutta(x0,y0,h,N)
    