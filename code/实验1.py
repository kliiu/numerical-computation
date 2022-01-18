'''
文件名：西华大学数值计算C++实验报告——实验一
当前版本：1.0
完成作者：李子涵
学号：3120190971401
完成日期：2021.11.10
'''

#Lagrange插值多项式
def Lang(val):#输入参数：预测点
    x=[pairs[i] for i in range(len(pairs)) if i%2==0]#x存放插值节点的x值
    y=[pairs[i] for i in range(len(pairs)) if i%2!=0]#y存放插值节点的y值
    ans = 0
    for i in range(n):
        temp=1#拉格朗日基函数
        for j in range(i-1):
            if j!=i:
                temp=temp*(val-x[j])/(x[i]-x[j])
        for j in range(i-1,n):
            if j>=0 and j!=i:
                temp=temp*(val-x[j])/(x[i]-x[j])

        ans=ans+temp*y[i]#计算插值多项式
    print("拉格朗日插值多项式：")
    print('f({}) = {} (保留三位小数)'.format(val, format(ans, '.3f')))
    #print(ans)

#牛顿插值多项式
def Newton(val):#输入参数：预测点
   #构造差商表，每一列表示一组
    x=[pairs[i] for i in range(len(pairs)) if i%2==0]#x存放插值节点的x值
    y=[pairs[i] for i in range(len(pairs)) if i%2!=0]#y存放插值节点的y值
    i = 1
    g=[pairs[i] for i in range(len(pairs)) if i%2!=0]# g[k]表示f[0,1,2,...,k],首先将y值赋给g
    #构造差商表
    while(i<n):
        k=0
        j=i
        while(j<n):
            g[j]=(y[j]-y[j-1])/(x[j]-x[k])
            j+=1
            k+=1
        y=[g[i] for i in range(len(g))]#将此列差商存放进y
        i+=1
    t=1
    i=1
    ans=g[0]
    while(i<n):#计算插值多项式
        t=t*(val-x[i-1])
        ans=ans+t*g[i]
        i+=1
    print("牛顿插值多项式：")
    print('f({}) = {} (保留三位小数)'.format(val, format(ans, '.3f')))
    #print(ans)

if __name__=='__main__':
    pairs=[1920,105711,1930,123203,1940,131669,1950,150697,1960,179323]#将每个点的x,y值成对列出 eg:(-2,17)(0,1)(1,2)(2,19)
    if(len(pairs)%2!=0):#若有奇数个点则给出提示
        print("missing element")
    n=int(len(pairs)/2)#n表示插值点数量，插值多项式阶数
    Lang(1910)
    
