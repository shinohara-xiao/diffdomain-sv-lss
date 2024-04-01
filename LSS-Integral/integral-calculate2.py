import numpy as np
from scipy import integrate
from scipy.integrate import dblquad #dblquad用于二重积分
#--------------------------------------------------------------#
# 一重积分
from scipy import integrate
def f(x, a, b):
    return a * x + b
v, err = integrate.quad(f, 1, 2, args = (-1, 1))  # 积分上下限分别为1，2。对函数 y=-x+b 求积分，积分结果为v
print(v)

#--------------------------------------------------------------#
# 重积分方法 一
# def main():
#    print(dblquad(lambda t,x:np.sin(t)*np.exp(-x*t)/t**5,0.5,0.8,lambda x:0.2,lambda x:0.7))
    #被积函数是sint*exp(-xt)/t^5, 其中 t 的积分上下限是 0.5和0.8，x 的积分上下限是0.2和0.7     
    
# 调用函数
# if __name__ == "__main__":
#    main()

# 重积分方法
from scipy import integrate
import numpy as np
def f(x, y):
    return x * y
def h(x):
    return x
v, err = integrate.dblquad(f, 1, 2, lambda x: 1, h) # 从 1 到 h
print(v)


# 重积分方法 二
from scipy import integrate
import numpy as np
f = lambda x, y, z : x
g = lambda x : 0
h = lambda x : (1 - x) / 2
q = lambda x, y : 0
r = lambda x, y : 1 - x - 2 * y 
v, err = integrate.tplquad(f, 0, 1, g, h, q, r)  # g,h,q,r 分别是一重和二重积分的上下限
print(v)


# 重积分方法 三

# 引入需要的包
import scipy.integrate
from numpy import exp
from math import sqrt
import math

# 创建表达式
f = lambda x,y : exp(x**2-y**2)

# 计算二重积分：（p:积分值，err:误差）
# 这里注意积分区间的顺序
# 第二重积分的区间参数要以函数的形式传入
p,err= scipy.integrate.dblquad(f, 0, 2, lambda g : 0, lambda h : 1)	
print(p)

# 重积分方法 四
from scipy.integrate import dblquad
area = dblquad(lambda x, y: x*y, 0, 0.5, 
               lambda x: 0, lambda x: 1-2*x)
print(area)

# 重积分方法 五         使用integrate.nquad
from scipy import integrate
def f(x, y):
    return x*y

def bounds_y():
    return [0, 0.5]

def bounds_x(y):
    return [0, 1-2*y]

print(integrate.nquad(f, [bounds_x, bounds_y]))
