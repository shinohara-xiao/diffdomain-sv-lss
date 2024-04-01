import numpy as np
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
v, err = integrate.dblquad(f, 1, 2, lambda x: 1, h)
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
