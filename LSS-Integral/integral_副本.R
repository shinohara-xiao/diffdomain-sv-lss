# 加载包
#install.packages('cubature')
library(cubature)
library(e1071)
# ------------- 函数定义 ----------------- #

# 计算均值
mean_carlo<-function(g) # con = 1 -- x ; con = 2 -- x^2 ; con = 3 --x^3
{
  # 计算积分
  #ff <- function(x){g(2*x)*((pi)^-1)/sqrt(1-x^2)}
  ff <- function(x) {g(2*x) /sqrt(1-x^2)}
  m <- 0.25*(g(-2) + g(2)) - 0.5 * pi^-1 * integrate(ff,-1,1)$value
  #cat('mean is:',m,'\n')  # -- 0
  return(m)
}

# 计算方差
sim_carlo<-function(g,x1,x2,y1,y2,n,k)
{
  arr = numeric()
  v = (x2-x1)*(y2-y1)
  for(i in 1:k)
  {
    r1 = runif(n,min = x1,max = x2) 
    r2 = runif(n,min = y1,max = y2) 
    h = g(r1,r2)
    out = v*mean(h)
    arr[i] = out
  }
  #cat('var is:',mean(arr))
  return(mean(arr))
}


# ------ 函数调用 ------ #

# 定义方差函数
# x,y
co1 = function(x,y) {((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}           # ----- 2
# 取x^2,y^2
co2 = function(x,y) {2*x*2*y*((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}      # ------4
# 取x^3,y^3
co3 = function(x,y) {3*x^2*3*y^2*((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}  # ------24
# 取 e^x,e^y
co4 = function(x,y) {exp(x+y)*((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}  # ------ 7.25
# 取 ln（1+x^2）
co5 = function(x,y) {(4*x*y/((1+x^2)(1+y^2)))*((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}

# 定义均值函数
f1 = function(x) {x}
f2 = function(x) {x^2}
f3 = function(x) {x^3}
f4 = function(x) {exp(x)}
f5 = function(x) {log(1+x^2)}

n = 100
k = 100

## 
c(mean_carlo(f1), sim_carlo(g=co1))
ls(n=n, k=k, g=f1)

c(mean_carlo(f2), sim_carlo(g=co2))
ls(n=n, k=k, g=f2)

c(mean_carlo(f3), sim_carlo(g=co3))
ls(n=n, k=k, g=f3)

c(mean_carlo(f4), sim_carlo(g=co4))
ls(n=n, k=k, g=f4)

c(mean_carlo(f5), sim_carlo(g=co5))
ls(n=n, k=k ,g=f5)

# integral test --二重积分验证 
t1 = function(x,y) { x + y^2 }                 # ----- 21.333
t2 = function(x,y) { 2*x + 3*y }               # ----- 0
t3 = function(x,y) { exp(x) + exp(y) }         # ----- 8*(exp(2)-exp(-2))   
t4 = function(x,y) { x^2 + y^2 }               # ----- 42.67
t5 = function(x,y) { x^2 + y^3 }               # ----- 21.33






