# 加载包
#install.packages('cubature')
library(cubature)
library(e1071)

# ------------------------- f(x)以及co(x,y)函数定义使用 ----------------------- #

# 定义方差函数
# x,y
co1 = function(x,y) {((2*pi^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}           # ----- 2
# 取x^2,y^2
co2 = function(x,y) {((2*pi^2)^-1)*2*x*2*y*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}      # ------4
# 取x^3,y^3
co3 = function(x,y) {((2*pi^2)^-1)*3*x^2*3*y^2*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}  # ------24
# 取 e^x,e^y
co4 = function(x,y) {((2*pi^2)^-1)*exp(x+y)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))}  # ------ 7.25
# 取 ln（1+x^2）
co5 = function(x,y) {((2*pi^2)^-1)*4*x*y*((1+x^2)^-1)*((1+y^2)^-1)*log((4-x*y+sqrt(4-x^2)*sqrt(4-y^2))/(4-x*y-sqrt(4-x^2)*sqrt(4-y^2)))} # ------0.63

# 定义均值函数
f1 = function(x) {x}           # ------ 0
f2 = function(x) {x^2}         # ------ 1
f3 = function(x) {x^3}         # ------ 0
f4 = function(x) {exp(x)}      # ------ 0.74
f5 = function(x) {log(1+x^2)}  # ------ 0.32


# -------------------------------- 使用函数定义 ------------------------------- #

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
sim_carlo<-function(co, n=10000, k=1000)  # g 根据实际选择选择co1，co2等
{
  x1 = -2; x2 = 2
  y1 = -2; y2 = 2
  arr = numeric()
  v = (x2-x1)*(y2-y1)
  for(i in 1:k)
  {
    r1 = runif(n,min = x1,max = x2) # 生成均匀分布的随机数 x(r1)
    r2 = runif(n,min = y1,max = y2) # 生成均匀分布的随机数 y(r2)
    h = co(r1,r2)
    out = v*mean(h)
    arr[i] = out
  }
  #cat('var is:',mean(arr))
  return(mean(arr))
}

# 模拟生成矩阵验证计算均值/方差 
# 生成n*n 阶矩阵； 循环k次，检验k次, g是所选择的函数f(x)
ls <- function(n,k,g)       
{
  Ff <- function(x) { g(x)/2/pi*sqrt(4-x^2) }
  Gn.F = integrate(Ff,-2,2,subdivisions = 10000)$value
  Gn.F = Gn.F * n
  
  Gn.f<-c()
  for(i in 1:k)
  {
    mat <- matrix(rnorm(n^2,mean=0, sd=1),n,n)
    mat[lower.tri(mat)] = 0
    mat <- mat+t(mat)-2*diag(diag(mat))+diag(rnorm(n,mean=0, sd=sqrt(2)))
    mat <- mat/sqrt(n)  # 矩阵标准化
    a <- eigen(mat)$values # 提取特征值
    Gn.Fn = sum(g(a))
    Gn.f[i] = Gn.Fn - Gn.F 
  }
  # return(Gn)
  cat('mean:',mean(Gn.f),'var:',var(Gn.f))
}



# --------- integral test --二重积分验证 --------- #

t1 = function(x,y) { x + y^2 }                 # ----- 21.333
t2 = function(x,y) { 2*x + 3*y }               # ----- 0
t3 = function(x,y) { exp(x) + exp(y) }         # ----- 8*(exp(2)-exp(-2)) , 58.02977
t4 = function(x,y) { x^2 + y^2 }               # ----- 42.67
t5 = function(x,y) { x^2 + y^3 }               # ----- 21.33

sim_carlo(co=t1)
sim_carlo(co=t2)
sim_carlo(co=t3)
sim_carlo(co=t4)
sim_carlo(co=t5)



# -------- 计算不同函数对应的 均值E(G) 和 方差V(t,s) -------- # 
n = 100
k = 100
## 生成n*n 阶矩阵； 循环k次，模拟检验k次,最终取k次模拟检验的均值； g是所选择的函数f(x)
c(mean_carlo(f1), sim_carlo(co=co1))
ls(n=n, k=k, g=f1)

c(mean_carlo(f2), sim_carlo(co=co2))
ls(n=n, k=k, g=f2)

c(mean_carlo(f3), sim_carlo(co=co3))
ls(n=n, k=k, g=f3)

c(mean_carlo(f4), sim_carlo(co=co4))  
ls(n=n, k=k, g=f4)

c(mean_carlo(f5), sim_carlo(co=co5))
ls(n=n, k=k ,g=f5)



# 将函数f 对应的 均值E(G) 和 方差V(t,s) 的值 写成数据框的形式
# 此处使用 c(mean_carlo(f4), sim_carlo(co=co4))  计算所得的均值和方差，没有验证极限趋势
# （因为是模拟，存在误差，如何去验证收敛数值？）

mv <- data.frame(m = c(0,1,0,0.74,0.32),v=c(2,4,24,7.24,0.63))  # 分别对应x,x2,x3,exp(x)，ln(1+x^2)


# lss test 函数定义 
# 计算Gn(f)
lss_test <- function(matrix,g,num) # 返回矩阵matrix检验的p值
{ # 写成模块/数据框 mv
  # g = x :mean = 0,var = 2 ; 
  # g = x^2 :mean = 1,var = 4;
  # g = x^3 :mean = 0,var = 24 ; 
  # g = exp(x): mean = 0.74,var = 7.24 
  # g = ln(1+x^2): mean = 0.32,var = 0.63 
  # usage ： num 对应选择1，2，3，4, 5
  Ff <- function(x) { g(x)/2/pi*sqrt(4-x^2) }
  Gn.F = integrate(Ff,-2,2,subdivisions = 10000)$value
  Gn.F = Gn.F * n    # n为矩阵的阶
  
  a = eigen(mat)$values 
  a = sum(g(a))
  Gn.f = a - Gn.F  # 此处计算Gn(f)
  
  m = mv$m[num]    # 根据选择的函数，选择mv数据框中概函数对应的值
  v = mv$v[num]
  # m = mean_carlo(g) # con = 1 -- x ; con = 2 -- x^2 ; con = 3 --x^3
  # v = sim_carlo(g=cox)
  lambdan = (Gn.f - m)/sqrt(v)
  p = pnorm(lambdan) # 积分，求取累计密度函数值
  
  if((Gn.f-m) > 0)
    p = 2*(1-p)
  else
    p = 2*p
  return(p)
}


# lss test 函数使用
## 生成随机矩阵

# --------------------------------- 检验type 1 error ------------------------- #
  k = 100
  n = 50    # n 选取20，100，200，800，1000
    
  
  p <- c()
  for(i in 1:k)
  {
    mat <- matrix(rnorm(n^2,mean=0, sd=1),n,n)
    mat[lower.tri(mat)] = 0
    mat <- mat+t(mat)-2*diag(diag(mat))+diag(rnorm(n,mean=0, sd=sqrt(2)))
    mat <- mat/sqrt(n)  # 矩阵标准化
    
    # lss_test(mat,g = f1,cox = co1)
    #p <- lss_test(mat,g = f1,cox = 1)
    #if(p < 0.05)
    #  num = num+1
    p <- append(p,lss_test(mat,g = f5,num = 5))
  }
   print(p)
  
  vol = 0
  for(j in 1:k)
  {
    if(p[j]<0.05)
      vol = vol+1
    else
      vol = vol
  }
print(vol)


# ------------------------------- 检验 type2 error ---------------------------- #
#
f1.2 <- function(x) { sin(x) }
f2.2 <- function(x) { x^2 -1 }
f3.2 <- function(x) { x^3 +1 }
f4.2 <- function(x) { exp(x) - 0.74 }

c(mean_carlo(f1.2), sim_carlo(g=co1))
ls(n=n, k=k, g=f1.2)

c(mean_carlo(f2.2), sim_carlo(g=co2))
ls(n=n, k=k, g=f2.2)

c(mean_carlo(f3.2), sim_carlo(g=co3))
ls(n=n, k=k, g=f3.2)


c(mean_carlo(f4.2), sim_carlo(g=co4))
ls(n=n, k=k, g=f4.2)

mv.2 <- data.frame(m = c(0,1,0,0.74),v=c(2,4,24,7.24))

k = 1000
n = 1000   # n 选取20，100，200，800，1000

p <- c()
for(i in 1:k)
{
  matrix(rep(1:50),nrow = n , ncol = n ,byrow =T)
  mat <- mat/sqrt(n)  # 矩阵标准化
  # lss_test(mat,g = f1,cox = co1)
  #p <- lss_test(mat,g = f1,cox = 1)
  
  #if(p < 0.05)
  #  num = num+1
  p <- append(p,lss_test(mat,g = f5,cox = 5))
}
# print(p)

num = 0
for(j in 1:k)
{
  if(p[j]>=0.05)
    num = num+1
  else
    num = num
}
print(num)
