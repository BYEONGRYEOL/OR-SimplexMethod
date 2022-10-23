# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as lin


M= np.loadtxt("C:/Users/son.d.w/example.txt", delimiter= ",", dtype= float)
c = M[0, :-1]
c = c.reshape(-1,1)
A = M[1:,:-1]
b= M[1:,-1]
b = b.reshape(-1,1)

numrow, numcol = np.shape(A)

#example 에서 row 3개

# Basis와 Non-basis 변수의 결정, Simplex method에서 계산에 활용될 행렬들 계산
def initialBasis(A):
    B = A[:,numcol-numrow:]
    N = A[:,:numcol-numrow]
    C_B = c[numcol-numrow:]
    C_B = C_B.reshape(-1,1)
    C_N = c[:numcol-numrow]
    C_N = C_N.reshape(-1,1)
    B_index = []
    for q in range(numcol-numrow+1,numcol+1):
        B_index.append(q)
    N_index = []
    for j in range(1,numcol-numrow+1):
        N_index.append(j)
    return B,N,C_B,C_N,B_index,N_index
    
def RelativeCost(C_N,C_B,B,N):
    R = C_N.T - C_B.T@lin.inv(B)@N
    return R
# 기존 Critical Point를 계산하는데에 사용되지 않은 Non-Basis변수들 중 Basis로 편입될 변수 선택
def SelectEntry(R):
    for i in range(0,numcol-numrow):
        if R[0,i] < 0:
            return i    
    if R[0,i] >= 0 :
        return -1
# 기존 Critical Point를 계산하는데에 사용된 Basis 변수들 중 Non-Basis 변수가 될(탈락할) 변수 선택
def SelectExit(i,N,B):
    N_i = N[:,i]
    B_1b = lin.inv(B)@b
    B_1N = lin.inv(B)@N_i
    B_1b1 = B_1b > 0
    B_1N1 = B_1N > 0
    j = 0
    ratio = [-333333]*(numrow)
    for j in range(0,numrow):
        if B_1b1[j] and B_1N1[j] == True:
            ratio[j] = float(B_1b[j] / B_1N[j] )
    ratio_1 = np.array(ratio)
    ratio_2 = ratio_1[ratio_1>0]
    t = ratio.index(min(ratio_2))
    
    return t
    
# 현재 iteration의 Basis의 해당하는 변수의 값들을 이용하여,
# 현재 도출할수 있는 해, Y값, 어떤 변수가 얼마나 해에 영향을 끼치는지 파악할 수 있는 sensitivity 값 들을 
# 계산하는 함수
def Basis(i,t,B,N,c,B_index,N_index):
    B_enter = B_index[t]
    N_enter = N_index[i]
    del B_index[t]
    del N_index[i]
    N_index.append(B_enter)
    B_index.append(N_enter)
    list.sort(N_index)
    list.sort(B_index)  
    B = np.zeros((numrow,numrow))
    N = np.zeros((numrow,numcol-numrow))
    for z in range(0, numrow):  
        B[:,[z]] = A[:,[B_index[z]-1]]
    for k in range(0, numcol-numrow):
        N[:,[k]] = A[:,[N_index[k]-1]] 
    cBlist = [None]*numrow
    for l in range(0,numrow):
        cBlist[l] = c[B_index[l]-1,0]
    C_B = np.array(cBlist)
    C_B = C_B.reshape(-1,1)
    cNlist = [None]*(numcol-numrow)
    for x in range(0, numcol-numrow):
        cNlist[x] = c[N_index[x]-1,0] 
    C_N = np.array(cNlist)
    C_N = C_N.reshape(-1,1)
    return B,N,C_B,C_N,B_index,N_index
    
def Result(B,C_B,b,B_index):
    X = (lin.inv(B)@b).T
    Z = C_B.T@lin.inv(B)@b
    Y = C_B.T@lin.inv(B)
    print("basis의 값들은 " , X) 
    print("basis의 종류는 ", B_index)
    print("Z의 값은 " , Z)
    print("y,sensitivity의 값들은", Y)
    

# 최적해가 결정될 수 없는 문제이거나, 같은 최적해를 도출하는 Critical Point가 유일해가 아닌 무한 해(선형)를 가질 경우를
# 검사하는 것을 구현하기에 어려움을 겪어 최대 iteration 횟수를 매번 다르게 직접 설정하도록 하였다.
R =[]
u = int(input("iteration 수를 정해주세요: "))
print("")
for u in range(1, u+1):
    print("iteration = %d " %u,)
    if u == 1:
        B,N,C_B,C_N,B_index,N_index = initialBasis(A)
        R = RelativeCost(C_N,C_B,B,N)
        i = SelectEntry(R)
        if SelectEntry(R) == -1:
            print("모든 nonbasis들의 relativecost가 양수여서 iteration을 stop한다.")
            print("iteration %d 의 값이 Optimal solution" %(u))
            break
        t = SelectExit(i,N,B)
        if SelectExit(i,N,B) == -1:
            print("basis의 값이 음수가 되서 iteration을 stop 한다.")
            print("iteration %d 의 값이 Optimal solution" %(u-1))
            break
        Result(B,C_B,b,B_index)
        print("")
    else:
        B,N,C_B,C_N,B_index,N_index = Basis(i,t,B,N,c,B_index,N_index)
        R = RelativeCost(C_N,C_B,B,N)
        i = SelectEntry(R)
        if SelectEntry(R) == -1:
            print("iteration %d  의 값이 Optimal solution" %(u))
            print("모든 nonbasis들의 relativecost가 양수여서 iteration을 stop한다.")
            break
        t = SelectExit(i,N,B)
        if SelectExit(i,N,B) == -1:
            print("basis의 값이 음수가 되서 iteration을 stop 한다.")
            print("iteration %d 의 값이 Optimal solution" %(u-1))
            break
        Result(B,C_B,b,B_index)
        print( " " )
Result(B,C_B,b,B_index)    


