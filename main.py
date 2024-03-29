import numpy 
import math
import matplotlib.pyplot as plt
#输入材料参数
#auther:郭相隆
def Get_value():
    print('imput the value of E1')
    E1 = float(input())
    print('imput the value of E2')
    E2 = float(input())
    print('imput the value of nu12')
    nu12 = float(input())
    print('imput the value of nu21')
    nu21 = float(input())
    print('imput the value of G12')
    G12 = float(input())
    return E1,E2,nu12,nu21,G12
#输入层厚度，角度
def Get_dif_layer():
    print('imput the number of layers')
    n = int(input())
    h = []
    theta = []
    for i in range(n):
        print('imput the thickness of layer',i+1)
        h.append(float(input()))
        print('imput the angle of layer',i+1)
        theta.append(float(input()))
    return h,theta
def Get_same_layer():
    print('imput the number of layers')
    n = int(input())
    theta = []
    print('imput the thickness of layer')
    h=float(input())
    for i in range(n):
        print('imput the angle of layer',i+1)
        theta.append(float(input()))
    return h,theta
#计算每层Q矩阵
def Q_matrix(E1,E2,nu12,nu21,G12,theta):
    Q = numpy.zeros((3,3))
    #theta = math.radians(theta)
    Q[0,0] = E1/(1-nu12*nu21)
    Q[0,1] = E2*nu12/(1-nu12*nu21)
    Q[0,2] = 0
    Q[1,0] = E2*nu21/(1-nu12*nu21)
    Q[1,1] = E2/(1-nu12*nu21)
    Q[1,2] = 0
    Q[2,0] = 0
    Q[2,1] = 0
    Q[2,2] = G12
    T = numpy.zeros((3,3))
    T[0,0] = math.cos(theta)**2
    T[0,1] = math.sin(theta)**2
    T[0,2] = 2*math.cos(theta)*math.sin(theta)
    T[1,0] = math.sin(theta)**2
    T[1,1] = math.cos(theta)**2
    T[1,2] = -2*math.cos(theta)*math.sin(theta)
    T[2,0] = -math.cos(theta)*math.sin(theta)
    T[2,1] = math.cos(theta)*math.sin(theta)
    T[2,2] = math.cos(theta)**2-math.sin(theta)**2
    Q_bar = numpy.dot(numpy.dot(numpy.linalg.inv(T),Q),numpy.linalg.inv(T.T))
    return Q_bar
def Q_matrix_1():
    Q = numpy.zeros((3,3))
    #theta = math.radians(theta)
    print('imput the value of Q11')
    Q[0,0] = float(input())
    print('imput the value of Q12')
    Q[0,1] = float(input())
    print('imput the value of Q22')
    Q[1,1] = float(input())
    print('imput the value of Q66')
    Q[2,2] = float(input())
    return Q
def Q_matrix_2(Q,theta):
    T = numpy.zeros((3,3))
    T[0,0] = math.cos(theta)**2
    T[0,1] = math.sin(theta)**2
    T[0,2] = 2*math.cos(theta)*math.sin(theta)
    T[1,0] = math.sin(theta)**2
    T[1,1] = math.cos(theta)**2
    T[1,2] = -2*math.cos(theta)*math.sin(theta)
    T[2,0] = -math.cos(theta)*math.sin(theta)
    T[2,1] = math.cos(theta)*math.sin(theta)
    T[2,2] = math.cos(theta)**2-math.sin(theta)**2
    Q_bar = numpy.dot(numpy.dot(numpy.linalg.inv(T),Q),numpy.linalg.inv(T.T))
    return Q_bar
#计算转变A矩阵
def A_matrix(E1,E2,nu12,nu21,G12,h,theta,Qmode=1):
    list = []
    for i in theta:
        i=math.radians(i)
        if Qmode.__class__==int:
        
            list.append(Q_matrix(E1,E2,nu12,nu21,G12,i))
        else:
            list.append(Q_matrix_2(Qmode,i))
    A = numpy.zeros((3,3))
    if h.__class__ == float:
        for i in list:
            A = A + i*h
    else:
        for i in range(len(h)):
            A = A + list[i]*h[i]
    return A
#计算B矩阵
def B_matrix(E1,E2,nu12,nu21,G12,h,theta,Qmode=1):
    list = []
    for i in theta:
        i=math.radians(i)
        if Qmode.__class__==int:
        
            list.append(Q_matrix(E1,E2,nu12,nu21,G12,i))
        else:
            list.append(Q_matrix_2(Qmode,i))
    B = numpy.zeros((3,3))
    if h.__class__ == float:
        t=h*len(theta)
        for i in range(len(list)):
            z0=-t/2+i*h
            z1=z0+h
            B = B + (z1**2-z0**2)/2*(list[i])
    else:
        t=0
        for i in h:
            t=t+i
        for k in range(len(h)):
            z0=-t/2
            z1=-t/2
            for i in range(k):
                z0=z0+h[i]
            for i in range(k+1):
                z1=z1+h[i]
            B = B + (z1**2-z0**2)/2*(list[k])
    return B
#计算D矩阵
def D_matrix(E1,E2,nu12,nu21,G12,h,theta,Qmode=1):
    list = []
    for i in theta:
        i=math.radians(i)
        if Qmode.__class__==int:
        
            list.append(Q_matrix(E1,E2,nu12,nu21,G12,i))
        else:
            list.append(Q_matrix_2(Qmode,i))
    D = numpy.zeros((3,3))
    if h.__class__ == float:
        t=h*len(theta)
        for i in range(len(list)):
            z0=-t/2+i*h
            z1=z0+h
            D = D + (z1**3-z0**3)/3*(list[i])
    else:
        t=0
        for i in h:
            t=t+i
        for k in range(len(h)):
            z0=-t/2
            z1=-t/2
            for i in range(k):
                z0=z0+h[i]
            for i in range(k+1):
                z1=z1+h[i]
            D = D + (z1**3-z0**3)/3*(list[k])
    return D

if __name__ == '__main__':
    print('使用E，Q计算','E:1,Q:2')
    choice = int(input())
    if choice == 1:
        E1,E2,nu12,nu21,G12 = Get_value()
        print('imput 1 for same layer, 2 for different layer',"输入1是相同层厚，2是不同层厚")
        choice = int(input())
        if choice == 1:
            h,theta = Get_same_layer()
        else:
            h,theta = Get_dif_layer()
        print('输入偏转角度')
        alpha = float(input())
        theta = [i+alpha for i in theta]
        A = A_matrix(E1,E2,nu12,nu21,G12,h,theta)
        numpy.around(A,decimals=4,out=A)
        print(A)
        B=B_matrix(E1,E2,nu12,nu21,G12,h,theta)
        numpy.around(B,decimals=4,out=B)
        print(B)
        D=D_matrix(E1,E2,nu12,nu21,G12,h,theta)
        numpy.around(D,decimals=4,out=D)
        print(D)
    else:
        Q=Q_matrix_1()
        print('imput 1 for same layer, 2 for different layer',"输入1是相同层厚，2是不同层厚")
        choice = int(input())
        if choice == 1:
            h,theta = Get_same_layer()
        else:
            h,theta = Get_dif_layer()
        print('输入偏转角度')
        alpha = float(input())
        theta = [i+alpha for i in theta]
        A = A_matrix(0,0,0,0,0,h,theta,Q)
        numpy.around(A,decimals=4,out=A)
        print(A)
        B=B_matrix(0,0,0,0,0,h,theta,Q)
        numpy.around(B,decimals=4,out=B)
        print(B)
        D=D_matrix(0,0,0,0,0,h,theta,Q)
        numpy.around(D,decimals=4,out=D)
        print(D)
"""[[ 5. -1.  0.]
 [-1.  5. -0.]
 [ 0. -0.  3.]]
[[ 0.1736 -0.1736  0.9848]
 [-0.1736  0.1736 -0.9848]
 [ 0.9848 -0.9848 -0.1736]]
[[ 6.6667 -1.3333  0.    ]
 [-1.3333  6.6667  0.    ]
 [ 0.     -0.      4.    ]]"""