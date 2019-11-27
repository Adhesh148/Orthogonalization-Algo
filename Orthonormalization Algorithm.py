#!/usr/bin/env python
# coding: utf-8

# In[610]:


import numpy as np
import pandas as pd
import sys
SIZE = 183


# ### Reading fs_183_6 Matrix

# In[611]:


df = pd.read_table('fs_183_6.mtx')
df.values[1:]


# In[612]:


a = np.zeros((SIZE,SIZE))
a


# In[613]:


c = np.array(df.values[1:])
for j in range(c.shape[0]):
    for i in c[j]:
        k = i.strip().split()
        a[int(k[0])-1][int(k[1])-1] = float(k[2])


# In[614]:


a = a.transpose()
n = a.shape[0]


# ### Matrix Conditions:

# In[615]:


tr = [1,2,3]
tr[0:0]


# In[616]:


a_norm = np.linalg.norm(a,ord=2)
a_cond = a_norm/np.linalg.norm(a,ord=-2)
print("norm(A): ",a_norm)
print("condition_No(A): ",a_cond)
print("float_pt_precision: ",np.finfo(float).eps)


# ### Initalising Values

# In[617]:


qt = np.zeros((n,n)) ## Q - orthogonal matrix is initialised as zero matrix
y =[] ## list of loss of orthogonality in CGS


# ## Classical Graham Schmidt Algorithm

# In[618]:


for i in range(0,n):
    w = a[i]
    for j in range(0,i):
        w = w - np.dot(np.dot(a[i],qt[j]),qt[j])
    qt[i] = w/(np.linalg.norm(w,ord=2))
    I = np.identity(i+1)
    Qj = qt[0:i+1]
    r = np.dot(Qj,np.transpose(Qj)) ## calculating transpose(Q).Q at each iteration i
    y.append(np.linalg.norm(I-r)) ## calculating loss of orthogonality at each iteration i


# In[619]:


y ## printing value of L.O for CGS


# In[620]:


#Printing Orthogonalized matrix Q
q = np.transpose(qt)
q


# ## Modified Graham Schmidt Algorithm 

# In[621]:


qt = np.zeros((n,n))
z = [] ## list of loss of orthogonality in MCGS/


# In[622]:


for i in range(0,n):
    w = a[i]
    for j in range(0,i):
        proj = np.dot(w,qt[j]) ## In MCGS a[i] is replaced with w
        w = w - proj*(qt[j])
    qt[i] = w/(np.linalg.norm(w))
    I = np.identity(i+1)
    r = np.dot(qt[0:i+1],np.transpose(qt[0:i+1]))
    z.append(np.linalg.norm(I-r))  ## Loss of orthogonality at each iteration i


# In[623]:


z


# In[624]:


#Printing Orthogonalized matrix Q
q = np.transpose(qt)
q


# ## Classic CGS with Reorthogonalization

# In[625]:


qt = np.zeros((n,n))
w = np.zeros((3,183))
u = [] ## list of loss of orthogonality in CGS with Reortho


# In[626]:


for j in range(0,n):
    w[0] = a[j]
    for i in range(1,3):
        w[i] = w[i-1]
        for k in range(0,j):
            w[i] = w[i]-np.dot(np.dot(w[i-1],qt[k]),qt[k])
    qt[j] =w[i]/np.linalg.norm(w[i],ord=2)
    I = np.identity(j+1)
    r = np.dot(qt[0:j+1],np.transpose(qt[0:j+1]))
    u.append(np.linalg.norm(I-r))  ## Loss of orthogonality at each iteration i


# In[627]:


u


# In[628]:


#Printing Orthogonalized matrix Q
q = np.transpose(qt)
q


# ## MCGS with Reorthogonalization

# In[629]:


qt = np.zeros((n,n))
w = np.zeros((3,183))
v = [] ## list of loss of orthogonality in CGS with Reortho


# In[630]:


for j in range(0,n):
    w[0] = a[j]
    for i in range(1,3):
        w[i] = w[i-1]
        for k in range(0,j):
            w[i] = w[i]-np.dot(np.dot(w[i],qt[k]),qt[k])
    qt[j] =w[i]/np.linalg.norm(w[i],ord=2)
    I = np.identity(j+1)
    r = np.dot(qt[0:j+1],np.transpose(qt[0:j+1]))
    v.append(np.linalg.norm(I-r))  ## Loss of orthogonality at each iteration i


# In[631]:


v


# In[632]:


#Printing Orthogonalized matrix Q
q = np.transpose(qt)
q


# ### Setting x-axis values

# In[633]:


x = np.arange(1,184)


# ## Ploting Loss of Orthogonality

# In[634]:


import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.yscale('log') ## y is scaled to log axis for visualization
plt.plot(x,y,'k--',label = 'CGS') ## CGS
plt.plot(x,z,'b:',label = 'MCGS') ## MGS
plt.plot(x,u,'g-',label = 'CGS with reortho') ##CGS with reortho
plt.plot(x,v,'r-.',label = 'MCGS with reortho') ## MCGS with reortho
plt.ylabel('Loss of Orthogonality')
plt.xlabel('Iteration Step')
legend = plt.legend(loc='upper left', shadow=True)
legend.get_frame().set_facecolor('C0')
plt.show


# #### The loss of orthogonality in the QR factorization for different Gram-Schmidt orthogonalization variants: 
#     CGS algorithm (dashed line) (black)
#     MGS algorithm (dotted line) (blue)
#     CGS algorithm with reorthogonalization (solid line) (green)
#     MGS algorithm with reorthogonalization (dotted-solid line) (red)
